##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################
#
# This module is used to perform energy level and Judd-Ofelt fits to
# determine optimised radial integrals and Judd-Ofelt parameters matching
# measured absorption lines through the class `Fits`.
#
##########################################################################
import math
from dataclasses import dataclass
import logging

import numpy as np
from scipy.optimize import least_squares

from . import Coupling
from .spectrum import jo_factors, Transition
from .matrix import normalize_radial
from .levels import Levels

logger = logging.getLogger("yalip.fit")


##########################################################################
# Formatter functions
##########################################################################

def format_significant(value, num):
    """ Return floating point string representation of the given value rounded to the given number of significant
    digits."""

    value = f"{value:.{num}g}"
    if 'e' in value:
        value = f"{float(value):.15f}".rstrip('0').rstrip('.')
    return value


def format_params(params, num):
    """ Return string representation of the given parameter set rounded to the given number of significant digits. """

    return ", ".join([f"'{k}': {format_significant(params[k], num)}" for k in sorted(params.keys())])


def format_fixed(params, num):
    """ Return string representation of the given parameter set rounded to the given number of decimal places. """

    return ", ".join([f"'{k}': {params[k]:.{num}f}" for k in sorted(params.keys())])


##########################################################################
# Result table
##########################################################################

@dataclass
class MeasLine:
    """ Dataclass containing measured values. """

    level: int | tuple[int, ...] | None
    name: str
    value: float | str
    delta: float | str


class MeasLines:
    """ This class provides level-based access to measured absorption lines. """

    def __init__(self, lines, num_states):
        """ Extract and store measurement data for each energy level. """

        assert isinstance(lines, list)
        assert isinstance(num_states, int)
        self.num_states = num_states

        # Extract measured data for each energy level
        self.lines = [MeasLine(None, "", "", "") for i in range(num_states)]
        for line in lines:
            idx, name, value, delta = line[:4]
            if isinstance(idx, tuple):
                assert all(a - b == 1 for a, b in zip(idx[1:], idx[:-1]))
                self.lines[idx[0]] = MeasLine(idx, name[0], value, delta)
                for i, n in zip(idx[1:], name[1:]):
                    self.lines[i] = MeasLine(None, n, "...", "")
            else:
                self.lines[idx] = MeasLine(idx, name, value, delta)

    def __getitem__(self, i):
        """ Return measurement data of the given level as MeasLine object. """

        return self.lines[i]

    def __iter__(self):
        """ Generate measurement data of each level as MeasLine object. """

        for line in self.lines:
            yield line


class MeasBase:
    """ Base class for level-based comparison classes. """

    val_fmt: dict
    row_fmt: str
    levels: list

    def __getitem__(self, i):
        """ Return measured and calculated data of the given level as MeasEnergy object. """

        return self.levels[i]

    def __iter__(self):
        """ Generate measured and calculated data of each level as MeasEnergy object. """

        for level in self.levels:
            yield level

    def get_sigma(self, name_delta, name_diff):
        """ Return mean error margin and weighted average deviation of measured and calculated values. """

        levels = [level for level in self.levels if not isinstance(getattr(level, name_delta), str)]
        delta, diff = zip(*[(getattr(level, name_delta), getattr(level, name_diff)) for level in levels])
        weights = [1 / value for value in delta]

        mean = sum(delta) / len(delta)
        sigma = math.sqrt(sum((w * v) ** 2 for w, v in zip(weights, diff)) / sum(w ** 2 for w in weights))
        return mean, sigma

    def str_values(self, name, fmt):
        """ Returns formatted level values as strings and their maximum size. """

        values = [getattr(level, name) for level in self.levels]
        values = [value if isinstance(value, str) else format(value, fmt) for value in values]
        size = max(len(value) for value in values)
        return values, size

    def row_str(self, row, sizes):
        """ Return given row elements as string. """

        return self.row_fmt.format(*[format(v, f">{s}s") for v, s in zip(row, sizes)])

    def iter_table(self, head, foot):
        """ Generate table lines with given header and footer. """

        assert isinstance(head, list)
        assert isinstance(foot, list)

        # Data lines and columns widths
        data, sizes = zip(*[self.str_values(name, fmt) for name, fmt in self.val_fmt.items()])

        # Separation line
        sep = [size * "-" for size in sizes]
        fmt = self.row_fmt.replace("|", "+").replace(" ", "-")
        sep = fmt.format(*sep)

        # Generate text table
        yield self.row_str(head, sizes)
        yield sep
        for row in list(zip(*data)):
            yield self.row_str(list(row), sizes)
        yield sep
        yield self.row_str(foot, sizes)


@dataclass
class MeasEnergy:
    """ Dataclass containing measured and calculated wavenumbers. """

    level: int
    line_name: str
    meas: float | str
    delta: float | str
    calc: float
    diff: float | str
    name: str


class MeasEnergies(MeasBase):
    """ This class provides a level-based comparison of measured and calculated wavenumbers. """

    val_fmt = {
        "level": "d", "line_name": "s",
        "meas": ".0f", "delta": ".0f", "calc": ".0f", "diff": ".1f",
        "name": "s",
    }
    row_fmt = "{}  {} | {}  {}  {}  {} | {}"

    def __init__(self, lines, names, energies, mult):
        assert isinstance(lines, list)
        assert isinstance(names, list)
        assert isinstance(energies, list)
        assert isinstance(mult, list)
        assert len(energies) == len(mult) == len(names)

        self.names = names
        self.energies = energies
        self.mult = mult

        # Prepare line measurement objects
        self.num_states = len(names)
        lines = [line[:4] for line in lines]
        self.lines = MeasLines(lines, self.num_states)

        # Prepare measured and calculated data for all energy levels
        self.levels = []
        for i in range(self.num_states):
            meas = self.lines[i]
            calc = self.energies[i]
            if isinstance(meas.level, int):
                diff = self.energies[i] - meas.value
            elif isinstance(meas.level, tuple):
                weight = sum(self.mult[j] for j in meas.level)
                diff = sum(self.energies[j] * self.mult[j] for j in meas.level) / weight - meas.value
            else:
                diff = meas.value
            name = self.names[i]
            self.levels.append(MeasEnergy(i, meas.name, meas.value, meas.delta, calc, diff, name))

        # Mean error margin and weighted average deviation of measured and calculated values
        self.mean, self.sigma = self.get_sigma("delta", "diff")

    def table(self):
        """ Generate line strings of result table. """

        # Header elements
        head = ["", "", "kmeas", "", "kcalc", "", ""]

        # Footer elements
        sigma = format(self.sigma, self.val_fmt["diff"])
        mean = format(self.mean, self.val_fmt["delta"])
        foot = ["", "", "", mean, "", sigma, ""]

        # Generate table lines
        yield from self.iter_table(head, foot)


@dataclass
class MeasStrength:
    """ Dataclass containing measured and calculated oscillator strengths. """

    level: int
    line_name: str
    meas: float | str
    delta: float | str
    ed: float
    md: float
    diff: float | str
    name: str


class MeasStrengths(MeasBase):
    """ This class provides a level-based comparison of measured and calculated oscillator strengths. """

    val_fmt = {
        "level": "d", "line_name": "s",
        "meas": ".1f", "delta": ".1f", "ed": ".1f", "md": ".1f", "diff": ".1f",
        "name": "s",
    }
    row_fmt = "{}  {} | {}  {}  {}  {}  {} | {}"

    def __init__(self, lines, names, strengths):
        assert isinstance(lines, list)
        assert isinstance(names, list)
        assert isinstance(strengths, Transition)
        assert len(strengths.ed) == len(strengths.md) == len(names)

        self.names = names
        self.strengths = strengths

        # Prepare line measurement objects
        self.num_states = len(names)
        lines = [line[:2] + line[4:6] for line in lines]
        self.lines = MeasLines(lines, self.num_states)

        # Prepare measured and calculated data for all energy levels
        self.levels = []
        for i in range(self.num_states):
            meas = self.lines[i]
            ed = self.strengths.ed[i]
            md = self.strengths.md[i]
            if isinstance(meas.level, int):
                diff = (self.strengths.ed[i] + self.strengths.md[i]) - meas.value
            elif isinstance(meas.level, tuple):
                diff = sum(self.strengths.ed[j] + self.strengths.md[j] for j in meas.level) - meas.value
            else:
                diff = meas.value
            name = self.names[i]
            self.levels.append(MeasStrength(i, meas.name, meas.value, meas.delta, ed, md, diff, name))

        # Mean error margin and weighted average deviation of measured and calculated values
        self.mean, self.sigma = self.get_sigma("delta", "diff")

    def table(self):
        """ Generate line strings of result table. """

        # Header elements
        head = ["", "", "fmeas", "", "fed", "fmd", "", ""]

        # Footer_elements
        sigma = format(self.sigma, self.val_fmt["diff"])
        mean = format(self.mean, self.val_fmt["delta"])
        foot = ["", "", "", mean, "", "", sigma, ""]

        # Generate table lines
        yield from self.iter_table(head, foot)


@dataclass
class MeasLevel:
    """ Dataclass containing measured and calculated wavenumbers and oscillator strengths. """

    level: int
    line_name: str
    k_meas: float | str
    k_delta: float | str
    k_calc: float
    k_diff: float | str
    f_meas: float | str
    f_delta: float | str
    f_ed: float
    f_md: float
    f_diff: float | str
    name: str


class MeasLevels(MeasBase):
    """ This class provides a level-based comparison of measured and calculated wavenumbers and oscillator
    strengths. """

    val_fmt = {
        "level": "d", "line_name": "s",
        "k_meas": ".0f", "k_delta": ".0f", "k_calc": ".0f", "k_diff": ".1f",
        "f_meas": ".1f", "f_delta": ".1f", "f_ed": ".1f", "f_md": ".1f", "f_diff": ".1f",
        "name": "s",
    }
    row_fmt = "{}  {} | {}  {}  {}  {} | {}  {}  {}  {}  {} | {}"

    def __init__(self, lines, names, energies, mult, strengths):
        self.num_states = len(names)

        # Comparison objects for wavenumbers and oscillator strengths
        self.energies = MeasEnergies(lines, names, energies, mult)
        self.strengths = MeasStrengths(lines, names, strengths)

        # Prepare measured and calculated data for all energy levels
        self.levels = []
        for i in range(self.num_states):
            k_meas = self.energies[i]
            f_meas = self.strengths[i]
            k_data = (k_meas.meas, k_meas.delta, k_meas.calc, k_meas.diff)
            f_data = (f_meas.meas, f_meas.delta, f_meas.ed, f_meas.md, f_meas.diff)
            self.levels.append(MeasLevel(i, k_meas.line_name, *k_data, *f_data, f_meas.name))

        # Mean error margin and weighted average deviation of measured and calculated values
        self.k_mean, self.k_sigma = self.get_sigma("k_delta", "k_diff")
        self.f_mean, self.f_sigma = self.get_sigma("f_delta", "f_diff")

    def table(self):
        """ Generate line strings of result table. """

        # Header elements
        head = ["", "", "kmeas", "", "kcalc", "", "fmeas", "", "fed", "fmd", "", ""]

        # Footer elements
        k_sigma = format(self.k_sigma, self.val_fmt["k_diff"])
        k_mean = format(self.k_mean, self.val_fmt["k_delta"])
        f_sigma = format(self.f_sigma, self.val_fmt["f_diff"])
        f_mean = format(self.f_mean, self.val_fmt["f_delta"])
        foot = ["", "", "", k_mean, "", k_sigma, "", f_mean, "", "", f_sigma, ""]

        # Generate table lines
        yield from self.iter_table(head, foot)


##########################################################################
# Energy level fit
##########################################################################

class LevelFit:
    """ This class is used to perform an energy level fit using the Levenberg-Marquardt algorithm. """

    def __init__(self, matrices, mult, lines):
        """ Store operator matrices, state multiplicities, and measured absorption lines. """

        self.matrices = matrices
        self.mult = mult
        self.lines = lines

        self.params = {}
        self.num_states = self.matrices[next(iter(matrices.keys()))].shape[0]

    def set_params(self, params):
        """ Store a set of radial parameters. """

        assert isinstance(params, dict)
        self.params = params

    def get_energies(self):
        """ Diagonalize the total perturbation Hamiltonian and return the energies and multiplicities of all states
         in intermediate coupling. """

        # Total perturbation Hamiltonian
        H = np.zeros((self.num_states, self.num_states), dtype=float)
        for name in self.params:
            if name != "base":
                H += self.params[name] * self.matrices[name]

        # Diagonalize Hamiltonian
        energies, transform = np.linalg.eigh(H)

        # Multiplicity of each energy level in intermediate coupling
        weight = np.abs(transform ** 2)
        mult = self.mult[np.argmax(weight, axis=0)]

        # Return energies and multiplicities of all states
        return energies, mult

    def compare(self):
        """ Return measured and calculated wavenumbers, measurement accuracies, and multiplicities of all measured
        absorption lines. """

        # Calculated energies and multiplicities of all energy levels
        energies, mult = self.get_energies()

        # Collect comparison data
        results = []
        for idx, k_meas, dk_meas in self.lines:
            if isinstance(idx, int):
                m = mult[idx]
                k_calc = energies[idx]
            else:
                idx = np.array(idx)
                m = np.sum(mult[idx])
                k_calc = np.sum(energies[idx] * mult[idx]) / m
            results.append((k_meas, k_calc, dk_meas, m))
        k_meas, k_calc, dk_meas, m = np.array(results).T

        # Equalise barycenter of measured and calculated energy levels
        k0 = np.sum((k_meas - k_calc) * m) / sum(m)
        k_calc += k0
        self.params["base"] = float(energies[0] + k0)

        # Return results
        return k_meas, k_calc, dk_meas, m

    def get_residuals(self):
        """ Return residuals of all measured absorption lines weighted by their multiplicities and measurement
        accuracies. """

        k_meas, k_calc, dk_meas, m = self.compare()
        return (k_meas - k_calc) / dk_meas

    def get_chi2(self):
        """ Return square sum of residuals. """

        residuals = self.get_residuals()
        return float(np.sum(residuals ** 2))

    def get_sigma(self):
        """ Return weighted average deviation of measured and calculated absorption lines. """

        k_meas, k_calc, dk_meas, m = self.compare()
        return float(np.sqrt(np.sum(((k_meas - k_calc) / dk_meas) ** 2) / np.sum(1 / dk_meas ** 2)))

    def update_params(self, names, values):
        """ Update given radial integrals. """

        self.params = self.params | dict(zip(names, values))

    def run(self, opt_names):
        """ Perform an energy level fit by optimizing the given radial integrals to match measured energy levels. """

        def calculate(values):
            self.update_params(opt_names, values)
            return self.get_residuals()

        # Perform the optimization
        initial = [self.params[n] for n in opt_names]
        res = least_squares(calculate, initial, method='lm')
        self.update_params(opt_names, res.x.tolist())


##########################################################################
# Judd-Ofelt fit
##########################################################################

def judd_ofelt_fit(ion, lines):
    """ Perform a Judd-Ofelt fit using a linear least squares operation. """

    assert isinstance(ion, Levels)
    assert isinstance(lines, list)

    # Prepare oscillator strength factors
    factor_ed, factor_md = jo_factors(ion.mult[0], ion.energies, ion.material)
    fed = np.column_stack((ion.dipole.U2[:, 0], ion.dipole.U4[:, 0], ion.dipole.U6[:, 0])) * factor_ed[:, None]
    fmd = ion.dipole.LS[:, 0] * factor_md

    # Build matrix A of calculated values and result vector b of measured values
    A = np.zeros((len(lines), 3), dtype=float)
    b = np.zeros(len(lines), dtype=float)
    for i, (idx, f_meas, df_meas) in enumerate(lines):
        if isinstance(idx, int):
            b[i] = (f_meas - fmd[idx]) / df_meas
            A[i, :] = fed[idx, :] / df_meas
        else:
            idx = np.array(idx)
            b[i] = (f_meas - np.sum(fmd[idx])) / df_meas
            A[i, :] = np.sum(fed[idx, :], axis=0) / df_meas

    # Perform linear least squares fit resulting in the Judd-Ofelt parameters omega
    omega, residuals, rank, _ = np.linalg.lstsq(A, b, rcond=None)
    chi2 = residuals[0]
    assert rank == 3

    # Return parameter dictionary and weighted mean deviation of measured and calculated oscillator strengths
    df_meas = np.array([line[2] for line in lines])
    sigma = float(np.sqrt(chi2 / np.sum(1 / df_meas ** 2))) * 1e8
    judd_ofelt = {f"JO/{2 * i + 2}": value for i, value in enumerate(omega)}
    return judd_ofelt, sigma


##########################################################################
# Class Fits
##########################################################################

class Fits:
    """ This class is used to perform energy level and Judd-Ofelt fits to determine optimised radial integrals and
    Judd-Ofelt parameters matching measured absorption lines. """

    def __init__(self, config, coupling, radial, material=None):

        assert isinstance(config, str)
        assert isinstance(coupling, Coupling)
        assert isinstance(radial, dict)

        # Electron configuration
        self.config = config

        # Coupling scheme
        self.coupling = coupling

        # Radial integrals
        self.radial = normalize_radial(radial)

        # Intermediate coupling object
        self.ion = Levels(config, coupling, self.radial, None, material)

        # Material object providing spectral refractive indices
        self.material = material
        self.has_strengths = material is not None

        # State multiplicities
        self.mult = np.array(self.base_states.mult)

        # Perturbation energy matrices
        self.matrices = {name: self.base_states.matrix(name) for name in self.radial.keys() if name != "base"}

        # No fit yet
        self.levels = None
        self.sigma_k = None
        self.sigma_f = None

    @property
    def base_states(self):
        """ Basis states. """

        return self.ion.base_states

    @property
    def radial_integrals(self):
        """ Radial integrals. """

        return self.ion.radial_integrals

    @property
    def judd_ofelt(self):
        """ Judd-Ofelt parameters. """

        return self.ion.judd_ofelt

    def run(self, lines, stages=None):
        """ Multi-stage combined energy level and Judd_ofelt fit to measured absorption lines. """

        # Measured absorption lines
        assert isinstance(lines, list)
        size = set(len(line) for line in lines)
        assert len(size) == 1
        size = size.pop()
        assert size == 6
        self.lines = lines

        if stages is None:
            self.ion = Levels(self.config, self.coupling, self.radial, None, self.material)

        else:
            # Handle single optimisation stage
            if not isinstance(stages[0], (list, tuple)):
                stages = [stages]

            k_lines = [[line[0], line[2], line[3]] for line in lines]
            opt = LevelFit(self.matrices, self.mult, k_lines)
            for i, names in enumerate(stages):
                raw_names = [n[1:] if n.startswith(":") else n for n in names]
                assert len(set(raw_names)) == len(raw_names)
                opt.set_params({n: self.radial[n] for n in raw_names})

                p = format_params(opt.params, 6)
                dk = opt.get_sigma()
                logger.info(f"Stage {i}: Initial dk: {dk:.2f}, parameters: {p}")

                names = [n for n in names if n != "base" and not n.startswith(":")]
                opt.run(names)
                self.radial |= opt.params

                p = format_params(opt.params, 6)
                self.sigma_k = opt.get_sigma()
                logger.info(f"Stage {i}: Final dk: {self.sigma_k:.2f}, parameters: {p}")

                self.ion = Levels(self.config, self.coupling, opt.params, None, self.material)

        # Judd-Ofelt fit
        if self.has_strengths:
            f_lines = [[line[0], line[4], line[5]] for line in lines]
            judd_ofelt, self.sigma_f = judd_ofelt_fit(self.ion, f_lines)
            self.ion.judd_ofelt = judd_ofelt
            p = format_fixed(judd_ofelt, 3)
            logger.info(f"Judd-Ofelt fit: df: {self.sigma_f:.2f}, parameters: {p}")

        # Return optimised Levels object
        return self.ion

    def table(self):
        """ Generate text table comparing measured and calculated absorption lines. """

        assert self.lines is not None, "Run an energy level fit first!"

        # Scale oscillator strengths
        lines = [line.copy() for line in self.lines]
        if self.has_strengths:
            strengths = self.ion.oscillator_strengths()
            strengths.ed = strengths.ed[:, 0] * 1e8
            strengths.md = strengths.md[:, 0] * 1e8
            for i in range(len(lines)):
                lines[i][4] *= 1e8
                lines[i][5] *= 1e8
        else:
            strengths = None

        # Generate line strings of table
        names = [state.short() for state in self.ion.states]
        energies = list(self.ion.states.energies)
        mult = list(self.ion.states.mult)

        # Select table type
        if len(set(line[2] for line in lines)) == 1:
            meas = MeasStrengths(lines, names, strengths)
        elif self.has_strengths:
            meas = MeasLevels(lines, names, energies, mult, strengths)
        else:
            meas = MeasEnergies(lines, names, energies, mult)

        # Generate table rows
        for line in meas.table():
            yield line
