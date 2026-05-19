##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################

import logging
import numpy as np
from scipy.optimize import least_squares

from .spectrum import jo_factors
from .levels import Levels

logger = logging.getLogger("yall.fit")


def format_significant(value, num):
    value = f"{value:.{num}g}"
    if 'e' in value:
        value = f"{float(value):.15f}".rstrip('0').rstrip('.')
    return value


def format_params(params, num):
    return ", ".join([f"'{k}': {format_significant(params[k], num)}" for k in sorted(params.keys())])


def format_fixed(params, num):
    return ", ".join([f"'{k}': {params[k]:.{num}f}" for k in sorted(params.keys())])


class LevelFit:
    def __init__(self, matrices, mult, lines):
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

    def fit(self, opt_names):
        """ Perform an energy level fit by optimizing the given radial integrals to match measured energy levels. """

        def calculate(values):
            self.update_params(opt_names, values)
            return self.get_residuals()

        # Perform the optimization
        initial = [self.params[n] for n in opt_names]
        res = least_squares(calculate, initial, method='lm')
        self.update_params(opt_names, res.x.tolist())


def judd_ofelt_fit(ion, lines):

    assert isinstance(ion, Levels)
    assert isinstance(lines, list)

    factor_ed, factor_md = jo_factors(ion.mult[0], ion.energies, ion.material)
    fed = np.column_stack((ion.dipole.U2[:, 0], ion.dipole.U4[:, 0], ion.dipole.U6[:, 0])) * factor_ed[:, None]
    fmd = ion.dipole.LS[:, 0] * factor_md

    A = np.zeros((len(lines), 3), dtype=float)
    b = np.zeros(len(lines), dtype=float)
    for i, (idx, f_meas, df_meas) in enumerate(lines):
        if isinstance(idx, int):
            b[i] = (f_meas - fmd[idx]) / df_meas
            A[i,:] = fed[idx, :] / df_meas
        else:
            idx = np.array(idx)
            b[i] = (f_meas - np.sum(fmd[idx])) / df_meas
            A[i,:] = np.sum(fed[idx, :], axis=0) / df_meas

    omega, residuals, rank, _ = np.linalg.lstsq(A, b, rcond=None)
    chi2 = residuals[0]
    assert rank == 3

    df_meas = np.array([line[2] for line in lines])
    sigma = float(np.sqrt(chi2 / np.sum(1 / df_meas ** 2))) * 1e8
    judd_ofelt = {f"JO/{2*i+2}": value for i, value in enumerate(omega)}
    return judd_ofelt, sigma


def str_compare(lines, states):
    meas = len(states) * [{"type": "empty"}]
    for line in lines:
        idx, name, k_meas, dk_meas = line[:4]
        if isinstance(idx, (list, tuple)):
            meas[idx[0]] = {"type": "overlapped", "range": idx, "k": k_meas, "dk": dk_meas, "name": name[0]}
            for i, n in zip(idx[1:], name[1:]):
                meas[i] = {"type": "continue", "name": n}
        else:
            meas[idx] = {"type": "normal", "k": k_meas, "dk": dk_meas, "name": name}

    result = []
    for i, state in enumerate(states):
        name_calc = state.short()
        level_type = meas[i]["type"]
        if level_type == "normal":
            name_meas = meas[i]["name"]
            k_meas = f"{meas[i]["k"]:.0f}"
            dk_meas = f"({meas[i]["dk"]:.0f})"
            k_calc = f"{state.energy:.0f}"
            dk_calc = f"{state.energy - meas[i]["k"]:.1f}"
        elif level_type == "overlapped":
            name_meas = meas[i]["name"]
            k_meas = f"{meas[i]["k"]:.0f}"
            dk_meas = f"({meas[i]["dk"]:.0f})"
            k = np.array([states[i].energy for i in meas[i]["range"]])
            m = np.array([states[i].mult for i in meas[i]["range"]])
            k = np.sum(k * m) / np.sum(m)
            k_calc = f"{state.energy:.0f}"
            dk_calc = f"{k - meas[i]["k"]:.1f}"
        elif level_type == "continue":
            name_meas = meas[i]["name"]
            k_meas = "..."
            dk_meas = ""
            k_calc = f"{state.energy:.0f}"
            dk_calc = "..."
        else:
            name_meas = ""
            k_meas = ""
            dk_meas = ""
            k_calc = f"{state.energy:.0f}"
            dk_calc = ""
        result.append((str(i), name_meas, k_meas, dk_meas, k_calc, dk_calc, name_calc))

    width = [max(map(len, values)) for values in zip(*result)]
    fmt = "  ".join([f"{{:>{width[i]}s}}" for i in range(len(width))])
    for values in result:
        yield fmt.format(*values)


class Fit:
    def __init__(self, config, radial, material=None):

        assert isinstance(config, str)
        assert isinstance(radial, dict)

        # Electron configuration
        self.config = config

        # Intermediate coupling object
        self.ion = Levels(config, radial, None, material)

        # Material object providing spectral refractive indices
        self.material = material

        # State multiplicities
        self.mult = np.array(self.base_states.mult)

        # Perturbation energy matrices
        self.matrices = {name: self.base_states.matrix(name) for name in radial.keys() if name != "base"}

        # No fit yet
        self.levels = None

    @property
    def coupling(self):
        """ Coupling scheme. """

        return self.ion.coupling

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

    def level_fit(self, lines, stages):
        """ Multi-stage energy level fit to measured absorption lines."""

        # Measured absorption lines
        assert isinstance(lines, list)
        size = set(len(line) for line in lines)
        assert len(size) == 1
        size = size.pop()
        assert size in (4, 6)
        self.lines = lines
        k_lines = [[line[0], line[2], line[3]] for line in lines]

        # Copy values of radial integrals
        radial = self.radial_integrals.copy()

        # Handle single optimisation stage
        if not isinstance(stages[0], (list, tuple)):
            stages = [stages]

        opt = LevelFit(self.matrices, self.mult, k_lines)
        for i, names in enumerate(stages):
            raw_names = [n[1:] if n.startswith(":") else n for n in names]
            assert len(set(raw_names)) == len(raw_names)
            opt.set_params({n: radial[n] for n in raw_names})

            p = format_params(opt.params, 6)
            dk = opt.get_sigma()
            logger.info(f"Stage {i}: Initial dk: {dk:.2f}, parameters: {p}")

            names = [n for n in names if n != "base" and not n.startswith(":")]
            opt.fit(names)
            radial |= opt.params

            p = format_params(opt.params, 6)
            dk = opt.get_sigma()
            logger.info(f"Stage {i}: Final dk: {dk:.2f}, parameters: {p}")

            self.ion = Levels(self.config, opt.params, None, self.material)

            # Judd-Ofelt fit
            if size == 6:
                f_lines = [[line[0], line[4], line[5]] for line in lines]
                judd_ofelt, df = judd_ofelt_fit(self.ion, f_lines)
                self.ion.judd_ofelt = judd_ofelt
                p = format_fixed(judd_ofelt, 3)
                logger.info(f"Stage {i}: Final df: {df:.2f}, parameters: {p}")

    def str_compare(self):
        assert self.lines is not None, "Run an energy level fit first!"

        yield from str_compare(self.lines, self.ion.states)
