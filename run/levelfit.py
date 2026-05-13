##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################

import logging
import numpy as np
from scipy.optimize import least_squares

from yall import Coupling, Lanthanide, RADIAL

logger = logging.getLogger("levelfit")

KMEAS = [
    [1, '3H_5', 2364.8685161861872, 4.46507163294938],
    [2, '3H_6', 4479.27317861957, 35.198082070458717],
    [3, '3F_2', 5103.3423575162715, 9.2446256526771649],
    [4, '3F_3', 6468.3053040103487, 13.304585365784103],
    [5, '3F_4', 6958.9422407794018, 20.517575554013366],
    [6, '1G_4', 9881.422924901186, 13.808791399350545],
    [7, '1D_2', 17035.775127768313, 40.903961426884479],
    [8, '3P_0', 20876.826722338206, 733.66201116333684],
    [(9, 10), ('3P_1', '1I_6'), 21505.376344086024, 231.95236384666148],
    [11, '3P_2', 22624.434389140268, 108.09826455155158]
]

STAGE = [
    ["base", "H1/2", "H1/4", "H1/6", "H2"],
    ["base", "H1/2", "H1/4", "H1/6", "H2", "H3/0", "H3/1", ":H3/2"],
    ["base", "H1/2", "H1/4", "H1/6", "H2", "H3/0", "H3/1", ":H3/2", "H5fix", "H6fix"],
]


def format_significant(value, num):
    value = f"{value:.{num}g}"
    if 'e' in value:
        value = f"{float(value):.15f}".rstrip('0').rstrip('.')
    return value


def format_params(params, num):
    return ", ".join([f"'{k}': {format_significant(params[k], num)}" for k in sorted(params.keys())])


class LevelFit:
    def __init__(self, params, matrices, mult, data_meas):
        self.params = params
        self.matrices = matrices
        self.mult = mult
        self.data_meas = data_meas

        self.num_states = self.matrices[next(iter(matrices.keys()))].shape[0]

    def get_energies(self):
        # Compute total Hamiltonian: H = sum_i(val_i * Matrix_i)
        H = np.zeros((self.num_states, self.num_states), dtype=float)
        for name in self.params:
            if name != "base":
                H += self.params[name] * self.matrices[name]

        # Diagonalize and shift energies
        energies, transform = np.linalg.eigh(H)
        weight = np.abs(transform ** 2)
        mult = self.mult[np.argmax(weight, axis=0)]
        return energies, mult

    def compare(self):
        # Calculated energies and multiplicities of all energy levels
        energies, mult = self.get_energies()

        # Collect comparison data
        results = []
        for idx, _, k_meas, dk_meas in data_meas:
            if isinstance(idx, (list, tuple)):
                idx = np.array(idx)
                m = np.sum(mult[idx])
                k_calc = np.sum(energies[idx] * mult[idx]) / m
            else:
                m = mult[idx]
                k_calc = energies[idx]
            results.append((k_meas, k_calc, dk_meas, m))
        k_meas, k_calc, dk_meas, m = np.array(results).T

        # Equalise barycenters of measured and calculated energy levels
        k0 = np.sum((k_meas - k_calc) * m) / sum(m)
        k_calc += k0
        self.params["base"] = float(energies[0] + k0)

        # Return results
        return k_meas, k_calc, dk_meas, m

    def get_residuals(self):
        k_meas, k_calc, dk_meas, m = self.compare()
        return ((k_meas - k_calc) * m) / (sum(m) * dk_meas)

    def get_chi2(self):
        residuals = self.get_residuals()
        return sum(residuals ** 2)

    def update_params(self, names, values):
        self.params |= dict(zip(names, values))

    def fit(self, opt_names):
        """ Perform an energy level fit by optimizing radial integrals to match measured energy levels. """

        def calculate(values):
            self.update_params(opt_names, values)
            return self.get_residuals()

        logger.info(format_params(self.params, 6))
        logger.info(f"Initial chi2: {self.get_chi2():.4f}")

        # Perform the optimization
        initial = [self.params[n] for n in opt_names]
        res = least_squares(calculate, initial, method='lm')
        self.update_params(opt_names, res.x.tolist())

        logger.info(format_params(self.params, 6))
        logger.info(f"Final chi2: {self.get_chi2():.4f}")


def init_logger(file_name=None, level=logging.INFO):
    root = logging.getLogger()
    root.setLevel(level)
    log_format = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

    if file_name is not None:
        file_h = logging.FileHandler(file_name, mode="a")
        file_h.setFormatter(log_format)
        file_h.setLevel(level)
        root.addHandler(file_h)

    console_h = logging.StreamHandler()
    console_h.setFormatter(log_format)
    console_h.setLevel(level)
    root.addHandler(console_h)


if __name__ == '__main__':
    init_logger(level=logging.DEBUG)

    num = 2
    coupling = Coupling.SLJ
    radial = RADIAL["Pr3+"]
    # radial = {"base": 327.39, "H1/2": 68576.05, "H1/4": 49972.76, "H1/6": 32415.29, "H2": 728.18,
    #     "H3/0": 16.99, "H3/1": -417.98, "H3/2": 1371, "H5fix": 0.19, "H6fix": 1.67}
    data_meas = KMEAS

    with Lanthanide(num, coupling, radial) as ion:
        matrix = ion.matrix("J2", coupling)
        J = (np.sqrt(4 * np.diag(matrix) + 1) - 1) / 2
        mult = 2 * J + 1
        matrices = {name: ion.matrix(name, coupling) for name in radial.keys() if name != "base"}

        for i, names in enumerate(STAGE):
            raw_names = [n[1:] if n.startswith(":") else n for n in names]
            assert len(set(raw_names)) == len(raw_names)

            opt_names = [n for n in names if n != "base" and not n.startswith(":")]
            params = {n: radial[n] for n in raw_names}

            opt = LevelFit(params, matrices, mult, data_meas)
            opt.fit(opt_names)
