## Class `Fits`

This class determines radial integrals and Judd-Ofelt parameters matching a
measured absorption spectrum.
Generation of a `Fits` object requires to specify the electron configuration,
the coupling scheme and an initial set of radial integrals (in cm<sup>-1</sup>).
Specification of the material is optional and activates the Judd_Ofelt fit.

A typical initialisation code for a Pr<sup>3+</sup> fit would look like:

```
from yalip import Cauchy, Coupling, Fits

config = "f2"
coupling = Coupling.SLJ
radial = {"base": 327.39, "H1/2": 68576.05, "H1/4": 49972.76, "H1/6": 32415.29, "H2": 728.18,
          "H3/0": 16.99, "H3/1": -417.98, "H3/2": 1371, "H5fix": 0.19, "H6fix": 1.67}
material = Cauchy(1.35123e-5, 2.94780e-3, 1.49985, -1.30933e-3, -3.23335e-6)

opt = Fits(config, coupling, radial, material)
```

### Absorption Measurement

The fitting algorithm expects line data derived from a measured absorption
spectrum in a certain list format.
The following code gives an example measurement representation for a
Pr<sup>3+</sup> ion:

```
meas = [
    [1,      '3H_5',             2365,  4, 196.2e-8,  2.8e-8],
    [2,      '3H_6',             4485, 18,  77.4e-8,  5.1e-8],
    [3,      '3F_2',             5105,  5, 283.3e-8,  6.2e-8],
    [4,      '3F_3',             6467,  7, 637.0e-8, 16.2e-8],
    [5,      '3F_4',             6958, 10, 244.2e-8, 13.1e-8],
    [6,      '1G_4',             9883,  7,  26.4e-8,  0.5e-8],
    [7,      '1D_2',            17026, 20, 199.2e-8, 10.7e-8],
    [8,      '3P_0',            20859, 12, 180.0e-8, 17.4e-8],
    [(9, 10), ('3P_1', '1I_6'), 21505, 20, 542.6e-8, 30.5e-8],
    [11,      '3P_2',           22645, 10, 926.8e-8, 28.7e-8],
]
```

Each element of this list represents an absorption line and is itself a list of
six elements:

1. Energy level index or tuple of overlapping level indices
2. Short name or tuple of short names of the respective levels
3. Measured barycenter energy of the line in cm<sup>-1</sup>
4. Error margin of the line energy in cm<sup>-1</sup>
5. Measured dimensionless oscillator strength of the line
6. Error margin of the oscillator strength

The level indices connect the measurement results to specific calculated energy
levels.
Level names currently are for information only.
They are not used by the algorithm yet, but may be verified to disentangle crossing
levels in the future.
For a level index `i` the level name should match `str(ion[i])`.

Note, that realistic error margins are quite important, because the algorithm uses
them as weight factors for the energy level fit as well as the Judd-Ofelt fit.
However, you may bypass this feature by setting the error margins of the energies
to 1 (absolute comparison) and the error margins of the oscillator strengths 
identical to the measurement values (relative comparison).

### Optimisation

The optimisation algorithm is started by calling the method
`Fits.run(lines, stages=None)` with the measured spectral data and optional
instructions for multiple stages of the energy level fit.
Without the second parameter only a Judd-Ofelt fit is performed keeping the
initial radial integrals.
The method returns the attribute `Fits.ion`, which is a
[`Levels` object](levels.md) with the optimised parameters.

In contrast to the linear Judd-Ofelt fit, the energy level fit is a nonlinear
fit which not only requires initial guesses, but also individual fitting
strategies.
Nonlinear fitting can always result in a run-off to unrealistic values or the 
algorithm may get stuck in a local minimum far from the global optimum.
Fortunately the parameter landscape of lanthanide energy levels is smooth and
fitting-friendly.
Therefore, the fast Levenberg-Marquardt algorithm used by YALIP
(`scipy.optimize.least_squares` with `method='lm'`) usually delivers reasonable
results.
Nevertheless it is advisable to perform the optimisation in several stages,
beginning with the most important first-order parameter groups `H1` (Coulomb)
and `H2` (spin-orbit) only.
In the next stage the second order Coulomb interactions `H3`and `H4`should be
added and in the final step you might check if the magnetic interactions `H5`
and `H6` can improve the result further.

For the Pr<sup>3+</sup> ion this approach would look like:

```
stages = [
    ["base", "H1/2", "H1/4", "H1/6", "H2"],
    ["base", "H1/2", "H1/4", "H1/6", "H2", "H3/0", "H3/1", ":H3/2"],
    ["base", "H1/2", "H1/4", "H1/6", "H2", "H3/0", "H3/1", ":H3/2", "H5fix", "H6fix"],
]
ion = opt.run(lines, stages)
```

The colon in `":H3/2"` tells the optimiser to use this parameter for the energy
level calculation, but not optimize it.
`"H5fix"` and `"H6fix"` are often used abbreviations for the fixed relationships
`{"H5/0": 1.0, "H5/2": 0.56, "H5/4": 0.38}` and 
`{"H6/2": 1.0, "H6/4": 0.75, "H6/6": 0.50}`, respectively.
If the `Fits` object was initialised with a materials object delivering spectral
refractive indices, the last stage of the energy-level fit will be followed by a
Judd-Ofelt fit.

In order to monitor the fitting progress in more detail, you can also call the
method run for each individual fitting stage separately: 

```
for stage in stages:
    opt.run(lines, stage)
```

### Fitting Algorithm

There is no universal solution to nonlinear fitting problems.
YALIP currently implements a basic, generic approach, and users should be prepared to make specific modifications
depending on their requirements.

The energy level fit in YALIP uses the energies $k^\mathrm{meas}_i$ of a set of measured absorption lines and their
associated error margins $\Delta k^\mathrm{meas}_i$ as references.
The initial set of radial integrals is used to construct the total perturbation Hamiltonian as a linear combination of
all individual matrices.
The eigenvalues of this total Hamiltonian yield the calculated energies of all energy levels.
Mapping these energy levels to their corresponding measured absorption lines results in a list of calculated energies
$k^\mathrm{calc}_i$, which the algorithm attempts to match to the measured energies as closely as possible.

Some absorption lines represent multiple overlapping levels.
In these cases, the calculated energy $k^\mathrm{calc}_i$ must be a weighted average of the constituent levels.
While the ideal weights would be the calculated absorption strengths of the respective levels, computing them
introduces significant computational overhead and is not yet implemented.
Consequently, YALIP takes the next best approach by using the level multiplicities as weight factors.

YALIP utilizes the Levenberg-Marquardt optimiser as a weighted nonlinear least-squares fitting algorithm to find the
optimum values for the radial integrals.
The error margin of a measured absorption line defines the weight factor for that line:

$$w_i = 1 / \Delta k^\mathrm{meas}_i$$

Because the calculated energies $k^\mathrm{calc}_i$ lack a fixed reference point, an offset adjustment is required.
However, instead of treating `base` (the ground level energy) as an explicit fit parameter, YALIP exploits the fact
that the optimum shift parameter can always be determined analytically.
The exact offset parameter is given by:

$$ k^\mathrm{off} = \frac{\sum_i w_i^2 (k^\mathrm{meas}_i - k^\mathrm{calc}_i)}{\sum_i w_i^2}$$ 

The residuals for the least-squares fit are then defined as:

$$r_i = w_i (k^\mathrm{calc}_i - k^\mathrm{meas}_i + k^\mathrm{off})$$

By definition, this choice of $k^\mathrm{off}$ ensures that the sum of the squared residuals is minimized with respect
to the offset:

$$\frac{\partial}{\partial k^\mathrm{off}} \sum_i r_i^2 = 0$$

Handling the offset parameter implicitly in this manner performs well because the Jacobian matrix is calculated
numerically.
Note, however, that an analytical Jacobian would require this parameter to be taken into account explicitly.

Finally, rather than returning the raw parameter $k^\mathrm{off}$, the fitting function automatically adjusts the
`base` parameter so that $k^\mathrm{off}$ becomes zero upon completion.

### Visualisation

For a quick visualisation of the fitting results, the method `Fits.table()` generates the lines of a text table:

```
for line in opt.table():
    print(line)
```

and the result looks like this:

```
         | kmeas        kcalc       | fmeas            fed   fmd         |     
---------+--------------------------+------------------------------------+-----
 0       |                330       |                  0.0   0.0         | 3H_4
 1  3H_5 |  2365   (4)   2364  -0.9 | 196.2   (2.8)  173.9  14.2    -8.2 | 3H_5
 2  3H_6 |  4485  (18)   4496  11.1 |  77.4   (5.1)   75.5   0.0    -1.9 | 3H_6
 3  3F_2 |  5105   (5)   5107   2.4 | 283.3   (6.2)  284.9   0.0     1.6 | 3F_2
 4  3F_3 |  6467   (7)   6463  -4.4 | 637.0  (16.2)  655.0   0.0    18.0 | 3F_3
 5  3F_4 |  6958  (10)   6955  -3.3 | 244.2  (13.1)  393.8   0.7   150.3 | 3F_4
 6  1G_4 |  9883   (7)   9883   0.3 |  26.4   (0.5)   27.0   0.4     1.0 | 1G_4
 7  1D_2 | 17026  (20)  17023  -3.0 | 199.2  (10.7)  117.6   0.0   -81.6 | 1D_2
 8  3P_0 | 20859  (12)  20857  -2.3 | 180.0  (17.4)  267.9   0.0    87.9 | 3P_0
 9  3P_1 | 21505  (20)  21472  -4.8 | 542.6  (30.5)  273.5   0.0  -127.7 | 3P_1
10  1I_6 |   ...        21507   ... |   ...          141.5   0.0     ... | 1I_6
11  3P_2 | 22645  (10)  22646   1.2 | 926.8  (28.7)  408.5   0.0  -518.3 | 3P_2
12       |              46454       |                 24.0   0.0         | 1S_0
---------+--------------------------+------------------------------------+-----
         |        11.3          2.7 |          13.1                 11.7 |     
```

The table consists of 12 columns:

1. Level index
2. Level name taken from the list of measured lines
3. Measured barycenter energy in cm<sup>-1</sup>. Overlapping levels
   are marked by ellipses (`...`).
4. Error margin of the measured energy in cm<sup>-1</sup>
5. Calculated level energy
6. Difference between calculated and measured energy
7. Measured oscillator strength in $10^{-8}$. Overlapping levels
   are marked by ellipses (`...`).
8. Error margin of the measured oscillator strength in $10^{-8}$
9. Calculated electric dipole oscillator strength in $10^{-8}$
10. Calculated magnetic dipole oscillator strength in $10^{-8}$
11. Difference between calculated and measured oscillator strength
12. Short level name

Columns 7-11 are only generated if the `Fits` object was initialised with a 
`materials` argument.
The final line contains four elements:

1. Mean error margin of all measured energies in cm<sup>-1</sup>
2. Weighted mean deviation of measured and calculated energies in cm<sup>-1</sup>
3. Mean error margin of all measured oscillator strengths in $10^{-8}$
4. Weighted mean deviation of measured and calculated oscillator strengths
   in $10^{-8}$

The weighted mean deviation of measured and calculated energies `Fits.sigma_k` is 
given by the expression 

$$
\sigma_k = \sqrt{\frac{
  \sum_i ((k^\mathrm{meas}_i - k^\mathrm{calc}_i) / \Delta k^\mathrm{meas}_i)^2
  }{\sum_i (1 / \Delta k^\mathrm{meas}_i)^2}}
$$

and the weighted mean deviation of measured and calculated oscillator strengths
`Fits.sigma_f` is

$$
\sigma_f = \sqrt{\frac{
  \sum_i ((f^\mathrm{meas}_i - f^\mathrm{ed}_i - f^\mathrm{md}_i) / \Delta f^\mathrm{meas}_i)^2
  }{\sum_i (1 / \Delta f^\mathrm{meas}_i)^2}}
$$

It is important to carefully monitor these values, because they can indicate
overfitting problems when they fall below the mean error margins.
This is the case for the energy level fit in the example above, which optimises
a matching of 9 measured lines to the calculation using 8 free parameters.
