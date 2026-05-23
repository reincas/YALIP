## Class `Levels`

This class is the interface to lanthanide states in intermediate coupling as
well as their energy levels and radiative transitions.
Generation of a `Levels` object requires to specify the electron configuration,
the coupling scheme and the set of radial integrals (in cm<sup>-1</sup>).
Specification of Judd-Ofelt parameters (in pm<sup>2</sup>) and material is
optional, but provides access to radiative dipole transition properties.

A typical initialisation code for a Pr<sup>3+</sup> ion would look like:

```
from yall import Cauchy, Coupling, Levels

config = "f2"
coupling = Coupling.SLJ
radial = {"base": 327.39, "H1/2": 68576.05, "H1/4": 49972.76, "H1/6": 32415.29, "H2": 728.18,
          "H3/0": 16.99, "H3/1": -417.98, "H3/2": 1371, "H5fix": 0.19, "H6fix": 1.67}
jo = {"JO/2": 1.981, "JO/4": 4.645, "JO/6": 6.972}
material = Cauchy(1.35123e-5, 2.94780e-3, 1.49985, -1.30933e-3, -3.23335e-6)

ion = Levels(config, coupling, radial, jo, material)
```

When the `Levels` object is initialised, it builds its total perturbation
Hamiltonian based on the given radial integrals and diagonalises it to obtain
the energy and linear combination of basis states for each state in intermediate
coupling.
The special key `"base"` fixes the energy of the ground state, all other energy
levels are shifted accordingly.

The unit of radial parameters is that of an energy, but YALL follows the standard
convention in the literature to use wavenumbers in cm<sup>-1</sup> instead.
A wavenumber is the inverse of the wavelength in vacuum

$$
k = \frac{1}{\lambda}
$$

which is proportional to the photon energy:

$$
E = h c k
$$

### Intermediate States

A `Levels` object acts as ordered list of intermediate states.
Therefore, `ion[i]` returns an `IntermediateState` object with `ion[0]` being the
ground state.
Instead of quantum numbers, an `IntermediateState` object holds a list of basis
states in its attribute `states` and the corresponding linear combination factors
in its attribute `values`.
In $SLJ$ coupling the list of basis states contains only states with the same
quantum number $J$ in contrast to all states in $SLJM$ coupling. 
The contribution weight of each base state is the square of its signed value and
the sum of the respective list attribute `weights` is therefore equal to 1.
All of these lists are ordered for decreasing weights.

The method `IntermediateState.short()` returns the short string representation of
the intermediate state's main basis state and
`IntermediateState.long(min_weight=0.0)` returns a string containing weight
factors together with the respective short string representations of all basis
states with weight factors greater or equal `min_weight`.
Furthermore the convenience method `Levels.str_levels(min_weight=0.0)` generates
a string with energy and basis state composition for each intermediate state.   

### Matrices

The methods `Levels.matrix(name)` and `Levels.reduced(name)` work exactly as they
do for a `States` object, but they return [operator matrices](operators.md) in
intermediate coupling instead.

### Radiative Transitions

If the optional arguments `jo` and `material` are provided when a `Levels` object
is initialised, it can be used to calculate the strength of radiative transitions.
The Judd-Ofelt parameters in the dictionary `jo` are expected in pm<sup>2</sup>
and the object `material` must provide a method `refractive_index(k)` which returns
the refractive index of the material for the given wavenumber in cm<sup>-1</sup>.
The YALL package provides the classes `Cauchy` and `Sellmeier` for this purpose.

There is no universally accepted definition of the line strength of a transition.
For electric dipole transitions, the YALL package uses the definition

$$
S^\mathrm{ed}_{ij} =
\frac{1}{3(2J_i+1)} \frac{e^2}{4\pi\varepsilon_0} \sum\limits_{\lambda=2,4,6}
\Omega_\lambda |\langle J_j\parallel \mathbf{U}^{(\lambda)}\parallel J_i\rangle|^2
$$

with the Judd-Ofelt parameters $\Omega_2$, $\Omega_4$, and $\Omega_6$.
For magnetic dipole transitions it uses

$$
S^\mathrm{md}_{ij} =
\frac{1}{3(2J_i+1)} \frac{1}{4\pi\varepsilon_0} \frac{\beta_m^2}{c^2 \hbar^2}
|\langle J_j\parallel \mathbf{L}+g_s\mathbf{S}\parallel J_i\rangle|^2
$$

with the Bohr magneton $\beta_m = \frac{e\hbar}{2m_e}$.
The physical unit of these line strength values is Jm<sup>3</sup>.
The matrices containing all electric and magnetic dipole transition line strengths 
are returned by the method `Levels.line_strengths()` as data class `Transition` 
with two attributes `ed` and `md`.
The following code prints the line strengths of all ground state absorptions:  

```
line = ion.line_strengths()
print(line.ed[1:, 0])
print(line.md[1:, 0])
```

The line strengths are typically used to calculate the dimensionless oscillator
strength of a transition $i\to j$ in absorption or emission:

$$
f_{ij} = \frac{4\pi\varepsilon_0}{e^2} \frac{4 \pi m_e c k_{ij}}{\hbar}
(\chi^\prime_\mathrm{ed}S^\mathrm{ed}_{ij} + \chi^\prime_\mathrm{md}S^\mathrm{md}_{ij})
$$

with the energy level difference $k_{ij}=|k_i-k_j|$ (in wavenumbers) and the local field correction factors

$$
\chi^\prime_\mathrm{ed} = \frac{(n^2+2)^2}{9n} \qquad \chi^\prime_\mathrm{md} = n
$$

The method `Levels.oscillator_strengths()` returns the oscillator strengths of
all transitions as `Transition` object and the following code prints the oscillator
strengths of all ground state absorptions with a convenient scaling factor:  

```
osc = ion.oscillator_strengths()
print(osc.ed[1:, 0] * 1e8)
print(osc.md[1:, 0] * 1e8)
```

Another useful quantity is the spontaneous radiative emission rate of a transition
in s<sup>-1</sup>:

$$
A_{ij} = \frac{32 \pi^3 k_{ij}^3}{\hbar}
(\chi_\mathrm{ed}S^\mathrm{ed}_{ij} + \chi_\mathrm{md}S^\mathrm{md}_{ij})
$$

with the local field correction factors

$$
\chi_\mathrm{ed} = \frac{n(n^2+2)^2}{9} \qquad \chi_\mathrm{md} = n^3
$$

The method `Levels.radiative_rates()` returns the radiative rates of all
transitions as `Transition` object and the following code prints the rates of all
emissions originating from the highest state:

```
A = ion.radiative_rates()
i = len(ion) - 1
print(A.ed[i-1::-1, i])
print(A.md[i-1::-1, i])
```

The radiative lifetime $\tau_i=1/\sum_jA_{ij}$ is returned by the method
`Levels.life_times()` in seconds and matrix of branching ratios
$\beta_{ij}=\tau_i A_{ij}$ by the method `Levels.branchig_ratios()`.

