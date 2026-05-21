# YALL 0.1.0 (Yet Another Lanthanide Library)

This is a Python 3 package to calculate the energy levels of multi-electron
systems populating the 4f configuration, which means the lanthanide or
rare-earth ions from Ce<sup>3+</sup> (4f<sup>1</sup>) to Yb<sup>3+</sup>
(4f<sup>13</sup>).
The calculation is based on Racah's tensor algebra using an approach
introduced in my PhD thesis [[1]](#ref1) and advanced in [[2]](#ref2).
You find a copy of the [thesis](docs/Dissertation.pdf) and
[corrections](docs/errata5.pdf) in the folder `docs`.
The original version of the software is available in the GitHub repository
[Lanthanide-0.3](https://github.com/reincas/Lanthanide-0.3) and the last
implementation of the original approach in the repository
[Lanthanide](https://github.com/reincas/Lanthanide).
However, the YALL package is using the matrix elements from the Zenodo repository
[AMELI](https://zenodo.org/communities/ameli) calculated and stored in exact
arithmetic according to [[2]](#ref2).

## Matrices of Tensor Operators

Three groups of spherical tensor operators are supported by the AMELI
repository: Unit and Coulomb tensor operators, angular momentum operators, and
perturbation Hamilton operators.
The matrix elements of each operator are provided in the product state space,
for full $SLJM$ coupling, stretched states in $SLJ$ coupling and reduced matrix
elements in the $SLJ$ space.

All matrix elements have actually been calculated in the space of determinantal
product states and then transformed to the $SLJM$ space exploiting the symmetry
properties of Racah's chain of tensor operators:

$$ \mathbf{S}^2 \to \mathbf{C}_2(\mathrm{SO}(7)) \to \mathbf{C}_2(\mathrm{G}_2)
\to \mathbf{L}^2 \to \mathbf{J}^2 \to \mathrm{J}_z $$

with the operators of total spin $\mathbf{S}$, total orbital angular momentum
$\mathbf{L}$, total angular momentum $\mathbf{J}$ and it's magnetic component
$\mathbf{J}_z$.
Furthermore, there is Casimir's quadratic operator
$\mathbf{C}_2$ of the special orthogonal group $\mathrm{SO}(7)$ in seven
dimensions $(2l+1)$ and the special group $\mathrm{G}_2$.

## Perturbation Hamiltonians

The centerpiece of AMELI and YALL are six free-ion perturbation Hamiltonians
resembling a total of 19 radial integrals.
They are used to calculate the energy levels of a given Lathanide ion.
In the order of their relative magnitude they are:

1. Coulomb interaction between the electrons (1st order):
   $\mathbf{H}_1 = F^2 \mathbf{f}_2 + F^4 \mathbf{f}_4 + F^6 \mathbf{f}_6$,
   with the radial integrals $F^2$, $F^4$, and $F^6$
   and the respective angular two-electron operators $\mathbf{f}_2$,
   $\mathbf{f}_4$, and $\mathbf{f}_6$.
   The YALL package uses the keys `"H1/2"`, `"H1/4"`, and `"H1/6"` for radial 
   parameters (integrals) and angular matrices (operators).
2. Magnetic spin-orbit interaction of each electron (1st order):
   $\mathbf{H}_2 = \zeta \mathbf{z}$,
   with the radial integral $\zeta$ and the angular one-electron operator
   $\mathbf{z}$.
   The YALL package uses the key `"H2"` for the radial parameter and the angular
   matrix.
3. Coulomb inter-configuration interactions (2nd order):
   $\mathbf{H}_3 = \alpha \mathbf{L}^2 + \beta \mathbf{C}_2(\mathrm{G}_2) + \gamma \mathbf{C}_2(\mathrm{SO}(7))$,
   with the radial integrals $\alpha$, $\beta$, and $\gamma$.
   The respective effective angular two-electron operators $\mathbf{L}^2$,
   $\mathbf{C}_2(\mathrm{G}_2)$, and $\mathbf{C}_2(\mathrm{SO}(7))$ are the squared operator of the total
   orbital angular momentum and the mentioned quadratic Casimir operators.
   The YALL package uses the keys `"H3/0"`, `"H3/1"`, and `"H3/2"` for radial
   parameters and angular matrices.
4. More Coulomb inter-configuration interactions (2nd order):
   $\mathbf{H}_4 = \sum_c T^c \mathbf{t}_c$, with $c = 2, 3, 4, 6, 7, 8$.
   The radial integrals are $T^c$ and the respective effective angular three-electron operators
   $\mathbf{t}_c$.
   The YALL package uses the keys `"H4/2"`, `"H4/3"`, `"H4/4"`, `"H4/6"`, `"H4/7"`, and `"H4/8"` for
   radial parameters and angular matrices.
5. Magnetic interactions of the spin of one electron with the spin (ss) or the orbital angular momentum (soo)
   of another electron (1st order):
   $\mathbf{H}_5 = M^0 \mathbf{m}_0 + M^2 \mathbf{m}_2 + M^4 \mathbf{m}_4$,
   with the Marvin integrals $M^0$, $M^2$, and $M^4$ and the respective angular two-electron operators
   $\mathbf{m}_0$, $\mathbf{m}_2$, and $\mathbf{m}_4$.
   The YALL package uses the keys `"H5/0"`, `"H5/2"`, and `"H5/4"` for radial parameters and angular matrices.
6. Magnetic inter-configuration spin-orbit interactions (2nd order):
   $\mathbf{H}_6 = P^2 \mathbf{p}_2 + P^4 \mathbf{p}_4 + P^6 \mathbf{p}_6$,
   with the radial integrals $P^2$, $P^4$, and $P^6$ and the respective angular two-electron operators
   $\mathbf{p}_2$, $\mathbf{p}_4$, and $\mathbf{p}_6$.
   The YALL package uses the keys `"H6/2"`, `"H6/4"`, and `"H6/6"` for radial parameters and angular matrices.

The first order perturbations are interactions inside the 4f configuration, while
second order perturbations take interactions with all other configurations into
account.
They are treated by effective operators mathematically operating inside the 4f
configuration.

Energy level calculations for lanthanides in crystalline hosts beyond the free-ion
approximation are supported by crystal field Hamiltonians `"Hcf/{k},{q}"` with the
radial integrals $B^k_q$.
Using these not only results in Stark splitting, because breaking the spherical
symmetry of the free ion also mixes states from different $J$-multiplets, which
means that $J$ is no longer a "good" quantum number.

## Installation

The package is available on PyPI and the installation therefore is possible using pip

```
python -m pip install yall
```

## Usage

The YALL package is designed for three use cases on different abstraction levels:

1. Access to raw numerical representations of states and matrices from the AMELI
repository using the class `States`
2. Calculation of intermediate states, energy levels and radiative transitions
based on given radial integrals and Judd-Ofelt parameters using the class `Levels`
3. Performing energy level and Judd-Ofelt fits to find optimised radial integrals
and Judd-Ofelt parameters matching a measured absorption spectrum using the class
`Fit`

## Class `States`

This class provides a convenient interface the all states and matrix elements
from the AMELI repository.
The following code initialises all states of the Pr<sup>3+</sup> ion in $SLJ$
coupling:

```
from yall import States, Coupling
config = "f2"
coupling = Coupling.SLJ
states = States(config, coupling)
```

`Coupling.SLJ` is the natural choice for lanthanide ions in amorphous
materials, `Coupling.SLJM` would be suitable for crystalline hosts, and
`Coupling.Product` gives access to the determinantal product states.

A `States` object acts as a list of states and `states[i]` returns an
individual`State` object.
The total number of states is obtained from `len(states)` and the following
code prints the full string representation of every state including the 
irreducible representations of all quantum numbers:

```
for state in states:
    print(state)
```

The default full string representation would explicitly be addressed by
`state.long()`, while `state.short()` returns a more compact string using
an integer term index to distinguish different states which share the same
set of quantum numbers $S$, $L$ and $J$.
This index aligns with the choice of Nielson and Koster.

Direct access to the individual quantum numbers of a state is provided by
the dictionary `state.quantum`.
The ìnteger attribute `state.mult` holds the multiplicity $2J+1$
of a state in $SLJ$ coupling and the value 1 for the other coupling schemes. 

The matrix of each tensor operator can be obtained from the `States` object as
`numpy.ndarray` object using its method `matrix`:

```
matrix = states.matrix("H1/2")
```

The rows and columns of this matrix correspond to the respective state indices,
which allows to identify each matrix element clearly.

In $SLJ$ coupling, the `States` object also provides access to reduced matrix
elements:

```
u2_sq = states.reduced("U/2") ** 2
```

Linear combinations of matrices are directly available from both method `matrix`
and `reduced` when the operator name is replaced by a dictionary of weight
factors:

```
from yall import CONST_gs
ls_sq = states.reduced({"L": 1.0, "S": CONST_gs}) ** 2
```

All matrix elements are returned as double-precision floating-point numbers.
The YALL package takes particular care to convert the exact integers from the
AMELI repository, which represent square roots of rationals, into their
floating-point representation, utilising the full resolution of the floating-point
format.

## Class `Levels`

This class is the interface to lanthanide states in intermediate coupling, their
energy levels and radiative transitions.

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

The method `short()` of an `IntermediateState` object returns the short string
representation of its main basis state and `long(min_weight=0.0)` returns a string
containing weight factors together with the respective short string representations
of all basis states with weight factors greater or equal `min_weight`.

The methods `matrix` and `reduced` of an `Levels` object work exactly as for a 
`States` object, but they return matrices in intermediate coupling instead.

If the optional arguments `jo` and `material` are used when a `Levels` object is 
initialised, it can be used to calculate the strength or radiative transitions.
The Judd-Ofelt parameters in the dictionary `jo` are expected in pm<sup>2</sup>
and the object `material` must provide a method `refractive_index(k)` which returns
the refractive index of the material for the given wavenumber in cm<sup>-1</sup>.
The YALL package provides the classes `Cauchy` and `Sellmeier` for this purpose.

There is no universally accepted definition of the line strength of a transition.
For electric dipole transitions, the YALL package uses the definition

$$
S_{ed} = \frac{1}{3(2J_i+1)} \frac{e^2}{4\pi\varepsilon_0} \sum\limits_{\lambda=2,4,6}
\Omega_\lambda |\langle J_j\parallel \mathbf{U}^{(\lambda)}\parallel J_i\rangle|^2
$$

with the Judd-Ofelt parameters $\Omega_2$, $\Omega_4$, and $\Omega_6$.
For magnetic dipole transitions it uses

$$
S_{md} = \frac{1}{3(2J_i+1)} \frac{1}{4\pi\varepsilon_0} \frac{\beta_m^2}{c^2 \hbar^2}
|\langle J_j\parallel \mathbf{L}+g_s\mathbf{S}\parallel J_i\rangle|^2
$$

with the Bohr magneton $\beta_m = \frac{e\hbar}{2m_e}$.
The physical unit of these line strength values is Jm<sup>3</sup>.
The matrices containing all electric and magnetic dipole transition line strengths 
are returned by the method `Levels.line_strengths()` as data class with two
attributes `ed` and `md`.
The following code prints the line strengths of all ground state absorptions:  

```
strength = ion.line_strengths()
print(strength.ed[1:, 0])
print(strength.md[1:, 0])
```

The line strengths are typically used to calculate the dimensionless oscillator strength of a transition $i\to j$
in absorption or emission:

$$
f_{ij} = \frac{4\pi\varepsilon_0}{e^2} \frac{4 \pi m_e c k_{ij}}{\hbar}
[\chi^\prime_{ed}S_{ed} + \chi^\prime_{md}S_{md}]
$$

with the energy level difference $k_{ij}=|k_i-k_j|$ (in wavenumbers) and the local field correction factors

$$
\chi^\prime_{ed} = \frac{(n^2+2)^2}{9n} \qquad \chi^\prime_{md} = n
$$

Another useful quantity is the spontaneous radiative emission rate (in s<sup>-1</sup>) of a transition:

$$
A_{ij} = \frac{32 \pi^3 k_{ij}^3}{\hbar}
\[\chi_{ed}S_{ed} + \chi_{md}S_{md}\]
$$

with the local field correction factors

$$
\chi_{ed} = \frac{n(n^2+2)^2}{9} \qquad \chi_{md} = n^3
$$

Please be aware that due to dispersion the refractive index $n(k_{ij})$ in general is a function of the
wavenumber $k_{ij}$ of the transition.


## Class Fit

## Parameter Sets

## Caching
...

where `"base"` fixes the energy of the ground state and the special parameters `"H5fix"` and `"H6fix"` are
abbreviations for the common choices $M^2 = 0.56 M^0$, $M^4 = 0.38 M^0$ and $P^4 = 0.75 P^2$, $P^6 = 0.50 P^2$.

You may choose operation in either the SLJ or the SLJM space. The latter is used to calculate the Stark splitting
and slows down the calculations significantly, see the section on crystal fields below. In SLJ coupling the quantum
number J of the total angular momentum is a good quantum number and the magnetic quantum number M of the total
angular momentum does not matter. Select SLJ coupling and initialize the Lanthanide ion by giving the number of
4f electrons in the range from 1 to 13:

```
coupling = Coupling.SLJ
ion = Lanthanide(2, coupling, radial)
print(ion)
```

If not specified, the default coupling is SLJ and the default radial parameters are those given by Carnall
for LaF<sub>3</sub>.

The Lanthanide object builds the matrix of the total perturbation Hamiltonian and diagonalises it, which results
in the energy level and the SLJ composition of each state in intermediate coupling.
Each state in intermediate coupling is a mixture of different SLJ states with the same total angular momentum J.
You can access all energies in the list `ion.energies` with the ground state in first position.

Using `Coupling.Intermediate` an object holding all intermediate states is delivered by the method
`ion.states(coupling)`.
The following code prints the energy of the first excited state, and then weight factor and state object of the most
important SLJ component:

```
states = ion.states(Coupling.Intermediate)
state = states[1]
print(state.energy)
print(state.weights[0], state.states[0])
```

For each state object there is a long string representation, which you can access by `str(state)` or
`state.long(min_weight=0.0)` and a short version by `state.short()`.
The parameter `min_weight` is useful for ions with a large number of states.
It gives the minimum weight of a SLJ state to appear in the list.
This allows to show the most important components only.
A shortcut to the list of energies of all states is the attribute `ion.energies`.
The method `ion.str_levels(min_weight=0.0)` provides a convenient way to display the energy level spectrum
on a terminal:

```
for line in ion.str_levels(0.05):
    print(line)
```

The calculation of a new energy level spectrum can be triggered at any time by providing a new set of
radial parameters:

```
from lanthanide import RADIAL
ion.set_radial(RADIAL["Pr3+/ZBLAN"])
```

## Crystal field splitting

The calculation of a full lanthanide spectrum including the Stark splitting is supported by the YALL package and
activated by selecting `Coupling.SLJM`:

```
ion = Lanthanide(2, Coupling.SLJM, radial)
```

The `Lanthanide` object supports real or complex crystal field parameters `Hcf/{k},{q}` with the rank $k=2,4,6$
and $q=0\ldots k$ besides the other radial parameters in the `radial` dictionary. The set of none-zero parameters
is determined by the site symmetry of the lanthanide ion.

Note that crystal field interactions release the J degeneracy of the free ion hamiltonian.
Diagonalisation of the full hamiltonian including crystal field matrices takes place in the full state space and is
thus much slower than in SLJ coupling.
States in intermediate SLJM coupling are linear combinations of all SLJM states, although usually a small number of
dominant SLJM states determines the character of an intermediate state.

## Radiative transitions

For radiative transitions inside the 4f configuration of Lanthanides only electric and magnetic
dipole moments are relevant.
The calculation of the respective transition strengths according to the Judd-Ofelt theory is based on
the reduced matrix elements
$\langle J_j\parallel \mathbf{U}^{(2)}\parallel J_i \rangle$,
$\langle J_j\parallel \mathbf{U}^{(4)}\parallel J_i \rangle$, and
$\langle J_j\parallel \mathbf{U}^{(6)}\parallel J_i \rangle$ for electric and
$\frac{1}{\hbar} \langle J_j\parallel \mathbf{L} + g_s \mathbf{S}\parallel J_i \rangle$
for magnetic dipole transitions.
The method `Lanthanide.line_reduced()` delivers a `Reduced` object, which contains all four squared reduced matrices
as attributes `U2`, `U4`, `U6`, and `LS` as required for the calculation of transition strengths.
The matrix element for a transition from an initial state `i` to a final state `j` is addressed by `array[j,i]`.
This code shows the squared elements
$|\langle J_j\parallel \mathbf{U}^{(4)}\parallel J_0 \rangle|^2$ for transitions from the ground state
to all excited states:

```
reduced = ion.line_reduced()
print(reduced.U4[1:, 0])
```

There is no universally accepted definition of the line strength of a transition. For electric dipole transitions,
the YALL package uses the definition

$$
S_{ed} = \frac{1}{3(2J_i+1)} \frac{e^2}{4\pi\varepsilon_0} \sum\limits_{\lambda=2,4,6}
\Omega_\lambda |\langle J_j\parallel \mathbf{U}^{(\lambda)}\parallel J_i\rangle|^2
$$

with the Judd-Ofelt parameters $\Omega_2$, $\Omega_4$, and $\Omega_6$. For magnetic dipole transitions we use

$$
S_{md} = \frac{1}{3(2J_i+1)} \frac{1}{4\pi\varepsilon_0} \frac{\beta_m^2}{c^2 \hbar^2}
|\langle J_j\parallel \mathbf{L}+g_s\mathbf{S}\parallel J_i\rangle|^2
$$

with the Bohr magneton $\beta_m = \frac{e\hbar}{2m_e}$.
The physical unit of these line strength values is Jm<sup>3</sup>.
For a given set of Judd-Ofelt parameters in pm<sup>2</sup> you can also get the radiative line strength of
each transition using the method `line_strengths()`:

```
judd_ofelt = { "JO/2": 1.981, "JO/4": 4.645, "JO/6": 6.972 }
strength = ion.line_strengths(judd_ofelt)
print(strength.Sed[1:, 0])
print(strength.Smd[1:, 0])
```

## Spectroscopic parameters

The line strengths are typically used to calculate the dimensionless oscillator strength of a transition $i\to j$
in absorption or emission:

$$
f_{ij} = \frac{4\pi\varepsilon_0}{e^2} \frac{4 \pi m_e c k_{ij}}{\hbar}
\[\chi^\prime_{ed}S_{ed} + \chi^\prime_{md}S_{md}\]
$$

with the energy level difference $k_{ij}=|k_i-k_j|$ (in wavenumbers) and the local field correction factors

$$
\chi^\prime_{ed} = \frac{(n^2+2)^2}{9n} \qquad \chi^\prime_{md} = n
$$

Another useful quantity is the spontaneous radiative emission rate (in s<sup>-1</sup>) of a transition:

$$
A_{ij} = \frac{32 \pi^3 k_{ij}^3}{\hbar}
\[\chi_{ed}S_{ed} + \chi_{md}S_{md}\]
$$

with the local field correction factors

$$
\chi_{ed} = \frac{n(n^2+2)^2}{9} \qquad \chi_{md} = n^3
$$

Please be aware that due to dispersion the refractive index $n(k_{ij})$ in general is a function of the
wavenumber $k_{ij}$ of the transition.

## Caching

Matrix elements of angular tensors are physical constants for a given electron configuration.
Therefore, the YALL package stores them after their first calculation in the sub-directory

```
.../Lib/site-packages/lathanide/vaults
```

of the package installation folder.
The calculation can be very time consuming especially for configurations with 4-10 electrons.
You should therefore consider running the the function `create_all()`:

```
from lanthanide import create_all
create_all()
```

It will trigger the calculation of all perturbation hamiltonians and reduced transition matrix elements
for all lanthanide ions.
Expect it to run for more than 24 hours on a typical office PC.
It will generate about 9 GB of data in the folder `vaults`.
There is one cache file for each of the 13 lanthanide configurations.

It is recommended to close the cache file gracefully by calling `Lanthanide.close()` when the work on
one configuration is finished:

```
ion = Lanthanide(2)
print(ion)
energies = ion.energies
reduced = ion.reduced()
ion.close()
```

As an alternative, the Lanthanide class also supports the context management protocol.
You can thus use it to start a `with` block.
The cache file will be closed automatically when the program leaves the block and no manual closing is required:

```
with Lanthanide(2) as ion:
    print(ion)
    energies = ion.energies
    reduced = ion.reduced()
```

The YALL package uses an additional cache for perturbation hamiltonians.
To accelerate energy level fits of radial parameters,
these matrices are kept in the memory after their first retrival from the cache file.
There is also a cache for the reduced transition matrix elements.
Consecutive calls of `Lanthanide.line_reduced()` will deliver the same object until a new set of radial parameters
is applied via `Lanthanide.set_radial()`.

## Usage examples

See Python scripts in the directory `test`.

## Internal structures

If you want to dig deeper into the package and use internal classes and functions, you might find some additional
[documentation](docs/internals.md) in the folder `docs` useful.

## License

This is free software under the MIT License.

## Reference

* <span id="ref1">[1]</span> Reinhard Caspary: "Applied Rare-Earth Spectroscopy for Fiber Laser Optimization", doctoral dissertation at
Technische Universität Braunschweig, published with Shaker, Aachen, 2002
* <span id="ref2">[2]</span> Reinhard Caspary: "AMELI: Angular Matrix Elements of Lanthanide Ions", arXiv, 2026,
https://doi.org/10.48550/arXiv.2603.21947