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
\to \mathbf{L}^2 \to \mathbf{J}^2 \to \mathbf{J}_z $$

with the operators of total spin $\mathbf{S}$, total orbital angular momentum
$\mathbf{L}$, total angular momentum $\mathbf{J}$ and it's magnetic component
$\mathbf{J}_z$.
Furthermore, there is Casimir's quadratic operator
$\mathbf{C}_2$ of the special orthogonal group $\mathrm{SO}(7)$ in seven
dimensions $(2l+1)$ and the special group $\mathrm{G}_2$.

### Available Operators

The following operator matrices are directly available from the AMELI repository:

| Filename       | Description                                                                                         | Tensor Operator                                   | Radial Integral | Rank   | Parameters                       |
|----------------|-----------------------------------------------------------------------------------------------------|---------------------------------------------------|-----------------|--------|----------------------------------|
| `U/{k},{q}`    | Component $q$ of the total unit tensor operator of rank $k$ in the orbital angular momentum space   | $\mathrm{{U}}^{{(k)}}_q$                          |                 | $k$    | $k=0\ldots 2l+1, q=-k\ldots k$   |
| `T/{k},{q}`    | Component $q$ of the total unit tensor operator of rank $k$ in the spin space                       | $\mathrm{{T}}^{{(k)}}_q$                          |                 | $k$    | $k=0, 1,\ q=-k\ldots k$          |
| `UU/{k}`       | Squared total unit tensor operator of rank $k$ in the orbital angular momentum space                | $(\mathrm{{U}}^{{(k)}}\cdot\mathrm{{U}}^{{(k)}})$ |                 | 0      | $k=0\ldots 2l+1$                 |
| `TT/{k}`       | Squared total unit tensor operator of rank $k$ in the spin space                                    | $(\mathrm{{T}}^{{(k)}}\cdot\mathrm{{T}}^{{(k)}})$ |                 | 0      | $k=0,1$                          |
| `UT/{k}`       | Scalar product of the total unit tensor operators of rank $k$ in the orbital and spin spaces        | $(\mathrm{{U}}^{{(k)}}\cdot\mathrm{{T}}^{{(k)}})$ |                 | 0      | $k=0,1$                          |
| `C/{k},{q}`    | Component $q$ of the Coulomb operator of rank $k$                                                   | $\mathrm{{C}}^{{(k)}}_q$                          |                 | $k$    | $k=0,2,4,6,\ q=-k\ldots k$       |
||
| `L/{q}`        | Component $q$ of the total orbital angular momentum operator                                        | $\mathrm{{L}}_q$                                  |                 | 1      | $q=-1,0,1$                       |
| `S/{q}`        | Component $q$ of the total spin angular momentum operator                                           | $\mathrm{{S}}_q$                                  |                 | 1      | $q=-1,0,1$                       |
| `J/{q}`        | Component $q$ of the total angular momentum operator                                                | $\mathrm{{J}}_q$                                  |                 | 1      | $q=-1,0,1$                       |
| `LL`           | Squared total orbital angular momentum operator                                                     | $(\mathrm{{L}}\cdot\mathrm{{L}})$                 | $\alpha$        | 0      |                                  |
| `SS`           | Squared total spin angular momentum operator                                                        | $(\mathrm{{S}}\cdot\mathrm{{S}})$                 |                 | 0      |                                  |
| `JJ`           | Squared total angular momentum operator                                                             | $(\mathrm{{J}}\cdot\mathrm{{J}})$                 |                 | 0      |                                  |
| `LS`           | Scalar product of the total orbital and spin angular momentum operators                             | $(\mathrm{{L}}\cdot\mathrm{{S}})$                 |                 | 0      |                                  |
||
| `C2`           | Casimir operator of the special group $G_2$                                                         | $\mathrm{{C}}_2(G_2)$                             | $\beta$         | 0      |                                  |
| `CR`           | Casimir operator of the special orthogonal (rotational) group in $2l+1$ dimensions                  | $\mathrm{{C}}_2(SO(2l+1))$                        | $\gamma$        | 0      |                                  |
| `H1/{k}`       | Coulomb first order perturbation Hamiltonian of rank $k$                                            | $\mathrm{{f}}_k$                                  | $F^k$           | 0      | $k=0,2,4,6$                      |
| `H2`           | Spin-orbit first order perturbation Hamiltonian                                                     | $\mathrm{{z}}$                                    | $\zeta$         | 0      |                                  |
| `H4/{c}`       | Effective Coulomb second order perturbation Hamiltonian                                             | $\mathrm{{t}}_c$                                  | $T^c$           | 0      | $c=1\ldots 9$                    |
| `Hss/{k}`      | Spin-spin first order perturbation Hamiltonian of rank $k$                                          | $\mathrm{{m}}_{{k,ss}}$                           |                 | 0      | $k=0,2,4$                        |
| `Hsoo/{k}`     | Spin-other-orbit first order perturbation Hamiltonian of rank $k$                                   | $\mathrm{{m}}_{{k,soo}}$                          |                 | 0      | $k=0,2,4$                        |
| `H5/{k}`       | Spin-spin and spin-other-orbit first order perturbation Hamiltonian of rank $k$                     | $\mathrm{{m}}_k$                                  | $M^k$           | 0      | $k=0,2,4$                        |
| `H6/{k}`       | Effective electrostatic spin-orbit second order perturbation Hamiltonian of rank $k$                | $\mathrm{{p}}_k$                                  | $P^k$           | 0      | $k=2,4,6$                        |
| `Hcf/{k},{q}`  | Component $q$ of the crystal field perturbation Hamiltonian of rank $k$                             | $\mathrm{{H}}_{{\mathrm{{cf}},q}}^{{(k)}}$        | $B^k_q$         | $k$    | $k=0,2,4,6,\ q=0\ldots k$        |

### Alternative Names

A couple of convenient alternative operator names can be used with the YALIP
package.
Their translation into AMELI names is given in the following table:

| YALIP   | AMELI                                       |
|---------|---------------------------------------------|
| `Lz`    | `L/0`                                       |
| `Sz`    | `S/0`                                       |
| `Jz`    | `J/0`                                       |
| `L2`    | `LL`                                        |
| `S2`    | `SS`                                        |
| `J2`    | `JJ`                                        |
| `C7`    | `CR`                                        |
| `C5`    | `CR`                                        |
| `H3/0`  | `LL`                                        |
| `H3/1`  | `C2`                                        |
| `H3/2`  | `CR`                                        |
| `H5fix` | `{"H5/0": 1.0, "H5/2": 0.56, "H5/4": 0.38}` |
| `H6fix` | `{"H6/2": 1.0, "H6/4": 0.75, "H6/6": 0.50}` |

### Perturbation Hamiltonians

The centerpiece of AMELI and YALIP are six free-ion perturbation Hamiltonians
resembling a total of 19 radial integrals.
They are used to calculate the energy levels of a given Lathanide ion.
In the order of their relative magnitude they are:

1. Coulomb interaction between the electrons (1st order):
   $\mathbf{H}_1 = F^2 \mathbf{f}_2 + F^4 \mathbf{f}_4 + F^6 \mathbf{f}_6$,
   with the radial integrals $F^2$, $F^4$, and $F^6$
   and the respective angular two-electron operators $\mathbf{f}_2$,
   $\mathbf{f}_4$, and $\mathbf{f}_6$.
   The YALIP package uses the keys `"H1/2"`, `"H1/4"`, and `"H1/6"` for radial 
   parameters (integrals) and angular matrices (operators).
2. Magnetic spin-orbit interaction of each electron (1st order):
   $\mathbf{H}_2 = \zeta \mathbf{z}$,
   with the radial integral $\zeta$ and the angular one-electron operator
   $\mathbf{z}$.
   The YALIP package uses the key `"H2"` for the radial parameter and the angular
   matrix.
3. Coulomb inter-configuration interactions (2nd order):
   $\mathbf{H}_3 = \alpha \mathbf{L}^2 + \beta \mathbf{C}_2(\mathrm{G}_2) + \gamma \mathbf{C}_2(\mathrm{SO}(7))$,
   with the radial integrals $\alpha$, $\beta$, and $\gamma$.
   The respective effective angular two-electron operators $\mathbf{L}^2$,
   $\mathbf{C}_2(\mathrm{G}_2)$, and $\mathbf{C}_2(\mathrm{SO}(7))$ are the squared operator of the total
   orbital angular momentum and the mentioned quadratic Casimir operators.
   The YALIP package uses the keys `"H3/0"`, `"H3/1"`, and `"H3/2"` for radial
   parameters and angular matrices.
4. More Coulomb inter-configuration interactions (2nd order):
   $\mathbf{H}_4 = \sum_c T^c \mathbf{t}_c$, with $c = 2, 3, 4, 6, 7, 8$.
   The radial integrals are $T^c$ and the respective effective angular three-electron operators
   $\mathbf{t}_c$.
   The YALIP package uses the keys `"H4/2"`, `"H4/3"`, `"H4/4"`, `"H4/6"`, `"H4/7"`, and `"H4/8"` for
   radial parameters and angular matrices.
5. Magnetic interactions of the spin of one electron with the spin (ss) or the orbital angular momentum (soo)
   of another electron (1st order):
   $\mathbf{H}_5 = M^0 \mathbf{m}_0 + M^2 \mathbf{m}_2 + M^4 \mathbf{m}_4$,
   with the Marvin integrals $M^0$, $M^2$, and $M^4$ and the respective angular two-electron operators
   $\mathbf{m}_0$, $\mathbf{m}_2$, and $\mathbf{m}_4$.
   The YALIP package uses the keys `"H5/0"`, `"H5/2"`, and `"H5/4"` for radial parameters and angular matrices.
6. Magnetic inter-configuration spin-orbit interactions (2nd order):
   $\mathbf{H}_6 = P^2 \mathbf{p}_2 + P^4 \mathbf{p}_4 + P^6 \mathbf{p}_6$,
   with the radial integrals $P^2$, $P^4$, and $P^6$ and the respective angular two-electron operators
   $\mathbf{p}_2$, $\mathbf{p}_4$, and $\mathbf{p}_6$.
   The YALIP package uses the keys `"H6/2"`, `"H6/4"`, and `"H6/6"` for radial parameters and angular matrices.

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

### Alternative Radial Integrals

The YALIP package supports some definitions of radial parameters which are often found in older literature.
They are automatically converted into the standard set when you use them in your dictionary of radial integrals.
The attribute `Levels.radial_integrals` always contains the standard set of converted values.

The radial integrals $F_0$, $F_2$, $F_4$, and $F_6$ with sub-script instead of super-script were used in the past to
avoid fractional coefficients in the numeric calculations.
You can specify them using the keys `F_0`, `F_2`, `F_4`, and `F_6`.
The same holds for $P_2$, $P_4$, and $P_6$ which you can use with the keys `P_2`, `P_4`, and `P_6`.

The first order Coulomb interaction parameters $E^0$, $E^1$, $E^2$, and $E^3$ are somewhat different, because each of
them is a linear combination of the standard set $F^0$, $F^2$, $F^4$, and $F^6$.
You can use them with the keys `E^0`, `E^1`, `E^2`, and `E^3`.