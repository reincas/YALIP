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

### Basis States

A `States` object acts as a list of basis states in the given coupling scheme
and `states[i]` returns an individual`State` object.
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

### Matrices

The matrix of each tensor operator can be obtained from the `States` object as
`numpy.ndarray` object by calling its method `matrix` with the respective operator
name:

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

Linear combinations of matrices are directly available from both methods `matrix`
and `reduced` when the operator name is replaced by a dictionary of weight
factors:

```
from yall import CONST_gs
ls_sq = states.reduced({"L": 1.0, "S": CONST_gs}) ** 2
```

All matrix elements are returned as double-precision floating-point numbers.
The YALL package takes particular care to convert the exact integer representation
of signed square roots of rationals from the AMELI repository into their
floating-point representation with maximum precision, utilising the full resolution
of the floating-point format.

