# YALL 0.9.0 (Yet Another Lanthanide Library)

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
However, instead of computing its numerical operator matrices, the YALL package
is using the matrix elements from the Zenodo repository
[AMELI](https://zenodo.org/communities/ameli) calculated and stored in exact
arithmetic according to [[2]](#ref2).

## Installation

The package is available on PyPI and the installation therefore is possible using `pip`

```
pip install yall
```

## Usage

The YALL package is designed for three use cases on different abstraction levels.

### 1. Class `States`

Access to raw numerical representations of states and
[operator matrices](docs/operators.md) from the AMELI repository is provided by
the [class `States`](docs/states.md).
The following code initialises all basis states of the Pr<sup>3+</sup> ion in $SLJ$
coupling:

```
from yall import States, Coupling

config = "f2"
coupling = Coupling.SLJ

states = States(config, coupling)
```

### 2. Class `Levels`

Calculation of states in intermediate coupling, their energy levels and radiative
transitions based on given radial integrals and Judd-Ofelt parameters is provided
by the [class `Levels`](docs/levels.md).
Typical initialisation code for a Pr<sup>3+</sup> ion looks like:

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

### 3. Class `Fits`

The [class `Fits`](docs/fits.md) is used to perform energy level and Judd-Ofelt
fits to determine optimised radial integrals and Judd-Ofelt parameters matching
a measured absorption spectrum.
Typical initialisation code for a Pr<sup>3+</sup> ion looks like:

```
from yall import Cauchy, Coupling, Fits

config = "f2"
coupling = Coupling.SLJ
radial = {"base": 327.39, "H1/2": 68576.05, "H1/4": 49972.76, "H1/6": 32415.29, "H2": 728.18,
          "H3/0": 16.99, "H3/1": -417.98, "H3/2": 1371, "H5fix": 0.19, "H6fix": 1.67}
material = Cauchy(1.35123e-5, 2.94780e-3, 1.49985, -1.30933e-3, -3.23335e-6)

opt = Fits(config, coupling, radial, material)
```

### Example Scripts

Some application example scripts are available in the folder [run](run).

## Logging

The YALL package uses the logger package for status messages.
The following code example provides a basic setup to make these messages visible
on the console:

```
import logging
from yall import Coupling, States

logger = logging.getLogger("my_script")

def init_logger(level=logging.INFO):
    root = logging.getLogger()
    root.setLevel(level)
    log_format = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

    console_h = logging.StreamHandler()
    console_h.setFormatter(log_format)
    console_h.setLevel(level)
    root.addHandler(console_h)

if __name__ == "__main__":
    init_logger(level=logging.DEBUG)

    config = "f2"
    coupling = Coupling.SLJ
    
    states = States(config, coupling)
```

## Caching

The YALL package uses `platformdirs.user_cache_dir()` to create a cache folder
for files from the AMELI repository and another one for converted floating point
matrices.
YALL regularly checks the repository for new versions, but at most twice a day.

The decorator `functools.lru_cache` is used to keep numerical matrices in the
memory. 

## License

This is free software under the MIT License.

## References

<span id="ref1">[1]</span> Reinhard Caspary: "Applied Rare-Earth Spectroscopy for Fiber Laser Optimization", doctoral dissertation at
Technische Universität Braunschweig, published with Shaker, Aachen, 2002

<span id="ref2">[2]</span> Reinhard Caspary: "AMELI: Angular Matrix Elements of Lanthanide Ions", arXiv, 2026,
https://doi.org/10.48550/arXiv.2603.21947