# `cgoertzel`
### a Goertzel algorithm written in C with a NumPy wrapper

This contains a simple implementation of a Discrete Time Fourier Transform (DTFT). The C implementation uses what is known as
> 'Generalized Goertzel algorithm with shortened iteration loop'

This term comes from the paper
> [Sysel and Rajmic: Goertzel algorithm generalized to non-integer multiples of fundamental frequency. EURASIP Journal on Advances in Signal Processing 2012 2012:56](https://asp-eurasipjournals.springeropen.com/track/pdf/10.1186/1687-6180-2012-56/)

The C implementation is contained in `cgoertzel.c`. The Python wrapper around that C implementation is contained in `cgoertzel.py`. The script to build the code is `cgoertzel_build.py`.

This is currently not a polished package, but rather a simple utility in the interest of flexibility and speed.

## Getting Started

To build the code, change to the `cgoertzel` directory and run
```
python3 cgoertzel_build.py
```

That should leave one or more binary files named something like `_cgoertzel.cpython-37m-darwin.so` in the same directory named `cgoertzel`.

There's one main function, which is `dtft_bins = goertzel1D(x, fs, hz)`
* `x` is the input vector to the DTFT/Goertzel Operation
* `fs` is the scalar sample rate of `x`
* `hz` is a vector of frequencies at which to calculate the DTFT/Goertzel coefficients

Once the code is built, you should be able to import it and use it. One example is
```
from cgoertzel import goertzel1D
from numpy.random import randn
x = randn(1000)
fs = 100
hz = [1, 2, 3, 4, 5]
g = goertzel1D(x, fs, hz)
```

### Prerequisites

To build and use `cgoertzel` requires the following libraries
* `python3` with
  * `cffi` with a configured compiler
  * `numpy`

## Authors

* **Dave Crist**

## License

Different files have different licenses, but everything is generally permissive; BSD or MIT.

## Acknowledgments

* Thanks to the authors of the Generalized Goertzel paper:
> [Sysel and Rajmic: Goertzel algorithm generalized to non-integer multiples of fundamental frequency. EURASIP Journal on Advances in Signal Processing 2012 2012:56](https://asp-eurasipjournals.springeropen.com/track/pdf/10.1186/1687-6180-2012-56/)

* Thanks to the developers of `cffi`
