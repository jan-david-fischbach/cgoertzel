# file "cgoertzel.py"

# MIT License
# Copyright (c) [2019] [David Crist]

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# file "cgoertzel.py"
import numpy as np

# if this fails, the import'ing file should handle the exception
from _cgoertzel import ffi, lib

"""
# To test goertzel1D, try the following
from numpy import arange, logspace, sin, pi, outer, abs
fs, f, nsec = 100, 5, 10; tv = arange(int(fs*nsec))*(1/fs)
hz_to_check = logspace(0, 1) # 1 to 10 Hz
x = sin(2*pi*f*tv)
from pygoertzel import goertzel1D
gabs = abs(goertzel1D(x, fs, hz_to_check))
figure(); plot(hz_to_check, gabs, 'o-'); grid()
"""
def goertzel1D(x, fs, hz, snap=None, reverse=False):
    # Make C compatible versions of NumPy Arrays
    cin = np.array(x).flatten().astype("float64")
    if reverse:
        cin = np.flip(cin, 0)
    cnsamp = ffi.cast("int", len(x))
    cfs = ffi.cast("double", fs)
    chz = np.array(hz).flatten().astype("float64")
    cnhz = ffi.cast("int", len(hz))
    if snap is None:
        snap = [len(cin), ]
    csnap = np.array(snap).flatten().astype("int")
    cnsnap = ffi.cast("int", len(snap))

    cout = np.zeros(len(hz)*len(snap)*2, dtype="float64")

    # Use numpy/ctypes trick to make valid C pointer from numpy array
    # https://stackoverflow.com/questions/16276268/how-to-pass-a-numpy-array-into-a-cffi-function-and-how-to-get-one-back-out
    p_cin   = ffi.cast("double *", ffi.from_buffer(cin))
    p_chz   = ffi.cast("double *", ffi.from_buffer(chz))
    p_csnap = ffi.cast("int    *", ffi.from_buffer(csnap))
    p_cout  = ffi.cast("double *", ffi.from_buffer(cout))

    # Call C function, goertzel1D
    ret = lib.goertzel1D(p_cin, cnsamp, cfs, p_chz, cnhz, p_csnap, cnsnap, p_cout)

    if 0 !=  ret:
        print("ERR: non-zero exit received from C extension: goertzel1D = "+str(ret))
        return ret

    # Create NumPy Array from cout
    pyout = np.reshape( cout[0::2] + (1j * cout[1::2]), (cnsnap, cnhz), order='F' )
    return pyout
# end goertzel1D

"""
# To test goertzel2D, try the following
from numpy import arange, sin, pi, outer, logspace, abs
fs, fv, nsec = 100, [1, 2, 3, 4, 5], 10; tv = arange(20*fs)*(1/fs)
hz_to_check = logspace(0, 1) # 1 to 10 Hz
x = sin(2*pi*outer(tv, fv))
from pygoertzel import goertzel2D
gabs = abs(goertzel2D(x, fs, hz_to_check))
figure(); plot(hz_to_check, gabs, 'o-'); grid()
"""
def goertzel2D(x, fs, hzvec):
    x = x.astype('float64')
    y = np.zeros( (len(hzvec), x.shape[1]), dtype='complex128' )
    for cix in range(x.shape[1]):
        y[:, cix] = goertzel1D(x[:, cix], fs, hzvec)

    return y
# end goertzel2D
