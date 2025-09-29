# file "pygoertzel.py"

# BSD 2-Clause License
# Copyright (c) 2012, Pavel Rajmic
# Copyright (c) 2019, David Crist
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.

# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# convenient numpy imports
im = 1j # FIXME: where should numpy imports be for other code?
from numpy import array, ndarray, zeros, zeros_like, ones, arange
from numpy import pi, sin, cos, exp, log10, clip, var, mean, argmax
from numpy import real, imag

def goertzel_general_shortened(x, indvec): #; timevec=[len(x)], flipx=false)
    # GOERTZEL_GENERAL_SHORTENED(X,INDVEC) computes DTFT of one-dimensional
    # signal X at 'indices' contained in INDVEC, using the generalized (and shortened)
    # second-order Goertzel algorithm.
    # Thanks to the generalization, the 'indices' can be non-integer valued
    # in the range 0 to N-1, where N is the len of X.
    # (Index 0 corresponds to the DC component.)
    # Integers in INDVEC result in the classical DFT coefficients.
    #
    # The output is a column complex vector of len len(INDVEC) containing
    # the desired DTFT values.
    #
    # See also: goertzel_classic.

    # (c) 2009-2012, Pavel Rajmic, Brno University of Technology, Czech Rep.

    # @show size(x)
    assert type(x) in (ndarray, array)
    assert x.ndim == 1
    x = x.astype('float32')
    lx = len(x)

    ## Initialization
    no_freq = len(indvec); #number of frequencies to compute
    y = zeros(no_freq, dtype='complex64'); #memory allocation for the output coefficients

    ## Computation via second-order system
    # loop over the particular frequencies
    for cnt_freq in range(no_freq):
        #for a single frequency
        #a/ precompute the constants
        pik_term = 2*pi*(indvec[cnt_freq])/(lx)
        cos_pik_term2 = cos(pik_term) * 2
        cc = exp(-im*pik_term) # complex constant
        #b/ state variables
        s0 = 0
        s1 = 0
        s2 = 0
        #c/ 'main' loop
        for ind in range(lx-1): #number of iterations is (by one) less than the len of signal
            #new state
            s0 = x[ind] + cos_pik_term2 * s1 - s2  # (*)
            #shifting the state variables
            s2 = s1
            s1 = s0
        # end for ind in range(lx-1)

        #d/ final computations
        s0 = x[int(lx-1)] + cos_pik_term2 * s1 - s2 #correspond to one extra performing of (*)

        # @show map(typeof, (s0, s1, cc))
        y[cnt_freq] = s0 - s1*cc #resultant complex coefficient
        # print("PYTHON, TC1, {}, {}, {}".format(cnt_freq, real(y[cnt_freq]), imag(y[cnt_freq])))

        #complex multiplication substituting the last iteration
        #and correcting the phase for (potentially) non-integer valued
        #frequencies at the same time
        y[cnt_freq] = y[cnt_freq] * exp(-im*pik_term*(lx-1))
    # end for cnt_freq in range(no_freq)

    return y # y[:]
# end goertzel_general_shortened


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
# goertzel1D = goertzel_general_shortened # alias
def goertzel1D(x, fs, hzvec):
    x = x.astype('float32')
    hz_per_bin = fs/len(x)
    bin_per_hz = 1/hz_per_bin
    indvec = bin_per_hz * hzvec
    return goertzel_general_shortened(x, indvec)



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
    x = x.astype('float32')
    y = zeros( (len(hzvec), x.shape[1]), dtype='complex64' )
    for cix in range(x.shape[1]):
        y[:, cix] = goertzel1D(x[:, cix], fs, hzvec)

    return y
# end goertzel2D

# END GOERTZEL PORT