// file cgoertzel.c

// MIT License
// Copyright (c) [2019] [David Crist]

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <math.h> // needed for sin/cos and M_PI for goertzel
#ifndef M_PI
#define M_PI (3.14159265359)
#endif

typedef struct _CPX {
    double r;
    double i;
} CPX;

typedef struct goertzel1D_t {
    int NBINS;
    double *BIN;
    double *A;
    double *B;
    CPX *C;
    CPX *D;
} goertzel1D_t;

int init_goertzel1D(goertzel1D_t *g,  double fs, int NSAMP, double *hzvec, int NBINS, *snapvec, int NSNAP);
int free_goertzel1D(goertzel1D_t *g);
int run_goertzel1D(goertzel1D_t *g, double *in, int NSAMP, CPX *out);

// below code written strictly according to figure 4 of 
// https://asp-eurasipjournals.springeropen.com/track/pdf/10.1186/1687-6180-2012-56/

int init_goertzel1D(goertzel1D_t *g,  double fs, int NSAMP, double *hzvec, int NBINS, *snapvec, int NSNAP) {
    int ix;
    double fNSAMP = (double)(NSAMP);

    // set NBINS
    g->NBINS = NBINS;
    g->NSNAP = NSNAP;
    
    // malloc buffers for pointers: HZ, BIN, A, B, C, D
    g->BIN  = (double *) (malloc( sizeof(double) * g->NBINS ));
    g->SNAP = (int    *) (malloc( sizeof(int   ) * g->NSNAP ));
    g->A    = (double *) (malloc( sizeof(double) * g->NBINS ));
    g->B    = (double *) (malloc( sizeof(double) * g->NBINS ));
    g->C    = (CPX    *) (malloc( sizeof(CPX   ) * g->NBINS ));
    g->D    = (CPX    *) (malloc( sizeof(CPX   ) * g->NBINS ));

    // initialize goertzel1D_t vectors
    for(ix=0; ix<g->NBINS; ix++) {
        g->BIN[ix] = hzvec[ix] / (fs/fNSAMP); // (Hz) / (Hz/Bin) => Bin
        g->A[ix] = 2 * M_PI * g->BIN[ix] / fNSAMP;
        g->B[ix] = 2 * cos(g->A[ix]);
        g->C[ix].r = cos(-(g->A[ix])); // real component of exp(-jA)
        g->C[ix].i = sin(-(g->A[ix])); // imag component of exp(-jA)
        g->D[ix].r = cos( -2 * M_PI * g->BIN[ix] * (fNSAMP-1.0) / fNSAMP); // real component of exp(-j [2*M_PI*k*(N-1)/N)] )
        g->D[ix].i = sin( -2 * M_PI * g->BIN[ix] * (fNSAMP-1.0) / fNSAMP); // imag component of exp(-j [2*M_PI*k*(N-1)/N)] )
        // printf("ix=%02d, hz=%2.3f, bin=%2.3f, A=%2.6f, B=%2.6f, C=[%2.6f, %2.6f], D=[%2.6f, %2.6f]\n", 
        //         ix, hzvec[ix], g->BIN[ix], g->A[ix], g->B[ix], g->C[ix].r, g->C[ix].i, g->D[ix].r, g->D[ix].i);
    }
    for(ix=0; ix<g->NSNAP; ix++) {
        g->SNAP[ix] = snapvec[ix];
        // check that g->SNAP is monotonically increasing
        if( (ix > 0) && (g->SNAP[ix-1] >= g->SNAP[ix]) ) {
            return 1; // g->SNAP violates monotonically increasing requirement
        }

        // check that g->SNAP is less than num. of samples, g->NSAMP
        if( g->SNAP[ix] >= g->NSAMP ) {
            return 2;
        }
    }

    return 0;
}

int free_goertzel1D(goertzel1D_t *g) {
    // free buffers for pointers: HZ, BIN, A, B, C, D
    free(g->BIN );
    free(g->SNAP);
    free(g->A   );
    free(g->B   );
    free(g->C   );
    free(g->D   );
    return 0;
}

int run_goertzel1D(goertzel1D_t *g, double *in, int NSAMP, CPX *out) {
    double s0, s1, s2;
    CPX tc1, tc2; // tcX = "temporary, complex" variable

    int bix, six, snpix=0, outix=0; // loop variables
    for(bix=0; bix<(g->NBINS); bix++) { // bix = "bin-index", added loop for multiple frequencies
        s0 = 0.0;
        s1 = 0.0;
        s2 = 0.0;

        // 'Main loop' from paper
        for(six=0; six<(NSAMP-1); six++) { // six = "samp-index"
            s0 = in[six] + (g->B[bix]*s1) - s2;
            s2 = s1;
            s1 = s0;

            if(six == g->SNAP[snpix]-1) {
                snpix++;

                // 'Finalizing calculations' from paper
                s0t = in[NSAMP-1] + (g->B[bix]*s1) - s2;   // 1-of-3
                
                tc1.r = s0t - s1*(g->C[bix].r); // 2-of-3, real
                tc1.i = -(    s1*(g->C[bix].i)); // 2-of-3, imag

                // Assign result directly to output        
                tc2.r = (tc1.r*(g->D[bix].r)) - (tc1.i*(g->D[bix].i)); // 3-of-3, real
                tc2.i = (tc1.r*(g->D[bix].i)) + (tc1.i*(g->D[bix].r)); // 3-of-3, imag

                out[outix].r = tc2.r;
                out[outix].i = tc2.i;
                outix++;
            }
        }
    }
    return 0; // FIXME, other error codes?
}


int goertzel1D(double *invec, int NSAMP, double fs, double *hzvec, int NBINS, double *cpx_out) {
    goertzel1D_t _g1d;
    goertzel1D_t *g1d = &(_g1d);

    int ret;

    ret = init_goertzel1D(g1d, fs, NSAMP, hzvec, NBINS);
    if( 0 != ret) {
        printf("ERR: goertzel1D call to init_goertzel1D = %d\n", ret);
        return ret;
    }

    ret = run_goertzel1D(g1d, invec, NSAMP, (CPX*)(cpx_out));
    if( 0 != ret) {
        printf("ERR: goertzel1D call to run_goertzel1D = %d\n", ret);
        return ret;
    }

    ret = free_goertzel1D(g1d);
    if( 0 != ret) {
        printf("ERR: goertzel1D call to free_goertzel1D = %d\n", ret);
        return ret;
    }

    return 0;
}
