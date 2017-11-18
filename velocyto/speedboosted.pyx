import cython
import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
from libc.string cimport memset
from libc.math cimport sqrt, fabs, log10
from cython.parallel import threadid, parallel, prange

# Pure C function (written in cython, just for syntactic sugar)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void x_colDeltaCor(double *e,
                       double *d,
                       double *rm,
                       int rows,
                       int cols,
                       int num_threads):
    cdef:
        int i, j, c
    
    with nogil:
        for c in prange(cols, schedule='guided', num_threads=num_threads):
            A = <double *>malloc(rows * cols * sizeof(double))
            tmp = <double *>malloc(1 * sizeof(double)) 
            # subtract the cth column
            for j in range(rows):
                for i in range(cols):
                    A[j*cols + i] = e[j*cols + i] - e[j*cols + c]
            
            #muA = A.mean(0)
            muA = <double *>malloc(cols * sizeof(double))
            memset(muA, 0, cols * sizeof(double))
            for j in range(rows):
                for i in range(cols):
                    muA[i] += A[j*cols + i]
            for i in range(cols):
                muA[i] = muA[i] / rows

            # A_mA = A - muA
            A_mA = <double *>malloc(rows * cols * sizeof(double))
            for j in range(rows):
                for i in range(cols):
                    A_mA[j*cols + i] = A[j*cols + i] - muA[i]

            # mub = b.mean()
            mub = <double *>malloc(1 * sizeof(double))
            mub[0] = 0
            for j in range(rows):
                mub[0] += d[j*cols + c]
            mub[0] = mub[0] / rows

            # b_mb = b - mub
            b_mb = <double *>malloc(rows * sizeof(double))
            for j in range(rows):
                b_mb[j] = d[j*cols + c] - mub[0]

            # ssA = (A_mA**2).sum(0)
            ssA = <double *>malloc(cols * sizeof(double))
            memset(ssA, 0, cols * sizeof(double))
            for j in range(rows):
                for i in range(cols):
                    ssA[i] += A_mA[j*cols + i] * A_mA[j*cols + i]
            for i in range(cols):
                ssA[i] = 1. / sqrt(ssA[i])

            # ssb = (b_mb**2).sum()
            ssb = <double *>malloc(1 * sizeof(double))
            ssb[0] = 0
            for j in range(rows):
                ssb[0] += b_mb[j] * b_mb[j] # **2
            ssb[0] = 1. / sqrt(ssb[0])

            # np.dot(b_mb, A_mA)/(np.sqrt(ssA) * np.sqrt(ssb))
            for j in range(rows):
                tmp[0] = b_mb[j] * ssb[0]
                for i in range(cols):
                    rm[c*cols + i] += (A_mA[j*cols + i] * ssA[i]) * tmp[0]
            
            free(A)
            free(tmp)
            free(muA)
            free(A_mA)
            free(mub)
            free(b_mb)
            free(ssA)
            free(ssb)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void x_colDeltaCorSqrt(double *e,
                            double *d,
                            double *rm,
                            int rows,
                            int cols,
                            int num_threads,
                            double psc):
    cdef:
        int i, j, c
    
    with nogil:
        for c in prange(cols, schedule='guided', num_threads=num_threads):
            A = <double *>malloc(rows * cols * sizeof(double))
            tmp = <double *>malloc(1 * sizeof(double)) 
            # subtract the cth column
            for j in range(rows):
                for i in range(cols):
                    tmp[0] = e[j*cols + i] - e[j*cols + c]
                    if tmp[0] > 0:
                        A[j*cols + i] = sqrt(tmp[0] + psc)
                    else:
                        A[j*cols + i] = -sqrt(-tmp[0] + psc)
            
            #muA = A.mean(0)
            muA = <double *>malloc(cols * sizeof(double))
            memset(muA, 0, cols * sizeof(double))
            for j in range(rows):
                for i in range(cols):
                    muA[i] += A[j*cols + i]
            for i in range(cols):
                muA[i] = muA[i] / rows

            # A_mA = A - muA
            A_mA = <double *>malloc(rows * cols * sizeof(double))
            for j in range(rows):
                for i in range(cols):
                    A_mA[j*cols + i] = A[j*cols + i] - muA[i]

            # mub = b.mean()
            mub = <double *>malloc(1 * sizeof(double))
            mub[0] = 0
            for j in range(rows):
                mub[0] += d[j*cols + c]
            mub[0] = mub[0] / rows

            # b_mb = b - mub
            b_mb = <double *>malloc(rows * sizeof(double))
            for j in range(rows):
                b_mb[j] = d[j*cols + c] - mub[0]

            # ssA = (A_mA**2).sum(0)
            ssA = <double *>malloc(cols * sizeof(double))
            memset(ssA, 0, cols * sizeof(double))
            for j in range(rows):
                for i in range(cols):
                    ssA[i] += A_mA[j*cols + i] * A_mA[j*cols + i]
            for i in range(cols):
                ssA[i] = 1. / sqrt(ssA[i])

            # ssb = (b_mb**2).sum()
            ssb = <double *>malloc(1 * sizeof(double))
            ssb[0] = 0
            for j in range(rows):
                ssb[0] += b_mb[j] * b_mb[j] # **2
            ssb[0] = 1. / sqrt(ssb[0])

            # np.dot(b_mb, A_mA)/(np.sqrt(ssA) * np.sqrt(ssb))
            for j in range(rows):
                tmp[0] = b_mb[j] * ssb[0]
                for i in range(cols):
                    rm[c*cols + i] += (A_mA[j*cols + i] * ssA[i]) * tmp[0]
            
            free(A)
            free(tmp)
            free(muA)
            free(A_mA)
            free(mub)
            free(b_mb)
            free(ssA)
            free(ssb)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void x_colDeltaCorLog10(double *e,
                      double *d,
                      double *rm,
                      int rows,
                      int cols,
                      int num_threads,
                      double psc):
    cdef:
        int i, j, c
    
    with nogil:
        for c in prange(cols, schedule='guided', num_threads=num_threads):
            A = <double *>malloc(rows * cols * sizeof(double))
            tmp = <double *>malloc(1 * sizeof(double)) 
            # subtract the cth column
            for j in range(rows):
                for i in range(cols):
                    tmp[0] = e[j*cols + i] - e[j*cols + c]
                    if tmp[0] > 0:
                        A[j*cols + i] = log10(tmp[0] + psc)
                    else:
                        A[j*cols + i] = -log10(-tmp[0] + psc)
            
            #muA = A.mean(0)
            muA = <double *>malloc(cols * sizeof(double))
            memset(muA, 0, cols * sizeof(double))
            for j in range(rows):
                for i in range(cols):
                    muA[i] += A[j*cols + i]
            for i in range(cols):
                muA[i] = muA[i] / rows

            # A_mA = A - muA
            A_mA = <double *>malloc(rows * cols * sizeof(double))
            for j in range(rows):
                for i in range(cols):
                    A_mA[j*cols + i] = A[j*cols + i] - muA[i]

            # mub = b.mean()
            mub = <double *>malloc(1 * sizeof(double))
            mub[0] = 0
            for j in range(rows):
                mub[0] += d[j*cols + c]
            mub[0] = mub[0] / rows

            # b_mb = b - mub
            b_mb = <double *>malloc(rows * sizeof(double))
            for j in range(rows):
                b_mb[j] = d[j*cols + c] - mub[0]

            # ssA = (A_mA**2).sum(0)
            ssA = <double *>malloc(cols * sizeof(double))
            memset(ssA, 0, cols * sizeof(double))
            for j in range(rows):
                for i in range(cols):
                    ssA[i] += A_mA[j*cols + i] * A_mA[j*cols + i]
            for i in range(cols):
                ssA[i] = 1. / sqrt(ssA[i])

            # ssb = (b_mb**2).sum()
            ssb = <double *>malloc(1 * sizeof(double))
            ssb[0] = 0
            for j in range(rows):
                ssb[0] += b_mb[j] * b_mb[j] # **2
            ssb[0] = 1. / sqrt(ssb[0])

            # np.dot(b_mb, A_mA)/(np.sqrt(ssA) * np.sqrt(ssb))
            for j in range(rows):
                tmp[0] = b_mb[j] * ssb[0]
                for i in range(cols):
                    rm[c*cols + i] += (A_mA[j*cols + i] * ssA[i]) * tmp[0]
            
            free(A)
            free(tmp)
            free(muA)
            free(A_mA)
            free(mub)
            free(b_mb)
            free(ssA)
            free(ssb)

# "partial" functions act on a random sample of the neighbours to reduce computation time
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void x_colDeltaCorpartial(double *e,
                                   double *d,
                                   double *rm,
                                   Py_ssize_t *ixs,
                                   int rows,
                                   int cols,
                                   int nrndm,
                                   int num_threads):
    cdef:
        int i, j, c, n
    
    with nogil:
        for c in prange(cols, schedule='guided', num_threads=num_threads):
            A = <double*> malloc(rows * nrndm * sizeof(double))
            tmp = <double*> malloc(sizeof(double)) 
            # subtract the cth column
            for j in range(rows):
                for n in range(nrndm):
                    i = ixs[c*nrndm + n]
                    A[j*nrndm + n] = e[j*cols + i] - e[j*cols + c]
            
            #muA = A.mean(0)
            muA = <double*> malloc(nrndm * sizeof(double))
            memset(muA, 0, nrndm * sizeof(double))
            for j in range(rows):
                for n in range(nrndm):
                    muA[n] += A[j*nrndm + n]
            for n in range(nrndm):
                muA[n] = muA[n] / rows

            # A_mA = A - muA
            A_mA = <double*> malloc(rows * nrndm * sizeof(double))
            for j in range(rows):
                for n in range(nrndm):
                    A_mA[j*nrndm + n] = A[j*nrndm + n] - muA[n]

            # mub = b.mean()
            mub = <double*> malloc(sizeof(double))
            mub[0] = 0
            for j in range(rows):
                mub[0] += d[j*cols + c]
            mub[0] = mub[0] / rows

            # b_mb = b - mub
            b_mb = <double*> malloc(rows * sizeof(double))
            for j in range(rows):
                b_mb[j] = d[j*cols + c] - mub[0]

            # ssA = (A_mA**2).sum(0)
            ssA = <double*> malloc(nrndm * sizeof(double))
            memset(ssA, 0, nrndm * sizeof(double))
            for j in range(rows):
                for n in range(nrndm):
                    ssA[n] += A_mA[j*nrndm + n] * A_mA[j*nrndm + n]
                    
            # one division and many multiplication (below) faster than many divisions
            for n in range(nrndm):
                ssA[n] = 1. / sqrt(ssA[n])

            # ssb = (b_mb**2).sum()
            ssb = <double*> malloc(sizeof(double))
            ssb[0] = 0
            for j in range(rows):
                ssb[0] += b_mb[j] * b_mb[j] # **2
            
            # one division and many multiplication (below) faster than many divisions
            ssb[0] = 1. / sqrt(ssb[0]) 

            # np.dot(b_mb, A_mA)/(np.sqrt(ssA) * np.sqrt(ssb))
            for j in range(rows):
                tmp[0] = b_mb[j] * ssb[0]
                for n in range(nrndm):
                    i = ixs[c*nrndm + n]
                    rm[c*cols + i] += (A_mA[j*nrndm + n] * ssA[n]) * tmp[0]
                    
            # Cleanup
            free(A)
            free(tmp)
            free(muA)
            free(A_mA)
            free(mub)
            free(b_mb)
            free(ssA)
            free(ssb)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void x_colDeltaCorSqrtpartial(double *e,
                                   double *d,
                                   double *rm,
                                   Py_ssize_t *ixs,
                                   int rows,
                                   int cols,
                                   int nrndm,
                                   int num_threads,
                                   double psc):
    cdef:
        int i, j, c, n
    
    with nogil:
        for c in prange(cols, schedule='guided', num_threads=num_threads):
            A = <double*> malloc(rows * nrndm * sizeof(double))
            tmp = <double*> malloc(sizeof(double)) 
            # subtract the cth column
            for j in range(rows):
                for n in range(nrndm):
                    i = ixs[c*nrndm + n]
                    tmp[0] = e[j*cols + i] - e[j*cols + c]
                    if tmp[0] >= 0:
                        A[j*nrndm + n] = sqrt(tmp[0] + psc)
                    else:
                        A[j*nrndm + n] = -sqrt(-tmp[0] + psc)
                    
            
            #muA = A.mean(0)
            muA = <double*> malloc(nrndm * sizeof(double))
            memset(muA, 0, nrndm * sizeof(double))
            for j in range(rows):
                for n in range(nrndm):
                    muA[n] += A[j*nrndm + n]
            for n in range(nrndm):
                muA[n] = muA[n] / rows

            # A_mA = A - muA
            A_mA = <double*> malloc(rows * nrndm * sizeof(double))
            for j in range(rows):
                for n in range(nrndm):
                    A_mA[j*nrndm + n] = A[j*nrndm + n] - muA[n]

            # mub = b.mean()
            mub = <double*> malloc(sizeof(double))
            mub[0] = 0
            for j in range(rows):
                mub[0] += d[j*cols + c]
            mub[0] = mub[0] / rows

            # b_mb = b - mub
            b_mb = <double*> malloc(rows * sizeof(double))
            for j in range(rows):
                b_mb[j] = d[j*cols + c] - mub[0]

            # ssA = (A_mA**2).sum(0)
            ssA = <double*> malloc(nrndm * sizeof(double))
            memset(ssA, 0, nrndm * sizeof(double))
            for j in range(rows):
                for n in range(nrndm):
                    ssA[n] += A_mA[j*nrndm + n] * A_mA[j*nrndm + n]
                    
            # one division and many multiplication (below) faster than many divisions
            for n in range(nrndm):
                ssA[n] = 1. / sqrt(ssA[n])

            # ssb = (b_mb**2).sum()
            ssb = <double*> malloc(sizeof(double))
            ssb[0] = 0
            for j in range(rows):
                ssb[0] += b_mb[j] * b_mb[j] # **2
            
            # one division and many multiplication (below) faster than many divisions
            ssb[0] = 1. / sqrt(ssb[0]) 

            # np.dot(b_mb, A_mA)/(np.sqrt(ssA) * np.sqrt(ssb))
            for j in range(rows):
                tmp[0] = b_mb[j] * ssb[0]
                for n in range(nrndm):
                    i = ixs[c*nrndm + n]
                    rm[c*cols + i] += (A_mA[j*nrndm + n] * ssA[n]) * tmp[0]
                    
            # Cleanup
            free(A)
            free(tmp)
            free(muA)
            free(A_mA)
            free(mub)
            free(b_mb)
            free(ssA)
            free(ssb)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void x_colDeltaCorLog10partial(double *e,
                                   double *d,
                                   double *rm,
                                   Py_ssize_t *ixs,
                                   int rows,
                                   int cols,
                                   int nrndm,
                                   int num_threads,
                                   double psc):
    cdef:
        int i, j, c, n
    
    with nogil:
        for c in prange(cols, schedule='guided', num_threads=num_threads):
            A = <double*> malloc(rows * nrndm * sizeof(double))
            tmp = <double*> malloc(sizeof(double)) 
            # subtract the cth column
            for j in range(rows):
                for n in range(nrndm):
                    i = ixs[c*nrndm + n]
                    tmp[0] = e[j*cols + i] - e[j*cols + c]
                    if tmp[0] >= 0:
                        A[j*nrndm + n] = log10(tmp[0] + psc)
                    else:
                        A[j*nrndm + n] = -log10(-tmp[0] + psc)
                    
            
            #muA = A.mean(0)
            muA = <double*> malloc(nrndm * sizeof(double))
            memset(muA, 0, nrndm * sizeof(double))
            for j in range(rows):
                for n in range(nrndm):
                    muA[n] += A[j*nrndm + n]
            for n in range(nrndm):
                muA[n] = muA[n] / rows

            # A_mA = A - muA
            A_mA = <double*> malloc(rows * nrndm * sizeof(double))
            for j in range(rows):
                for n in range(nrndm):
                    A_mA[j*nrndm + n] = A[j*nrndm + n] - muA[n]

            # mub = b.mean()
            mub = <double*> malloc(sizeof(double))
            mub[0] = 0
            for j in range(rows):
                mub[0] += d[j*cols + c]
            mub[0] = mub[0] / rows

            # b_mb = b - mub
            b_mb = <double*> malloc(rows * sizeof(double))
            for j in range(rows):
                b_mb[j] = d[j*cols + c] - mub[0]

            # ssA = (A_mA**2).sum(0)
            ssA = <double*> malloc(nrndm * sizeof(double))
            memset(ssA, 0, nrndm * sizeof(double))
            for j in range(rows):
                for n in range(nrndm):
                    ssA[n] += A_mA[j*nrndm + n] * A_mA[j*nrndm + n]
                    
            # one division and many multiplication (below) faster than many divisions
            for n in range(nrndm):
                ssA[n] = 1. / sqrt(ssA[n])

            # ssb = (b_mb**2).sum()
            ssb = <double*> malloc(sizeof(double))
            ssb[0] = 0
            for j in range(rows):
                ssb[0] += b_mb[j] * b_mb[j] # **2
            
            # one division and many multiplication (below) faster than many divisions
            ssb[0] = 1. / sqrt(ssb[0]) 

            # np.dot(b_mb, A_mA)/(np.sqrt(ssA) * np.sqrt(ssb))
            for j in range(rows):
                tmp[0] = b_mb[j] * ssb[0]
                for n in range(nrndm):
                    i = ixs[c*nrndm + n]
                    rm[c*cols + i] += (A_mA[j*nrndm + n] * ssA[n]) * tmp[0]
                    
            # Cleanup
            free(A)
            free(tmp)
            free(muA)
            free(A_mA)
            free(mub)
            free(b_mb)
            free(ssA)
            free(ssb)
    

# Functions accessible from python, adapt the C function to python thorugh cPython API
def _colDeltaCor(double[:, ::1] e,
                 double[:, ::1] d,
                 double[:, ::1] rm,
                 int num_threads):
    cdef:
        int rows = e.shape[0]
        int cols = e.shape[1]
        
    x_colDeltaCor(&e[0,0], &d[0,0], &rm[0,0], rows, cols, num_threads)

def _colDeltaCorSqrt(double[:, ::1] e,
                     double[:, ::1] d,
                     double[:, ::1] rm,
                     int num_threads,
                     double psc):
    cdef:
        int rows = e.shape[0]
        int cols = e.shape[1]
        
    x_colDeltaCorSqrt(&e[0,0], &d[0,0], &rm[0,0], rows, cols, num_threads, psc)

def _colDeltaCorLog10(double[:, ::1] e,
                     double[:, ::1] d,
                     double[:, ::1] rm,
                     int num_threads,
                     double psc):
    cdef:
        int rows = e.shape[0]
        int cols = e.shape[1]
        
    x_colDeltaCorLog10(&e[0,0], &d[0,0], &rm[0,0], rows, cols, num_threads, psc)

def _colDeltaCorpartial(double[:, ::1] e,
                        double[:, ::1] d,
                        double[:, ::1] rm,
                        Py_ssize_t[:, ::1] ixs,
                        int num_threads,
                        double psc):
    cdef:
        int rows = e.shape[0]
        int cols = e.shape[1]
        int nrndm = ixs.shape[1]
        
    x_colDeltaCorpartial(&e[0,0], &d[0,0], &rm[0,0], &ixs[0,0], rows, cols, nrndm, num_threads)

def _colDeltaCorSqrtpartial(double[:, ::1] e,
                            double[:, ::1] d,
                            double[:, ::1] rm,
                            Py_ssize_t[:, ::1] ixs,
                            int num_threads,
                            double psc):
    cdef:
        int rows = e.shape[0]
        int cols = e.shape[1]
        int nrndm = ixs.shape[1]
        
    x_colDeltaCorSqrtpartial(&e[0,0], &d[0,0], &rm[0,0], &ixs[0,0], rows, cols, nrndm, num_threads, psc)

def _colDeltaCorLog10partial(double[:, ::1] e,
                            double[:, ::1] d,
                            double[:, ::1] rm,
                            Py_ssize_t[:, ::1] ixs,
                            int num_threads,
                            double psc):
    cdef:
        int rows = e.shape[0]
        int cols = e.shape[1]
        int nrndm = ixs.shape[1]
        
    x_colDeltaCorLog10partial(&e[0,0], &d[0,0], &rm[0,0], &ixs[0,0], rows, cols, nrndm, num_threads, psc)