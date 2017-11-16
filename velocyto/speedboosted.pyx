import cython
import numpy as np
cimport numpy as np
from libc.math cimport sqrt, fabs, log10
from cython.parallel import threadid, parallel, prange

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _colDeltaCor(double[:, ::1] e,
                 double[:, ::1] d,
                 double[:, ::1] rm,
                 int num_threads):
    cdef:
        Py_ssize_t rows = e.shape[0]
        Py_ssize_t cols = e.shape[1]
        Py_ssize_t i, j, c, t
        double[:, :, ::1] t_A = np.zeros((num_threads, e.shape[0], e.shape[1])) # np.tile(e, (3,1,1))
        double[:, ::1] t_b = np.array(d.T, order="C")
        double[:, ::1] t_out = np.zeros((num_threads, e.shape[1]))
        double[:, ::1] t_muA = np.zeros((num_threads, e.shape[1]))
        double[:, :, ::1] t_A_mA = np.zeros((num_threads, e.shape[0], e.shape[1]))
        double[:, ::1] t_b_mb = np.zeros((num_threads, d.shape[0]))
        double[:, ::1] t_ssA = np.zeros((num_threads, e.shape[1]))
        double[::1] t_mub = np.zeros(num_threads)
        double[::1] t_ssb = np.zeros(num_threads)
        double[::1] t_tmp = np.zeros(num_threads)
        int thread_id
    
    with nogil, cython.boundscheck(False), cython.wraparound(False), cython.cdivision(True):
        for c in prange(cols, schedule='static', num_threads=num_threads):
            t = threadid() # or 
            
            # subtract the cth column
            for j in range(rows):
                for i in range(cols):
                    t_A[t, j, i] = e[j, i] - e[j, c]
            
            #muA = A.mean(0)
            for j in range(rows):
                for i in range(cols):
                    t_muA[t, i] += t_A[t, j, i]
            for i in range(cols):
                t_muA[t, i] = t_muA[t, i] / rows

            # A_mA = A - muA
            for j in range(rows):
                for i in range(cols):
                    t_A_mA[t, j, i] = t_A[t, j, i] - t_muA[t, i]

            # mub = b.mean()
            for j in range(rows):
                t_mub[t] += t_b[c, j]
            t_mub[t] = t_mub[t] / rows

            # b_mb = b - mub
            for j in range(rows):
                t_b_mb[t, j] = t_b[c, j] - t_mub[t]

            # ssA = (A_mA**2).sum(0)
            for j in range(rows):
                for i in range(cols):
                    t_ssA[t, i] += t_A_mA[t, j, i] * t_A_mA[t, j, i]
            for i in range(cols):
                t_ssA[t, i] = 1. / sqrt(t_ssA[t, i])

            # ssb = (b_mb**2).sum()
            for j in range(rows):
                t_ssb[t] += t_b_mb[t, j] * t_b_mb[t, j] # **2
            t_ssb[t] = 1. / sqrt(t_ssb[t])

            # np.dot(b_mb, A_mA)/(np.sqrt(ssA) * np.sqrt(ssb))
            for j in range(rows):
                t_tmp[t] = t_b_mb[t, j] * t_ssb[t]
                for i in range(cols):
                    rm[c, i] += (t_A_mA[t, j, i] * t_ssA[t, i]) * t_tmp[t]
                    
            # Cleanup
            for i in range(cols):
                t_muA[t, i] = 0
            
            t_mub[t] = 0
                
            for i in range(cols):
                t_ssA[t, i] = 0
                
            t_ssb[t] = 0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _colDeltaCorLog10(double[:, ::1] e,
                      double[:, ::1] d,
                      double[:, ::1] rm,
                      int num_threads,
                      double psc):
    cdef:
        Py_ssize_t rows = e.shape[0]
        Py_ssize_t cols = e.shape[1]
        Py_ssize_t i, j, c, t
        double[:, :, ::1] t_A = np.zeros((num_threads, e.shape[0], e.shape[1])) # np.tile(e, (3,1,1))
        double[:, ::1] t_b = np.array(d.T, order="C")
        double[:, ::1] t_out = np.zeros((num_threads, e.shape[1]))
        double[:, ::1] t_muA = np.zeros((num_threads, e.shape[1]))
        double[:, :, ::1] t_A_mA = np.zeros((num_threads, e.shape[0], e.shape[1]))
        double[:, ::1] t_b_mb = np.zeros((num_threads, d.shape[0]))
        double[:, ::1] t_ssA = np.zeros((num_threads, e.shape[1]))
        double[::1] t_mub = np.zeros(num_threads)
        double[::1] t_ssb = np.zeros(num_threads)
        double[::1] t_tmp = np.zeros(num_threads)
        int thread_id
    
    with nogil, cython.boundscheck(False), cython.wraparound(False), cython.cdivision(True):
        for c in prange(cols, schedule='static', num_threads=num_threads):
            t = threadid() # or 
            
            # subtract the cth column
            for j in range(rows):
                for i in range(cols):
                    t_tmp[t] = e[j, i] - e[j, c]
                    if t_tmp[t] > 0:
                        t_A[t, j, i] = log10(fabs(t_tmp[t]) + psc)
                    else:
                        t_A[t, j, i] = -log10(fabs(t_tmp[t]) + psc)
                    
            
            #muA = A.mean(0)
            for j in range(rows):
                for i in range(cols):
                    t_muA[t, i] += t_A[t, j, i]
            for i in range(cols):
                t_muA[t, i] = t_muA[t, i] / rows

            # A_mA = A - muA
            for j in range(rows):
                for i in range(cols):
                    t_A_mA[t, j, i] = t_A[t, j, i] - t_muA[t, i]

            # mub = b.mean()
            for j in range(rows):
                t_mub[t] += t_b[c, j]
            t_mub[t] = t_mub[t] / rows

            # b_mb = b - mub
            for j in range(rows):
                t_b_mb[t, j] = t_b[c, j] - t_mub[t]

            # ssA = (A_mA**2).sum(0)
            for j in range(rows):
                for i in range(cols):
                    t_ssA[t, i] += t_A_mA[t, j, i] * t_A_mA[t, j, i]
            for i in range(cols):
                t_ssA[t, i] = 1. / sqrt(t_ssA[t, i])

            # ssb = (b_mb**2).sum()
            for j in range(rows):
                t_ssb[t] += t_b_mb[t, j] * t_b_mb[t, j] # **2
            t_ssb[t] = 1. / sqrt(t_ssb[t])

            # np.dot(b_mb, A_mA)/(np.sqrt(ssA) * np.sqrt(ssb))
            for j in range(rows):
                t_tmp[t] = t_b_mb[t, j] * t_ssb[t]
                for i in range(cols):
                    rm[c, i] += (t_A_mA[t, j, i] * t_ssA[t, i]) * t_tmp[t]
                    
            # Cleanup
            for i in range(cols):
                t_muA[t, i] = 0
            
            t_mub[t] = 0
                
            for i in range(cols):
                t_ssA[t, i] = 0
                
            t_ssb[t] = 0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _colDeltaCorSqrt(double[:, ::1] e,
                      double[:, ::1] d,
                      double[:, ::1] rm,
                      int num_threads,
                      double psc):
    cdef:
        Py_ssize_t rows = e.shape[0]
        Py_ssize_t cols = e.shape[1]
        Py_ssize_t i, j, c, t
        double[:, :, ::1] t_A = np.zeros((num_threads, e.shape[0], e.shape[1])) # np.tile(e, (3,1,1))
        double[:, ::1] t_b = np.array(d.T, order="C")
        double[:, ::1] t_out = np.zeros((num_threads, e.shape[1]))
        double[:, ::1] t_muA = np.zeros((num_threads, e.shape[1]))
        double[:, :, ::1] t_A_mA = np.zeros((num_threads, e.shape[0], e.shape[1]))
        double[:, ::1] t_b_mb = np.zeros((num_threads, d.shape[0]))
        double[:, ::1] t_ssA = np.zeros((num_threads, e.shape[1]))
        double[::1] t_mub = np.zeros(num_threads)
        double[::1] t_ssb = np.zeros(num_threads)
        double[::1] t_tmp = np.zeros(num_threads)
        int thread_id
    
    with nogil, cython.boundscheck(False), cython.wraparound(False), cython.cdivision(True):
        for c in prange(cols, schedule='static', num_threads=num_threads):
            t = threadid() # or 
            
            # subtract the cth column
            for j in range(rows):
                for i in range(cols):
                    t_tmp[t] = e[j, i] - e[j, c]
                    if t_tmp[t] > 0:
                        t_A[t, j, i] = log10(fabs(t_tmp[t]) + psc)
                    else:
                        t_A[t, j, i] = -log10(fabs(t_tmp[t]) + psc)
                    
            
            #muA = A.mean(0)
            for j in range(rows):
                for i in range(cols):
                    t_muA[t, i] += t_A[t, j, i]
            for i in range(cols):
                t_muA[t, i] = t_muA[t, i] / rows

            # A_mA = A - muA
            for j in range(rows):
                for i in range(cols):
                    t_A_mA[t, j, i] = t_A[t, j, i] - t_muA[t, i]

            # mub = b.mean()
            for j in range(rows):
                t_mub[t] += t_b[c, j]
            t_mub[t] = t_mub[t] / rows

            # b_mb = b - mub
            for j in range(rows):
                t_b_mb[t, j] = t_b[c, j] - t_mub[t]

            # ssA = (A_mA**2).sum(0)
            for j in range(rows):
                for i in range(cols):
                    t_ssA[t, i] += t_A_mA[t, j, i] * t_A_mA[t, j, i]
            for i in range(cols):
                t_ssA[t, i] = 1. / sqrt(t_ssA[t, i])

            # ssb = (b_mb**2).sum()
            for j in range(rows):
                t_ssb[t] += t_b_mb[t, j] * t_b_mb[t, j] # **2
            t_ssb[t] = 1. / sqrt(t_ssb[t])

            # np.dot(b_mb, A_mA)/(np.sqrt(ssA) * np.sqrt(ssb))
            for j in range(rows):
                t_tmp[t] = t_b_mb[t, j] * t_ssb[t]
                for i in range(cols):
                    rm[c, i] += (t_A_mA[t, j, i] * t_ssA[t, i]) * t_tmp[t]
                    
            # Cleanup
            for i in range(cols):
                t_muA[t, i] = 0
            
            t_mub[t] = 0
                
            for i in range(cols):
                t_ssA[t, i] = 0
                
            t_ssb[t] = 0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _colDeltaCorpartial(double[:, ::1] e,
                        double[:, ::1] d,
                        double[:, ::1] rm,
                        Py_ssize_t[:, ::1] ixs,
                        int num_threads):
    cdef:
        Py_ssize_t rows = e.shape[0]
        Py_ssize_t cols = e.shape[1]
        Py_ssize_t nrndm = ixs.shape[1]
        Py_ssize_t i, j, c, t, n
        double[:, :, ::1] t_A = np.zeros((num_threads, e.shape[0], nrndm)) # np.tile(e, (3,1,1))
        double[:, ::1] t_b = np.array(d.T, order="C")
        double[:, ::1] t_out = np.zeros((num_threads, nrndm))
        double[:, ::1] t_muA = np.zeros((num_threads, nrndm))
        double[:, :, ::1] t_A_mA = np.zeros((num_threads, e.shape[0], nrndm))
        double[:, ::1] t_b_mb = np.zeros((num_threads, d.shape[0]))
        double[:, ::1] t_ssA = np.zeros((num_threads, nrndm))
        double[::1] t_mub = np.zeros(num_threads)
        double[::1] t_ssb = np.zeros(num_threads)
        double[::1] t_tmp = np.zeros(num_threads)
        int thread_id
    
    with nogil, cython.boundscheck(False), cython.wraparound(False), cython.cdivision(True):
        for c in prange(cols, schedule='static', num_threads=num_threads):
            t = threadid() # or
            
            # subtract the cth column
            for j in range(rows):
                for n in range(nrndm):
                    i = ixs[c, n]
                    t_A[t, j, n] = e[j, i] - e[j, c]
                    
            
            #muA = A.mean(0)
            for j in range(rows):
                for n in range(nrndm):
                    i = ixs[c, n]
                    t_muA[t, n] += t_A[t, j, n]
            for n in range(nrndm):
                t_muA[t, n] = t_muA[t, n] / rows

            # A_mA = A - muA
            for j in range(rows):
                for n in range(nrndm):
                    t_A_mA[t, j, n] = t_A[t, j, n] - t_muA[t, n]

            # mub = b.mean()
            for j in range(rows):
                t_mub[t] += t_b[c, j]
            t_mub[t] = t_mub[t] / rows

            # b_mb = b - mub
            for j in range(rows):
                t_b_mb[t, j] = t_b[c, j] - t_mub[t]

            # ssA = (A_mA**2).sum(0)
            for j in range(rows):
                for n in range(nrndm):
                    t_ssA[t, n] += t_A_mA[t, j, n] * t_A_mA[t, j, n]
            for n in range(nrndm):
                t_ssA[t, n] = 1. / sqrt(t_ssA[t, n])

            # ssb = (b_mb**2).sum()
            for j in range(rows):
                t_ssb[t] += t_b_mb[t, j] * t_b_mb[t, j] # **2
            t_ssb[t] = 1. / sqrt(t_ssb[t])

            # np.dot(b_mb, A_mA)/(np.sqrt(ssA) * np.sqrt(ssb))
            for j in range(rows):
                t_tmp[t] = t_b_mb[t, j] * t_ssb[t]
                for n in range(nrndm):
                    i = ixs[c, n]
                    rm[c, i] += (t_A_mA[t, j, n] * t_ssA[t, n]) * t_tmp[t]
                    
            # Cleanup
            for n in range(nrndm):
                t_muA[t, n] = 0
            
            t_mub[t] = 0
                
            for n in range(nrndm):
                t_ssA[t, n] = 0
                
            t_ssb[t] = 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _colDeltaCorLog10partial(double[:, ::1] e,
                             double[:, ::1] d,
                             double[:, ::1] rm,
                             Py_ssize_t[:, ::1] ixs,
                             int num_threads,
                             double psc):
    cdef:
        Py_ssize_t rows = e.shape[0]
        Py_ssize_t cols = e.shape[1]
        Py_ssize_t nrndm = ixs.shape[1]
        Py_ssize_t i, j, c, t, n
        double[:, :, ::1] t_A = np.zeros((num_threads, e.shape[0], nrndm)) # np.tile(e, (3,1,1))
        double[:, ::1] t_b = np.array(d.T, order="C")
        double[:, ::1] t_out = np.zeros((num_threads, nrndm))
        double[:, ::1] t_muA = np.zeros((num_threads, nrndm))
        double[:, :, ::1] t_A_mA = np.zeros((num_threads, e.shape[0], nrndm))
        double[:, ::1] t_b_mb = np.zeros((num_threads, d.shape[0]))
        double[:, ::1] t_ssA = np.zeros((num_threads, nrndm))
        double[::1] t_mub = np.zeros(num_threads)
        double[::1] t_ssb = np.zeros(num_threads)
        double[::1] t_tmp = np.zeros(num_threads)
    
    with nogil, cython.boundscheck(False), cython.wraparound(False), cython.cdivision(True):
        for c in prange(cols, schedule='static', num_threads=num_threads):
            t = threadid() # or
            
            # subtract the cth column
            for j in range(rows):
                for n in range(nrndm):
                    i = ixs[c, n]
                    t_tmp[t] = e[j, i] - e[j, c]
                    if t_tmp[t] > 0:
                        t_A[t, j, n] = log10(fabs(t_tmp[t]) + psc)
                    else:
                        t_A[t, j, n] = -log10(fabs(t_tmp[t]) + psc)
                    
            
            #muA = A.mean(0)
            for j in range(rows):
                for n in range(nrndm):
                    i = ixs[c, n]
                    t_muA[t, n] += t_A[t, j, n]
            for n in range(nrndm):
                t_muA[t, n] = t_muA[t, n] / rows

            # A_mA = A - muA
            for j in range(rows):
                for n in range(nrndm):
                    t_A_mA[t, j, n] = t_A[t, j, n] - t_muA[t, n]

            # mub = b.mean()
            for j in range(rows):
                t_mub[t] += t_b[c, j]
            t_mub[t] = t_mub[t] / rows

            # b_mb = b - mub
            for j in range(rows):
                t_b_mb[t, j] = t_b[c, j] - t_mub[t]

            # ssA = (A_mA**2).sum(0)
            for j in range(rows):
                for n in range(nrndm):
                    t_ssA[t, n] += t_A_mA[t, j, n] * t_A_mA[t, j, n]
            for n in range(nrndm):
                t_ssA[t, n] = 1. / sqrt(t_ssA[t, n])

            # ssb = (b_mb**2).sum()
            for j in range(rows):
                t_ssb[t] += t_b_mb[t, j] * t_b_mb[t, j] # **2
            t_ssb[t] = 1. / sqrt(t_ssb[t])

            # np.dot(b_mb, A_mA)/(np.sqrt(ssA) * np.sqrt(ssb))
            for j in range(rows):
                t_tmp[t] = t_b_mb[t, j] * t_ssb[t]
                for n in range(nrndm):
                    i = ixs[c, n]
                    rm[c, i] += (t_A_mA[t, j, n] * t_ssA[t, n]) * t_tmp[t]
                    
            # Cleanup
            for n in range(nrndm):
                t_muA[t, n] = 0
            
            t_mub[t] = 0
                
            for n in range(nrndm):
                t_ssA[t, n] = 0
                
            t_ssb[t] = 0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _colDeltaCorSqrtpartial(double[:, ::1] e,
                             double[:, ::1] d,
                             double[:, ::1] rm,
                             Py_ssize_t[:, ::1] ixs,
                             int num_threads,
                             double psc):
    cdef:
        Py_ssize_t rows = e.shape[0]
        Py_ssize_t cols = e.shape[1]
        Py_ssize_t nrndm = ixs.shape[1]
        Py_ssize_t i, j, c, t, n
        double[:, :, ::1] t_A = np.zeros((num_threads, e.shape[0], nrndm)) # np.tile(e, (3,1,1))
        double[:, ::1] t_b = np.array(d.T, order="C")
        double[:, ::1] t_out = np.zeros((num_threads, nrndm))
        double[:, ::1] t_muA = np.zeros((num_threads, nrndm))
        double[:, :, ::1] t_A_mA = np.zeros((num_threads, e.shape[0], nrndm))
        double[:, ::1] t_b_mb = np.zeros((num_threads, d.shape[0]))
        double[:, ::1] t_ssA = np.zeros((num_threads, nrndm))
        double[::1] t_mub = np.zeros(num_threads)
        double[::1] t_ssb = np.zeros(num_threads)
        double[::1] t_tmp = np.zeros(num_threads)
    
    with nogil, cython.boundscheck(False), cython.wraparound(False), cython.cdivision(True):
        for c in prange(cols, schedule='static', num_threads=num_threads):
            t = threadid() # or
            
            # subtract the cth column
            for j in range(rows):
                for n in range(nrndm):
                    i = ixs[c, n]
                    t_tmp[t] = e[j, i] - e[j, c]
                    if t_tmp[t] > 0:
                        t_A[t, j, n] = sqrt(fabs(t_tmp[t]) + psc)
                    else:
                        t_A[t, j, n] = -sqrt(fabs(t_tmp[t]) + psc)
                    
            
            #muA = A.mean(0)
            for j in range(rows):
                for n in range(nrndm):
                    i = ixs[c, n]
                    t_muA[t, n] += t_A[t, j, n]
            for n in range(nrndm):
                t_muA[t, n] = t_muA[t, n] / rows

            # A_mA = A - muA
            for j in range(rows):
                for n in range(nrndm):
                    t_A_mA[t, j, n] = t_A[t, j, n] - t_muA[t, n]

            # mub = b.mean()
            for j in range(rows):
                t_mub[t] += t_b[c, j]
            t_mub[t] = t_mub[t] / rows

            # b_mb = b - mub
            for j in range(rows):
                t_b_mb[t, j] = t_b[c, j] - t_mub[t]

            # ssA = (A_mA**2).sum(0)
            for j in range(rows):
                for n in range(nrndm):
                    t_ssA[t, n] += t_A_mA[t, j, n] * t_A_mA[t, j, n]
            for n in range(nrndm):
                t_ssA[t, n] = 1. / sqrt(t_ssA[t, n])

            # ssb = (b_mb**2).sum()
            for j in range(rows):
                t_ssb[t] += t_b_mb[t, j] * t_b_mb[t, j] # **2
            t_ssb[t] = 1. / sqrt(t_ssb[t])

            # np.dot(b_mb, A_mA)/(np.sqrt(ssA) * np.sqrt(ssb))
            for j in range(rows):
                t_tmp[t] = t_b_mb[t, j] * t_ssb[t]
                for n in range(nrndm):
                    i = ixs[c, n]
                    rm[c, i] += (t_A_mA[t, j, n] * t_ssA[t, n]) * t_tmp[t]
                    
            # Cleanup
            for n in range(nrndm):
                t_muA[t, n] = 0
            
            t_mub[t] = 0
                
            for n in range(nrndm):
                t_ssA[t, n] = 0
                
            t_ssb[t] = 0