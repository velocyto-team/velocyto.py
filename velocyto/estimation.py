import numpy as np
import scipy.optimize
from scipy import sparse
import logging
from typing import *
from sklearn.neighbors import NearestNeighbors
from .speedboosted import _colDeltaCor, _colDeltaCorLog10, _colDeltaCorSqrt
from .speedboosted import _colDeltaCorpartial, _colDeltaCorLog10partial, _colDeltaCorSqrtpartial


def colDeltaCor(emat: np.ndarray, dmat: np.ndarray, threads: int=None) -> np.ndarray:
    """Calculate the correlation between the displacement (d[:,i])
    and the difference between a cell and every other (e - e[:, i])
    
    Parallel cython+OpenMP implemetation

    Arguments
    ---------
    emat: np.ndarray (ngenes, ncells)
        gene expression matrix
    dmat: np.ndarray (ngenes, ncells)
        gene velocity/displacement matrix
    threads: int
        number of parallel threads to use
    """
    import multiprocessing
    if threads is None:
        num_threads = int(multiprocessing.cpu_count() / 2)
    else:
        num_threads = max(threads, multiprocessing.cpu_count())
    out = np.zeros((emat.shape[1], emat.shape[1]))
    _colDeltaCor(emat, dmat, out, num_threads)
    return out


def colDeltaCorpartial(emat: np.ndarray, dmat: np.ndarray, ixs: np.ndarray, threads: int=None) -> np.ndarray:
    """Calculate the correlation between the displacement (d[:,i])
    and the difference between a cell and every other (e - e[:, i])
    
    Parallel cython+OpenMP implemetation

    Arguments
    ---------
    emat: np.ndarray (ngenes, ncells)
        gene expression matrix
    dmat: np.ndarray (ngenes, ncells)
        gene velocity/displacement matrix
    ixs: the neighborhood matrix (ncells, nneighbours)
        ixs[i, k] is the kth neighbour to the cell i
    threads: int
        number of parallel threads to use
    """
    import multiprocessing
    if threads is None:
        num_threads = int(multiprocessing.cpu_count() / 2)
    else:
        num_threads = max(threads, multiprocessing.cpu_count())
    out = np.zeros((emat.shape[1], emat.shape[1]))
    _colDeltaCorpartial(emat, dmat, out, ixs.astype("int32"), num_threads)
    return out


def colDeltaCorLog10(emat: np.ndarray, dmat: np.ndarray, threads: int=None, psc: float=1.0) -> np.ndarray:
    """Calculate the correlation between the displacement (d[:,i])
    and the difference between a cell and every other (e - e[:, i])
    
    Parallel cython+OpenMP implemetation

    Arguments
    ---------
    emat: np.ndarray (ngenes, ncells)
        gene expression matrix
    dmat: np.ndarray (ngenes, ncells)
        gene velocity/displacement matrix
    threads: int
        number of parallel threads to use
    """
    import multiprocessing
    if threads is None:
        num_threads = int(multiprocessing.cpu_count() / 2)
    else:
        num_threads = max(threads, multiprocessing.cpu_count())
    out = np.zeros((emat.shape[1], emat.shape[1]))
    _colDeltaCorLog10(emat, dmat, out, num_threads, psc)
    return out


def colDeltaCorLog10partial(emat: np.ndarray, dmat: np.ndarray, ixs: np.ndarray, threads: int=None, psc: float=1.0) -> np.ndarray:
    """Calculate the correlation between the displacement (d[:,i])
    and the difference between a cell and every other (e - e[:, i])
    
    Parallel cython+OpenMP implemetation

    Arguments
    ---------
    emat: np.ndarray (ngenes, ncells)
        gene expression matrix
    dmat: np.ndarray (ngenes, ncells)
        gene velocity/displacement matrix
    ixs: the neighborhood matrix (ncells, nneighbours)
        ixs[i, k] is the kth neighbour to the cell i
    threads: int
        number of parallel threads to use
    """
    import multiprocessing
    if threads is None:
        num_threads = int(multiprocessing.cpu_count() / 2)
    else:
        num_threads = max(threads, multiprocessing.cpu_count())
    out = np.zeros((emat.shape[1], emat.shape[1]))
    emat = np.require(emat, requirements="C")
    ixs = np.require(ixs, requirements="C").astype(np.intp)
    _colDeltaCorLog10partial(emat, dmat, out, ixs, num_threads, psc)
    return out


def colDeltaCorSqrt(emat: np.ndarray, dmat: np.ndarray, threads: int=None, psc: float=0.0) -> np.ndarray:
    """Calculate the correlation between the displacement (d[:,i])
    and the difference between a cell and every other (e - e[:, i])
    
    Parallel cython+OpenMP implemetation

    Arguments
    ---------
    emat: np.ndarray (ngenes, ncells)
        gene expression matrix
    dmat: np.ndarray (ngenes, ncells)
        gene velocity/displacement matrix
    threads: int
        number of parallel threads to use
    """
    import multiprocessing
    if threads is None:
        num_threads = int(multiprocessing.cpu_count() / 2)
    else:
        num_threads = max(threads, multiprocessing.cpu_count())
    out = np.zeros((emat.shape[1], emat.shape[1]))
    _colDeltaCorSqrt(emat, dmat, out, num_threads, psc)
    return out


def colDeltaCorSqrtpartial(emat: np.ndarray, dmat: np.ndarray, ixs: np.ndarray, threads: int=None, psc: float=0.0) -> np.ndarray:
    """Calculate the correlation between the displacement (d[:,i])
    and the difference between a cell and every other (e - e[:, i])
    
    Parallel cython+OpenMP implemetation

    Arguments
    ---------
    emat: np.ndarray (ngenes, ncells)
        gene expression matrix
    dmat: np.ndarray (ngenes, ncells)
        gene velocity/displacement matrix
    ixs: the neighborhood matrix (ncells, nneighbours)
        ixs[i, k] is the kth neighbour to the cell i
    threads: int
        number of parallel threads to use
    """
    import multiprocessing
    if threads is None:
        num_threads = int(multiprocessing.cpu_count() / 2)
    else:
        num_threads = max(threads, multiprocessing.cpu_count())
    out = np.zeros((emat.shape[1], emat.shape[1]))
    emat = np.require(emat, requirements="C")
    ixs = np.require(ixs, requirements="C").astype(np.intp)
    _colDeltaCorSqrtpartial(emat, dmat, out, ixs, num_threads, psc)
    return out


def _fit1_slope(y: np.ndarray, x: np.ndarray) -> float:
    """Simple function that fit a linear regression model without intercept
    """
    if not np.any(x):
        m = np.NAN  # It is definetelly not at steady state!!!
    elif not np.any(y):
        m = 0
    else:
        result, rnorm = scipy.optimize.nnls(x[:, None], y)  # Fastest but costrains result >= 0
        m = result[0]
        # Second fastest: m, _ = scipy.optimize.leastsq(lambda m: x*m - y, x0=(0,))
        # Third fastest: m = scipy.optimize.minimize_scalar(lambda m: np.sum((x*m - y)**2 )).x
        # Before I was doinf fastest: scipy.optimize.minimize_scalar(lambda m: np.sum((y - m * x)**2), bounds=(0, 3), method="bounded").x
        # Optionally one could clip m if high value make no sense
        # m = np.clip(m,0,3)
    return m


def _fit1_slope_weighted(y: np.ndarray, x: np.ndarray, w: np.ndarray, bounds: Tuple[float, float]=(0, 3)) -> float:
    """Simple function that fit a weighted linear regression model without intercept
    """
    if not np.any(x):
        m = np.NAN  # It is definetelly not at steady state!!!
    elif not np.any(y):
        m = 0
    else:
        m = scipy.optimize.minimize_scalar(lambda m: np.sum(w * (x * m - y)**2), bounds=(0, 3), method="bounded").x
    return m


def _fit1_slope_weighted_offset(y: np.ndarray, x: np.ndarray, w: np.ndarray, fixperc_q: bool=False, limit_gamma: bool=False) -> Tuple[float, float]:
    """Function that fits a weighted linear regression model with intercept
    with some adhoc
    """
    if not np.any(x):
        m = (np.NAN, 0)  # It is definetelly not at steady state!!!
    elif not np.any(y):
        m = (0, 0)
    else:
        if fixperc_q:
            m1 = np.percentile(y[x <= np.percentile(x, 1)], 50)
            m0 = scipy.optimize.minimize_scalar(lambda m: np.sum(w * (x * m - y + m1)**2), bounds=(0, 3), method="bounded").x
            m = (m0, m1)
        else:
            # m, _ = scipy.optimize.leastsq(lambda m: np.sqrt(w) * (-y + x * m[0] + m[1]), x0=(0, 0))  # This is probably faster but it can have negative slope
            # NOTE: The up_gamma is to deal with cases where consistently y > x. Those should have positive velocity everywhere
            if limit_gamma:
                if np.median(y) > np.median(x):
                    high_x = x > np.percentile(x, 90)
                    up_gamma = np.percentile(y[high_x], 10) / np.median(x[high_x])
                    up_gamma = np.maximum(1.5, up_gamma)
                else:
                    up_gamma = 1.5  # Just a bit more than 1
            else:
                up_gamma = 30
            up_q = 2 * np.sum(y * w) / np.sum(w)
            m = scipy.optimize.minimize(lambda m: np.sum(w * (-y + x * m[0] + m[1])**2),
                                        x0=(0.1, 1e-16), method="L-BFGS-B",
                                        bounds=[(1e-8, up_gamma), (0, up_q)]).x  # If speedup is needed either the gradient or numexpr could be used
    return m[0], m[1]


def _fit1_slope_offset(y: np.ndarray, x: np.ndarray, fixperc_q: bool=False) -> Tuple[float, float]:
    """Simple function that fit a linear regression model with intercept
    """
    if not np.any(x):
        m = (np.NAN, 0)  # It is definetelly not at steady state!!!
    elif not np.any(y):
        m = (0, 0)
    else:
        # result, rnorm = scipy.optimize.nnls(x[:, None], y)  # Fastest but costrains result >= 0
        # m = result[0]
        if fixperc_q:
            m1 = np.percentile(y[x <= np.percentile(x, 1)], 50)
            m0 = scipy.optimize.minimize_scalar(lambda m: np.sum((x * m - y + m1)**2), bounds=(0, 3), method="bounded").x
            m = (m0, m1)
        else:
            m, _ = scipy.optimize.leastsq(lambda m: -y + x * m[0] + m[1], x0=(0, 0))
        # Third fastest: m = scipy.optimize.minimize_scalar(lambda m: np.sum((x*m - y)**2 )).x
        # Before I was doinf fastest: scipy.optimize.minimize_scalar(lambda m: np.sum((y - m * x)**2), bounds=(0, 3), method="bounded").x
        # Optionally one could clip m if high value make no sense
        # m = np.clip(m,0,3)
    return m[0], m[1]


def fit_slope(Y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Loop through the genes and fits the slope

    Y: np.ndarray, shape=(genes, cells)
        the dependent variable (unspliced)
    X: np.ndarray, shape=(genes, cells)
        the independent variable (spliced)
    """
    # NOTE this could be easily parallelized
    slopes = np.fromiter((_fit1_slope(Y[i, :], X[i, :]) for i in range(Y.shape[0])),
                         dtype="float32",
                         count=Y.shape[0])
    return slopes


def fit_slope_offset(Y: np.ndarray, X: np.ndarray, fixperc_q: bool=False) -> Tuple[np.ndarray, np.ndarray]:
    """Loop through the genes and fits the slope

    Y: np.ndarray, shape=(genes, cells)
        the dependent variable (unspliced)
    X: np.ndarray, shape=(genes, cells)
        the independent variable (spliced)
    """
    # NOTE this could be easily parallelized
    slopes = np.zeros(Y.shape[0], dtype="float32")
    offsets = np.zeros(Y.shape[0], dtype="float32")
    for i in range(Y.shape[0]):
        m, q = _fit1_slope_offset(Y[i, :], X[i, :], fixperc_q)
        slopes[i] = m
        offsets[i] = q
    return slopes, offsets


def fit_slope_weighted(Y: np.ndarray, X: np.ndarray, W: np.ndarray, bounds: Tuple[float, float]=(0, 3)) -> np.ndarray:
    """Loop through the genes and fits the slope

    Y: np.ndarray, shape=(genes, cells)
        the dependent variable (unspliced)
    X: np.ndarray, shape=(genes, cells)
        the independent variable (spliced)
    W: np.ndarray, shape=(genes, cells)
        the weights that will scale the square residuals
    """
    # NOTE this could be easily parallelized
    slopes = np.fromiter((_fit1_slope_weighted(Y[i, :], X[i, :], W[i, :], bounds=bounds) for i in range(Y.shape[0])),
                         dtype="float32",
                         count=Y.shape[0])
    return slopes


def fit_slope_weighted_offset(Y: np.ndarray, X: np.ndarray, W: np.ndarray, fixperc_q: bool=False, return_R2: bool=True, limit_gamma: bool=False) -> Any:
    """Loop through the genes and fits the slope

    Y: np.ndarray, shape=(genes, cells)
        the dependent variable (unspliced)
    X: np.ndarray, shape=(genes, cells)
        the independent variable (spliced)
    """
    # NOTE this could be easily parallelized
    slopes = np.zeros(Y.shape[0], dtype="float32")
    offsets = np.zeros(Y.shape[0], dtype="float32")
    if return_R2:
        R2 = np.zeros(Y.shape[0], dtype="float32")
    for i in range(Y.shape[0]):
        m, q = _fit1_slope_weighted_offset(Y[i, :], X[i, :], W[i, :], fixperc_q, limit_gamma)
        slopes[i] = m
        offsets[i] = q
        if return_R2:
            # NOTE: the coefficient of determination is not weighted but the fit is
            SSres = np.sum((m * X[i, :] + q - Y[i, :])**2)
            SStot = np.sum((Y[i, :].mean() - Y[i, :])**2)
            R2[i] = 1 - (SSres / SStot)
    if return_R2:
        return slopes, offsets, R2
    return slopes, offsets


def clusters_stats(U: np.ndarray, S: np.ndarray,
                   clusters_uid: np.ndarray, cluster_ix: np.ndarray,
                   size_limit: int=40) -> Tuple[np.ndarray, np.ndarray]:
    """Calculate the averages per cluster
    
    If the cluster is too small (size<size_limit) the average of the toal is reported instead
    """
    U_avgs = np.zeros((S.shape[0], len(clusters_uid)))
    S_avgs = np.zeros((S.shape[0], len(clusters_uid)))
    avgU_div_avgS = np.zeros((S.shape[0], len(clusters_uid)))
    slopes_by_clust = np.zeros((S.shape[0], len(clusters_uid)))
    for i, uid in enumerate(clusters_uid):
        cluster_filter = cluster_ix == i
        n_cells = np.sum(cluster_filter)
        logging.info(f"Cluster: {uid} ({n_cells} cells)")
        if n_cells > size_limit:
            U_avgs[:, i], S_avgs[:, i] = U[:, cluster_filter].mean(1), S[:, cluster_filter].mean(1)
        else:
            U_avgs[:, i], S_avgs[:, i] = U.mean(1), S.mean(1)
            
    return U_avgs, S_avgs
