from copy import deepcopy
import warnings
import logging
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import norm as normal
import scipy.stats
from scipy import sparse
import matplotlib
import matplotlib.pyplot as plt
from sklearn.svm import SVR
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.neighbors import NearestNeighbors
from numba import jit
import loompy
from .neighbors import knn_distance_matrix, connectivity_to_weights, convolve_by_sparse_weights, BalancedKNN
from .estimation import fit_slope, fit_slope_offset, fit_slope_weighted, fit_slope_weighted_offset
from .estimation import clusters_stats
from .estimation import colDeltaCor, colDeltaCorSqrt, colDeltaCorLog10, colDeltaCorpartial, colDeltaCorSqrtpartial, colDeltaCorLog10partial
from .diffusion import Diffusion
from .serialization import dump_hdf5, load_hdf5
from typing import *


class VelocytoLoom:
    """A convenient object to store the data of a velocyto loom file.

    Data will be stored in memory

    Examples
    --------
    For usage examples consult the documentation

    Attributes
    ----------
    S: np.ndarray
        Expressed spliced molecules
    U: np.ndarray
        Unspliced molecule count
    A: np.ndarray
        Ambiguous molecule count
    ca: dict
        Column attributes of the loom file
    ra: dict
        Row attributes of the loom file
    loom_filepath: str
        The original path the loom files has been read from
    initial_cell_size: int
        The sum of spliced molecules
    initial_Ucell_size: int
        The sum of unspliced molecules

    """

    def __init__(self, loom_filepath: str) -> None:
        self.loom_filepath = loom_filepath
        ds = loompy.connect(self.loom_filepath)
        self.S = ds.layer["spliced"][:, :]
        self.U = ds.layer["unspliced"][:, :]
        self.A = ds.layer["ambiguous"][:, :]
        self.ca = dict(ds.col_attrs.items())
        self.ra = dict(ds.row_attrs.items())
        ds.close()

        self.initial_cell_size = self.S.sum(0)
        self.initial_Ucell_size = self.U.sum(0)

        try:
            if np.mean(self.ca["_Valid"]) < 1:
                logging.warn(f"fraction of _Valid cells is {np.mean(self.ca['_Valid'])} but all will be taken in consideration")
        except KeyError:
            pass
            # logging.debug("The file did not specify the _Valid column attribute")

    def to_hdf5(self, filename: str, **kwargs: Dict[str, Any]) -> None:
        """Serialize the VelocytoLoom object in its current state

        Arguments
        ---------
        filename:
            The name of the file that will be generated (the suffix hdf5 is suggested but not enforced)
        **kwargs:
            The function accepts the arguments of `dump_hdf5`

        Returns
        -------
        Nothing. It saves a file that can be used to recreate the object in another session.

        Note
        ----
        The object can be reloaded calling ``load_velocyto_hdf5(filename)``
        """
        dump_hdf5(self, filename, **kwargs)

    def plot_fractions(self, save2file: str=None) -> None:
        """Plots a barplot showing the abundance of spliced/unspliced molecules in the dataset

        Arguments
        ---------
        save2file: str (default: None)
            If not None specifies the file path to which plots get saved

        Returns
        -------
        Nothing, it plots a barplot
        """
        plt.figure(figsize=(3.2, 5))
        try:
            chips, chip_ix = np.unique(self.ca["SampleID"], return_inverse=1)
        except KeyError:
            chips, chip_ix = np.unique([i.split(":")[0] for i in self.ca["CellID"]], return_inverse=1)
        n = len(chips)
        for i in np.unique(chip_ix):
            tot_mol_cell_submatrixes = [X[:, chip_ix == i].sum(0) for X in [self.S, self.A, self.U]]
            total = np.sum(tot_mol_cell_submatrixes, 0)
            _mean = [np.mean(j / total) for j in tot_mol_cell_submatrixes]
            _std = [np.std(j / total) for j in tot_mol_cell_submatrixes]
            plt.ylabel("Fraction")
            plt.bar(np.linspace(-0.2, 0.2, n)[i] + np.arange(3), _mean, 0.5 / (n * 1.05), label=chips[i])
            plt.errorbar(np.linspace(-0.2, 0.2, n)[i] + np.arange(3), _mean, _std, c="k", fmt="none", lw=1, capsize=2)

            # Hide the right and top spines
            plt.gca().spines['right'].set_visible(False)
            plt.gca().spines['top'].set_visible(False)
            # Only show ticks on the left and bottom spines
            plt.gca().yaxis.set_ticks_position('left')
            plt.gca().xaxis.set_ticks_position('bottom')
            plt.gca().spines['left'].set_bounds(0, 0.8)
            plt.legend()
            
        plt.xticks(np.arange(3), ["spliced", "ambiguous", "unspliced"])
        plt.tight_layout()
        if save2file:
            plt.savefig(save2file, bbox_inches="tight")

    def filter_cells(self, bool_array: np.ndarray) -> None:
        """Filter cells using a boolean array.

        Arguments
        ---------
        bool_array: np.ndarray (size )
            array describing the cells to keep (True).

        Return
        ------
        Nothing but it removes some cells from S and U.
        """
        self.S, self.U, self.A = (X[:, bool_array] for X in (self.S, self.U, self.A))
        self.initial_cell_size = self.initial_cell_size[bool_array]
        self.initial_Ucell_size = self.initial_Ucell_size[bool_array]
        try:
            self.ts = self.ts[bool_array]  # type: np.ndarray
        except:
            pass
        try:
            self.size_factor = self.size_factor[bool_array]  # type: np.ndarray
        except:
            pass
        self.ca = {k: v[bool_array] for k, v in self.ca.items()}
        try:
            self.cluster_labels = self.cluster_labels[bool_array]  # type: np.ndarray
            self.colorandum = self.colorandum[bool_array, :]  # type: np.ndarray
        except AttributeError:
            pass

    def set_clusters(self, cluster_labels: np.ndarray, cluster_colors_dict: Dict[str, List[float]]=None, colormap: Any=None) -> None:
        """Utility function to set cluster labels, names and colormap

        Arguments
        ---------
        cluster_labels: np.ndarray
            A vector of strings containing the name of the cluster for each cells
        cluster_colors_dict: dict[str, List[float]]
            A mapping  cluster_name -> rgb_color_triplet for example "StemCell":[0.65,0.1,0.4]
        colormap:
            (optional)
            In alternative to cluster_colors_dict a colormap object (e.g. from matplotlib or similar callable) can be passed

        Returns
        -------
        Nothing, the attributes `cluster_labels, colorandum, cluster_ix, cluster_uid` are created.

        """
        self.cluster_labels = np.array(cluster_labels)
        if self.cluster_labels.dtype == "O":  # Fixes a bug when importing from pandas
            self.cluster_labels = self.cluster_labels.astype(np.string_)
        if cluster_colors_dict:
            self.colorandum = np.array([cluster_colors_dict[i] for i in cluster_labels])
            self.cluster_colors_dict = cluster_colors_dict
            self.colormap = None
        else:
            if colormap is None:
                self.colorandum = colormap_fun(self.cluster_ix)
                cluster_uid = self.cluster_uid
                self.cluster_colors_dict = {cluster_uid[i]: colormap_fun(i) for i in range(len(cluster_uid))}
            else:
                self.colormap = colormap
                self.colorandum = self.colormap(self.cluster_ix)
                cluster_uid = self.cluster_uid
                self.cluster_colors_dict = {cluster_uid[i]: self.colormap(i) for i in range(len(cluster_uid))}

    @property
    def cluster_uid(self) -> np.ndarray:
        clusters_uid = np.unique(self.cluster_labels)
        return clusters_uid

    @property
    def cluster_ix(self) -> np.ndarray:
        _, cluster_ix = np.unique(self.cluster_labels, return_inverse=True)
        return cluster_ix

    def score_cv_vs_mean(self, N: int=3000, min_expr_cells: int=2, max_expr_avg: float=20, min_expr_avg: int=0, svr_gamma: float=None,
                         winsorize: bool=False, winsor_perc: Tuple[float, float]=(1, 99.5), sort_inverse: bool=False, which: str="S",
                         plot: bool=False) -> np.ndarray:
        """Rank genes on the basis of a CV vs mean fit, it uses a nonparametric fit (Support Vector Regression)

        Arguments
        ---------
        N: int
            the number to select
        min_expr_cells: int, (default=2)
            minimum number of cells that express that gene for it to be considered in the fit
        min_expr_avg: int, (default=0)
            The minimum average accepted before discarding from the the gene as not expressed
        max_expr_avg: float, (default=20)
            The maximum average accepted before discarding from the the gene as house-keeping/outlier
        svr_gamma: float
            the gamma hyper-parameter of the SVR
        winsorize: bool
            Wether to winsorize the data for the cv vs mean model
        winsor_perc: tuple, default=(1, 99.5)
            the up and lower bound of the winsorization
        sort_inverse: bool, (default=False)
            if True it sorts genes from less noisy to more noisy (to use for size estimation not for feature selection)
        which: bool, (default="S")
            it performs the same cv_vs mean procedure on spliced "S" or unspliced "U" count
            "both" is NOT supported here because most often S the two procedure would have different parameters
            (notice that default parameters are good heuristics only for S)
        plot: bool, default=False
            whether to show a plot

        Returns
        -------
        Nothing but it creates the attributes
        cv_mean_score: np.ndarray
            How much the observed CV is higher than the one predicted by a noise model fit to the data
        cv_mean_selected: np.ndarray bool
            on the basis of the N parameter

        Note: genes excluded from the fit will have in the output the same score as the lowest scoring gene in the dataset.

        To perform the filtering use the method `filter_genes`
        """
        if which == "S":
            if winsorize:
                if min_expr_cells <= ((100 - winsor_perc[1]) * self.S.shape[1] * 0.01):
                    min_expr_cells = int(np.ceil((100 - winsor_perc[1]) * self.S.shape[0] * 0.01)) + 2
                    logging.debug(f"min_expr_cells is too low for winsorization with upper_perc ={winsor_perc[1]}, upgrading to min_expr_cells ={min_expr_cells}")
                    
            detected_bool = ((self.S > 0).sum(1) > min_expr_cells) & (self.S.mean(1) < max_expr_avg) & (self.S.mean(1) > min_expr_avg)
            Sf = self.S[detected_bool, :]
            if winsorize:
                down, up = np.percentile(Sf, winsor_perc, 1)
                Sfw = np.clip(Sf, down[:, None], up[:, None])
                mu = Sfw.mean(1)
                sigma = Sfw.std(1, ddof=1)
            else:
                mu = Sf.mean(1)
                sigma = Sf.std(1, ddof=1)

            cv = sigma / mu
            log_m = np.log2(mu)
            log_cv = np.log2(cv)

            if svr_gamma is None:
                svr_gamma = 150. / len(mu)
                logging.debug(f"svr_gamma set to {svr_gamma}")
            # Fit the Support Vector Regression
            clf = SVR(gamma=svr_gamma)
            clf.fit(log_m[:, None], log_cv)
            fitted_fun = clf.predict
            ff = fitted_fun(log_m[:, None])
            score = log_cv - ff
            if sort_inverse:
                score = - score
            nth_score = np.sort(score)[::-1][N]
            if plot:
                scatter_viz(log_m[score > nth_score], log_cv[score > nth_score], s=3, alpha=0.4, c="tab:red")
                scatter_viz(log_m[score <= nth_score], log_cv[score <= nth_score], s=3, alpha=0.4, c="tab:blue")
                mu_linspace = np.linspace(np.min(log_m), np.max(log_m))
                plt.plot(mu_linspace, fitted_fun(mu_linspace[:, None]), c="k")
                plt.xlabel("log2 mean S")
                plt.ylabel("log2 CV S")
            self.cv_mean_score = np.zeros(detected_bool.shape)
            self.cv_mean_score[~detected_bool] = np.min(score) - 1e-16
            self.cv_mean_score[detected_bool] = score
            self.cv_mean_selected = self.cv_mean_score >= nth_score
        else:
            if winsorize:
                if min_expr_cells <= ((100 - winsor_perc[1]) * self.U.shape[1] * 0.01):
                    min_expr_cells = int(np.ceil((100 - winsor_perc[1]) * self.U.shape[0] * 0.01)) + 2
                    logging.debug(f"min_expr_cells is too low for winsorization with upper_perc ={winsor_perc[1]}, upgrading to min_expr_cells ={min_expr_cells}")
                    
            detected_bool = ((self.U > 0).sum(1) > min_expr_cells) & (self.U.mean(1) < max_expr_avg) & (self.U.mean(1) > min_expr_avg)
            Uf = self.U[detected_bool, :]
            if winsorize:
                down, up = np.percentile(Uf, winsor_perc, 1)
                Ufw = np.clip(Uf, down[:, None], up[:, None])
                mu = Ufw.mean(1)
                sigma = Ufw.std(1, ddof=1)
            else:
                mu = Uf.mean(1)
                sigma = Uf.std(1, ddof=1)

            cv = sigma / mu
            log_m = np.log2(mu)
            log_cv = np.log2(cv)

            if svr_gamma is None:
                svr_gamma = 150. / len(mu)
                logging.debug(f"svr_gamma set to {svr_gamma}")
            # Fit the Support Vector Regression
            clf = SVR(gamma=svr_gamma)
            clf.fit(log_m[:, None], log_cv)
            fitted_fun = clf.predict
            ff = fitted_fun(log_m[:, None])
            score = log_cv - ff
            if sort_inverse:
                score = - score
            nth_score = np.sort(score)[::-1][N]
            if plot:
                scatter_viz(log_m[score > nth_score], log_cv[score > nth_score], s=3, alpha=0.4, c="tab:red")
                scatter_viz(log_m[score <= nth_score], log_cv[score <= nth_score], s=3, alpha=0.4, c="tab:blue")
                mu_linspace = np.linspace(np.min(log_m), np.max(log_m))
                plt.plot(mu_linspace, fitted_fun(mu_linspace[:, None]), c="k")
                plt.xlabel("log2 mean U")
                plt.ylabel("log2 CV U")
            self.Ucv_mean_score = np.zeros(detected_bool.shape)
            self.Ucv_mean_score[~detected_bool] = np.min(score) - 1e-16
            self.Ucv_mean_score[detected_bool] = score
            self.Ucv_mean_selected = self.Ucv_mean_score >= nth_score

    def robust_size_factor(self, pc: float=0.1, which: str="both") -> None:
        """Calculates a size factor in a similar way of Anders and Huber 2010

        Arguments
        --------
        pc: float, default=0.1
            The pseudocount to add to the expression before taking the log for the purpose of the size factor calculation
        which: str, default="both"
            For which counts estimate the normalization size factor. It can be "both", "S" or "U"

        Returns
        -------
        Nothing but it creates the attribute `self.size_factor` and `self.Usize_factor`
        normalization is self.S / self.size_factor and is performed by using `self.normalize(relative_size=self.size_factor)`

        Note
        ----
        Before running this method `score_cv_vs_mean` need to be run with sort_inverse=True, since only lowly variable genes are used for this size estimation
        """
        if which == "both":
            Y = np.log2(self.S[self.cv_mean_selected, :] + pc)
            Y_avg = Y.mean(1)
            self.size_factor: np.ndarray = np.median(2**(Y - Y_avg[:, None]), axis=0)
            self.size_factor = self.size_factor / np.mean(self.size_factor)

            Y = np.log2(self.U[self.Ucv_mean_selected, :] + pc)
            Y_avg = Y.mean(1)
            self.Usize_factor: np.ndarray = np.median(2**(Y - Y_avg[:, None]), axis=0)
            self.Usize_factor = self.Usize_factor / np.mean(self.Usize_factor)
        elif which == "S":
            Y = np.log2(self.S[self.cv_mean_selected, :] + pc)
            Y_avg = Y.mean(1)
            self.size_factor: np.ndarray = np.median(2**(Y - Y_avg[:, None]), axis=0)
            self.size_factor = self.size_factor / np.mean(self.size_factor)
        elif which == "U":
            Y = np.log2(self.U[self.Ucv_mean_selected, :] + pc)
            Y_avg = Y.mean(1)
            self.Usize_factor: np.ndarray = np.median(2**(Y - Y_avg[:, None]), axis=0)
            self.Usize_factor = self.Usize_factor / np.mean(self.Usize_factor)

    def score_cluster_expression(self, min_avg_U: float=0.02, min_avg_S: float=0.08) -> np.ndarray:
        """Prepare filtering genes on the basis of cluster-wise expression threshold

        Arguments
        ---------
        min_avg_U: float
            Include genes that have unspliced average bigger than `min_avg_U` in at least one of the clusters
        min_avg_S: float
            Include genes that have spliced average bigger than `min_avg_U` in at least one of the clusters
        Note: the two conditions are combined by and "&" logical operator

        Returns
        -------
        Nothing but it creates the attribute
        clu_avg_selected: np.ndarray bool
            The gene cluster that is selected
        To perform the filtering use the method `filter_genes`
        """
        self.U_avgs, self.S_avgs = clusters_stats(self.U, self.S, self.cluster_uid, self.cluster_ix, size_limit=40)
        self.clu_avg_selected = (self.U_avgs.max(1) > min_avg_U) & (self.S_avgs.max(1) > min_avg_S)

    def score_detection_levels(self, min_expr_counts: int= 50, min_cells_express: int= 20,
                               min_expr_counts_U: int= 0, min_cells_express_U: int= 0) -> np.ndarray:
        """Prepare basic filtering of genes on the basis of their detection levels

        Arguments
        ---------
        min_expr_counts: float
            The minimum number of spliced molecules detected considering all the cells
        min_cells_express: float
            The minimum number of cells that express spliced molecules of a gene
        min_expr_counts_U: float
            The minimum number of unspliced molecules detected considering all the cells
        min_cells_express_U: float
            The minimum number of cells that express unspliced molecules of a gene
        Note: the conditions are combined by and "&" logical operator

        Returns
        -------
        Nothing but an attribute self.detection_level_selected is created
        To perform the filtering by detection levels use the method `filter_genes`
        """
        # Some basic filtering
        S_sum = self.S.sum(1)
        S_ncells_express = (self.S > 0).sum(1)
        U_sum = self.U.sum(1)
        U_ncells_express = (self.U > 0).sum(1)
        filter_bool = (S_sum >= min_expr_counts) & (S_ncells_express >= min_cells_express) & (U_sum >= min_expr_counts_U) & (U_ncells_express >= min_cells_express_U)
        self.detection_level_selected = filter_bool

    def filter_genes(self, by_detection_levels: bool=False, by_cluster_expression: bool=False,
                     by_cv_vs_mean: bool=False, by_custom_array: Any=None, keep_unfiltered: bool=False) -> None:
        """Filter genes taking care that all the matrixes and all the connected annotation get filtered accordingly

        Attributes affected: .U, .S, .ra

        Arguments
        ---------
        by_detection_levels: bool, default=False
            filter genes by the score_detection_levels result

        by_cluster_expression: bool, default=False
            filter genes by the score_cluster_expression result

        by_cv_vs_mean: bool, default=False
            filter genes by the score_cluster_expression result

        by_custom_array, np.ndarray, default=None
            provide a boolean or index array

        keep_unfiltered: bool, default=False
            whether to create attributes self.S_prefilter, self.U_prefilter, self.ra_prefilter,
            (array will be made sparse to minimize memory footprint)
            or just overwrite the previous arrays

        Returns
        -------
        Nothing but it updates the self.S, self.U, self.ra attributes
        """
        assert np.any([by_detection_levels, by_cluster_expression,
                       by_cv_vs_mean, (type(by_custom_array) is np.ndarray)]), "At least one of the filtering methods needs to be True"
        tmp_filter = np.ones(self.S.shape[0], dtype=bool)
        if by_cluster_expression:
            assert hasattr(self, "clu_avg_selected"), "clu_avg_selected was not found"
            logging.debug("Filtering by cluster expression")
            tmp_filter = tmp_filter & self.clu_avg_selected
        if by_cv_vs_mean:
            assert hasattr(self, "cv_mean_selected"), "cv_mean_selected was not found"
            logging.debug("Filtering by cv vs mean")
            tmp_filter = tmp_filter & self.cv_mean_selected
        if by_detection_levels:
            assert hasattr(self, "detection_level_selected"), "detection_level_selected was not found"
            logging.debug("Filtering by detection level")
            tmp_filter = tmp_filter & self.detection_level_selected
        if type(by_custom_array) is np.ndarray:
            if by_custom_array.dtype == bool:
                logging.debug("Filtering by custom boolean array")
                tmp_filter = tmp_filter & by_custom_array
            elif by_custom_array.dtype == int:
                logging.debug("Filtering by custom index array")
                bool_negative = ~np.in1d(np.arange(len(tmp_filter)), by_custom_array)
                tmp_filter[bool_negative] = False

        if keep_unfiltered:
            if hasattr(self, "U_prefilter"):
                logging.debug("Attributes *_prefilter are already present and were overwritten")
            self.U_prefilter = sparse.csr_matrix(self.U)
            self.S_prefilter = sparse.csr_matrix(self.S)
            self.ra_prefilter = deepcopy(self.ra)
            
        self.U = self.U[tmp_filter, :]
        self.S = self.S[tmp_filter, :]
        self.ra = {k: v[tmp_filter] for k, v in self.ra.items()}

    def custom_filter_attributes(self, attr_names: List[str], bool_filter: np.ndarray) -> None:
        """Filter attributes given a boolean array. attr_names can be dictionaries or numpy arrays
        
        Arguments
        ---------
        attr_names: List[str]
            a list of the attributes to be modified. The can be
            1d arrays, dictionary of 1d arrays, ndarrays, will be filtered by axis=0
            if .T is specified by axis=-1
        bool_filter:
            the boolean filter to be applied

        Returns
        -------
        Nothing it filters the specified attributes
        """
        transpose_flag = False
        for attr in attr_names:
            if attr[-2:] == ".T":
                obj = getattr(self, attr[:-2])
                transpose_flag = True
            else:
                obj = getattr(self, attr)
                transpose_flag = False
            if type(obj) is dict:
                setattr(self, attr, {k: v[bool_filter] for k, v in obj.items()})
            elif type(obj) is np.ndarray:
                if len(obj.shape) > 1:
                    if transpose_flag:
                        setattr(self, attr, obj[..., bool_filter])
                    else:
                        setattr(self, attr, obj[bool_filter, :])
                else:
                    setattr(self, attr, obj[bool_filter])
            else:
                raise NotImplementedError(f"The filtering of an object of type {type(obj)} is not defined")

    def _normalize_S(self, size: bool=True, log: bool=True, pcount: float=1, relative_size: Any=None, target_size: Any=None) -> np.ndarray:
        """Internal function for the spliced molecule filtering. The `normalize` method should be used as a standard interface"""
        if size:
            if type(relative_size) is np.ndarray:
                self.cell_size = relative_size
            else:
                self.cell_size = self.S.sum(0)
            if target_size is None:
                self.avg_size = self.cell_size.mean()
            else:
                self.avg_size = target_size
            self.norm_factor = self.avg_size / self.cell_size
        else:
            self.norm_factor = 1
        self.S_sz = self.norm_factor * self.S
        if log:
            self.S_norm = np.log2(self.S_sz + pcount)  # np.sqrt(S_sz )# np.log2(S_sz + 1)

    def _normalize_U(self, size: bool=True, log: bool=True, pcount: float=1, use_S_size: bool=False, relative_size: np.ndarray=None, target_size: Any=None) -> np.ndarray:
        """Internal function for the unspliced molecule filtering. The `normalize` method should be used as a standard interface"""
        if size:
            if use_S_size:
                if hasattr(self, "cell_size"):
                    cell_size = self.cell_size
                else:
                    cell_size = self.S.sum(0)
            elif type(relative_size) is np.ndarray:
                cell_size = relative_size
            else:
                cell_size = self.U.sum(0)
            self.Ucell_size = cell_size
            if target_size is None:
                avg_size = cell_size.mean()
            else:
                avg_size = target_size
            self.Uavg_size = avg_size
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                norm_factor = avg_size / cell_size
        else:
            norm_factor = 1
        self.Unorm_factor = norm_factor
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.U_sz = norm_factor * self.U
        self.U_sz[~np.isfinite(self.U_sz)] = 0  # it happened only once but it is here as a precaution
        if log:
            self.U_norm = np.log2(self.U_sz + pcount)  # np.sqrt(S_sz )# np.log2(S_sz + 1)

    def _normalize_Sx(self, size: bool=True, log: bool=True, pcount: float=1, relative_size: Any=None, target_size: Any=None) -> np.ndarray:
        """Internal function for the smoothed spliced molecule filtering. The `normalize` method should be used as a standard interface"""
        if size:
            if relative_size:
                self.xcell_size = relative_size
            else:
                self.xcell_size = self.Sx.sum(0)
            if target_size is None:
                self.xavg_size = self.xcell_size.mean()
            else:
                self.xavg_size = target_size
            self.xnorm_factor = self.xavg_size / self.xcell_size
        else:
            self.xnorm_factor = 1
        self.Sx_sz = self.xnorm_factor * self.Sx
        if log:
            self.Sx_norm = np.log2(self.Sx_sz + pcount)  # np.sqrt(S_sz )# np.log2(S_sz + 1)

    def _normalize_Ux(self, size: bool=True, log: bool=True, pcount: float=1, use_Sx_size: bool=False, relative_size: Any=None, target_size: Any=None) -> np.ndarray:
        """Internal function for the smoothed unspliced molecule filtering. The `normalize` method should be used as a standard interface"""
        if size:
            if use_Sx_size:
                if hasattr(self, "cell_size"):
                    cell_size = self.xcell_size
                else:
                    cell_size = self.Sx.sum(0)
            elif type(relative_size) is np.ndarray:
                cell_size = relative_size
            else:
                cell_size = self.Ux.sum(0)
            self.xUcell_size = cell_size
            if target_size is None:
                avg_size = cell_size.mean()
            else:
                avg_size = target_size
            self.xUavg_size = avg_size
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                norm_factor = avg_size / cell_size
        else:
            norm_factor = 1
        self.xUnorm_factor = norm_factor
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.Ux_sz = norm_factor * self.Ux
        self.Ux_sz[~np.isfinite(self.Ux_sz)] = 0  # it happened only once but it is here as a precaution
        if log:
            self.Ux_norm = np.log2(self.Ux_sz + pcount)  # np.sqrt(S_sz )# np.log2(S_sz + 1)

    def normalize(self, which: str="both", size: bool=True, log: bool=True, pcount: float=1,
                  relative_size: np.ndarray=None, use_S_size_for_U: bool=False, target_size: Tuple[float, float]=(None, None)) -> None:
        """Normalization interface

        Arguments
        ---------
        which: either 'both', 'S', 'U', "imputed", "Sx", "Ux"
            which attributes to normalize.
            "both" corresponds to "S" and "U"
            "imputed" corresponds to "Sx" and "Ux"
        size: bool
            perform size normalization
        log: bool
            perform log normalization (if size==True, this comes after the size normalization)
        pcount: int, default: 1
            The extra count added when logging (log2)
        relative_size: np.ndarray, default=None
            if None it calculate the sums the molecules per cell (self.S.sum(0))
            if an array is provided it is used for the normalization
        use_S_size_for_U: bool
            U is size normalized using the sum of molecules of S
        target_size: float or Tuple[float, float] (depending if the which parameter implies 1 or more normalizations)
            the size of the cells after normalization will be set to.
            If tuple the order is (S, U) or (Sx, Ux)
            If None the target size is the average of the cell sizes
        Returns
        -------
        Nothing but creates the attributes `U_norm`, `U_sz` and `S_norm`, "S_sz"
        or `Ux_norm`, `Ux_sz` and `Sx_norm`, "Sx_sz"
        """
        if which == "both":
            self._normalize_S(size=size, log=log, pcount=pcount, relative_size=relative_size, target_size=target_size[0])
            self._normalize_U(size=size, log=log, pcount=pcount, use_S_size=use_S_size_for_U, relative_size=relative_size, target_size=target_size[1])
        if "S" == which:
            self._normalize_S(size=size, log=log, pcount=pcount, relative_size=relative_size, target_size=target_size[0])
        if "U" == which:
            self._normalize_U(size=size, log=log, pcount=pcount, use_S_size=use_S_size_for_U, relative_size=relative_size, target_size=target_size[1])
        if which == "imputed":
            self._normalize_Sx(size=size, log=log, pcount=pcount, relative_size=relative_size, target_size=target_size[0])
            self._normalize_Ux(size=size, log=log, pcount=pcount, use_Sx_size=use_S_size_for_U, relative_size=relative_size, target_size=target_size[1])
        if "Sx" == which:
            self._normalize_Sx(size=size, log=log, pcount=pcount, relative_size=relative_size, target_size=target_size[0])
        if "Ux" == which:
            self._normalize_Ux(size=size, log=log, pcount=pcount, use_Sx_size=use_S_size_for_U, relative_size=relative_size, target_size=target_size[1])

    def perform_PCA(self, which: str="S_norm", n_components: int=None, div_by_std: bool=False) -> None:
        """Perform PCA (cells as samples)

        Arguments
        ---------
        which: str, default="S_norm"
            The name of the attribute to use for the calculation (e.g. S_norm or Sx_norm)
        n_components: int, default=None
            Number of components to keep. If None all the components will be kept.
        div_by_std: bool, default=False
            Wether to divide by standard deviation

        Returns
        -------
        Returns nothing but it creates the attributes:
        pca: np.ndarray
            a numpy array of shape (cells, npcs)

        """
        X = getattr(self, which)
        self.pca = PCA(n_components=n_components)
        if div_by_std:
            self.pcs = self.pca.fit_transform(X.T / X.std(0))
        else:
            self.pcs = self.pca.fit_transform(X.T)

    def normalize_by_total(self, min_perc_U: float=0.5, plot: bool=False, skip_low_U_pop: bool=True, same_size_UnS: bool=False) -> None:
        """Normalize the cells using the (initial) total molecules as size estimate

        Arguments
        ---------
        min_perc_U: float
            the percentile to use as a minimum value allowed for the size normalization
        plot: bool, default=False
            whether
        skip_low_U_pop: bool, default=True
            population with very low unspliced will not be multiplied by the scaling factor to avoid predicting very strong
            velocity just as a consequence of low detection
        same_size_UnS: bool, default=False
            Each cell is set tot have the same total number of spliced and unspliced molecules

        Returns
        -------
        Returns nothing but it creates the attributes:
        small_U_pop: np.ndarray
            Cells with extremely low unspliced count

        """
        target_cell_size = np.median(self.initial_cell_size)
        min_Ucell_size = np.percentile(self.initial_Ucell_size, min_perc_U)
        if min_Ucell_size < 2:
            raise ValueError(f"min_perc_U={min_perc_U} corresponds to total Unspliced of 1 molecule of less. Please choose higher value or filter our these cell")
        bool_f = self.initial_Ucell_size < min_Ucell_size
        self.small_U_pop = bool_f
        if same_size_UnS:
            target_Ucell_size = target_cell_size  # 0.15 * target_cell_size
        else:
            target_Ucell_size = np.median(self.initial_Ucell_size[~self.small_U_pop])  # 0.15 * target_cell_size

        if plot:
            plt.figure(None, (12, 6))
            plt.subplot(121)
            
            plt.scatter(self.initial_cell_size, self.initial_Ucell_size, s=3, alpha=0.1)
            plt.xlabel("total spliced")
            plt.ylabel("total unspliced")
            plt.scatter(self.initial_cell_size[bool_f], self.initial_Ucell_size[bool_f], s=3, alpha=0.1)
            plt.subplot(122)
            plt.scatter(np.log2(self.initial_cell_size), np.log2(self.initial_Ucell_size), s=7, alpha=0.3)
            plt.scatter(np.log2(self.initial_cell_size)[bool_f], np.log2(self.initial_Ucell_size)[bool_f], s=7, alpha=0.3)
            plt.xlabel("log total spliced")
            plt.ylabel("log total unspliced")

        self._normalize_S(relative_size=self.initial_cell_size,
                          target_size=target_cell_size)
        if skip_low_U_pop:
            self._normalize_U(relative_size=np.clip(self.initial_Ucell_size, min_Ucell_size, None),
                              target_size=target_Ucell_size)
        else:
            self._normalize_U(relative_size=self.initial_Ucell_size,
                              target_size=target_Ucell_size)

    def normalize_by_size_factor(self, min_perc_U: float=0.5, plot: bool=False, skip_low_U_pop: bool=True, same_size_UnS: bool=False) -> None:
        """Normalize the cells using the (initial) size_factor

        Arguments
        ---------
        min_perc_U: float
            the percentile to use as a minimum value allowed for the size normalization
        plot: bool, default=False
            whether
        skip_low_U_pop: bool, default=True
            population with very low unspliced will not be multiplied by the scaling factor to avoid predicting very strong
            velocity just as a consequence of low detection
        same_size_UnS: bool, default=False
            Each cell is set tot have the same total number of spliced and unspliced molecules

        Returns
        -------
        Returns nothing but it creates the attributes:
        small_U_pop: np.ndarray
            Cells with extremely low unspliced count

        """
        cell_size = self.S.sum(0)
        Ucell_size = self.U.sum(0)
        target_cell_size = np.median(cell_size)
        min_Ucell_size = np.percentile(Ucell_size, min_perc_U)
        if min_Ucell_size < 2:
            raise ValueError(f"min_perc_U={min_perc_U} corresponds to total Unspliced of 1 molecule of less. Please choose higher value or filter our these cell")
        bool_f = Ucell_size < min_Ucell_size
        self.small_U_pop = bool_f
        if same_size_UnS:
            target_Ucell_size = target_cell_size  # 0.15 * target_cell_size
        else:
            target_Ucell_size = np.median(Ucell_size[~self.small_U_pop])
            
        if plot:
            plt.figure(None, (12, 6))
            plt.subplot(121)
            plt.scatter(cell_size, Ucell_size, s=3, alpha=0.1)
            plt.xlabel("S cell_size")
            plt.ylabel("U cell_size")
            plt.scatter(cell_size[bool_f], Ucell_size[bool_f], s=3, alpha=0.1)
            plt.subplot(122)
            plt.scatter(np.log2(cell_size), np.log2(Ucell_size), s=7, alpha=0.3)
            plt.scatter(np.log2(cell_size)[bool_f], np.log2(Ucell_size)[bool_f], s=7, alpha=0.3)
            plt.xlabel("log S cell_size")
            plt.ylabel("log U cell_size")

        self._normalize_S(relative_size=self.size_factor,
                          target_size=target_cell_size)
        if skip_low_U_pop:
            self._normalize_U(relative_size=np.clip(self.initial_Ucell_size, min_Ucell_size, None),
                              target_size=target_Ucell_size)
        else:
            self._normalize_U(relative_size=self.initial_Ucell_size,
                              target_size=target_Ucell_size)

    def adjust_totS_totU(self, skip_low_U_pop: bool=True, normalize_total: bool=False,
                         fit_with_low_U: bool=True,
                         svr_C: float=100, svr_gamma: float=1e-6, plot: bool=False) -> None:
        """Adjust the spliced count on the base of the relation S_sz_tot and U_sz_tot

        Arguments
        ---------
        skip_low_U_pop: bool, default=True
            Do not normalize the low unspliced molecules cell population to avoid overinflated values
        normalize_total: bool, default=False
            If this is True the function results in a normalization by median of both U and S.
            NOTE: Legacy compatibility, I might want to split this into a different function.
        fit_with_low_U: bool, default=True
            Wether to consider the low_U population for the fit
        svr_C: float
            The C parameter of scikit-learn Support Vector Regression
        svr_gamma: float
            The gamma parameter of scikit-learn Support Vector Regression
        plot: bool
            Whether to plot the results of the fit

        Returns
        -------
        Nothing but it modifies the attributes:
        U_sz: np.ndarray
        """

        svr = SVR(C=svr_C, kernel="rbf", gamma=svr_gamma)
        X, y = self.S_sz.sum(0), self.U_sz.sum(0)
        if fit_with_low_U:
            svr.fit(X[:, None], y)
            predicted = svr.predict(X[:, None])
        else:
            svr.fit(X[~self.small_U_pop, None], y[~self.small_U_pop])
            predicted = np.copy(y)
            predicted[~self.small_U_pop] = svr.predict(X[~self.small_U_pop, None])
            
        adj_factor = predicted / y
        adj_factor[~np.isfinite(adj_factor)] = 1
        if skip_low_U_pop:
            self.U_sz[:, ~self.small_U_pop] = self.U_sz[:, ~self.small_U_pop] * adj_factor[~self.small_U_pop]
        else:
            self.U_sz = self.U_sz * adj_factor

        if normalize_total:
            self.normalize_median(which="renormalize", skip_low_U_pop=skip_low_U_pop)

        if plot:
            plt.figure(None, (8, 8))
            plt.scatter(X, y, s=3, alpha=0.1)
            plt.scatter(X, predicted, c="k", s=5, alpha=0.1)

    def normalize_median(self, which: str="imputed", skip_low_U_pop: bool=True) -> None:
        """Normalize cell size to the median, for both S and U.

        Arguments
        ---------
        which: str, default="imputed"
            "imputed" or "renormalized"
        skip_low_U_pop: bool=True
            Whether to skip the low U population defined in normalize_by_total

        Returns
        -------
        Nothing but it modifies the attributes:
        S_sz: np.ndarray
        U_sz: np.ndarray
        or
        Sx_sz: np.ndarray
        Ux_sz: np.ndarray
        """
        if not hasattr(self, "small_U_pop") and skip_low_U_pop:
            self.small_U_pop = np.zeros(self.U_sz.shape[1], dtype=bool)
            logging.warning("object does not have the attribute `small_U_pop`, so all the unspliced will be normalized by relative size, this might cause the overinflation the unspliced counts of cells where only few unspliced molecules were detected")
        if which == "renormalize":
            self.S_sz = self.S_sz * (np.median(self.S_sz.sum(0)) / self.S_sz.sum(0))
            if skip_low_U_pop:
                self.U_sz[:, ~self.small_U_pop] = self.U_sz[:, ~self.small_U_pop] * (np.median(self.U_sz[:, ~self.small_U_pop].sum(0)) / self.U_sz[:, ~self.small_U_pop].sum(0))
            else:
                self.U_sz = self.U_sz * (np.median(self.U_sz.sum(0)) / self.U_sz.sum(0))
        elif which == "imputed":
            self.Sx_sz = self.Sx * (np.median(self.Sx.sum(0)) / self.Sx.sum(0))
            
            if skip_low_U_pop:
                self.Ux_sz = np.copy(self.Ux)
                self.Ux_sz[:, ~self.small_U_pop] = self.Ux[:, ~self.small_U_pop] * (np.median(self.Ux[:, ~self.small_U_pop].sum(0)) / self.Ux[:, ~self.small_U_pop].sum(0))
            else:
                self.Ux_sz = self.Ux * (np.median(self.Ux.sum(0)) / self.Ux.sum(0))

    def plot_pca(self, dim: List[int]=[0, 1, 2], elev: float=60, azim: float=-140) -> None:
        """Plot 3d PCA
        """
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(self.pcs[:, dim[0]],
                   self.pcs[:, dim[1]],
                   self.pcs[:, dim[2]],
                   c=self.colorandum)
        ax.view_init(elev=elev, azim=azim)

    def _perform_PCA_imputed(self, n_components: int=None) -> None:
        """Simply performs PCA of `Sx_norm` and save the result as  `pcax`"""
        self.pcax = PCA(n_components=n_components)
        self.pcsx = self.pcax.fit_transform(self.Sx_norm.T)

    def _plot_pca_imputed(self, dim: List[int]=[0, 1, 2], elev: float=60, azim: float=-140) -> None:
        """Plot 3d PCA of the smoothed data
        """
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(self.pcsx[:, dim[0]],
                   self.pcsx[:, dim[1]],
                   self.pcsx[:, dim[2]],
                   c=self.colorandum)
        ax.view_init(elev=elev, azim=azim)

    def knn_imputation(self, k: int=None, pca_space: float=True, metric: str="euclidean", diag: float=1,
                       n_pca_dims: int=None, maximum: bool=False, size_norm: bool=True,
                       balanced: bool=False, b_sight: int=None, b_maxl: int=None,
                       group_constraint: Union[str, np.ndarray]=None, n_jobs: int=8) -> None:
        """Performs k-nn smoothing of the data matrix

        Arguments
        ---------
        k: int
            number of neighbors. If None the default it is chosen to be `0.025 * Ncells`
        pca_space: bool, default=True
            if True the knn will be performed in PCA space (`pcs`)
            otherwise it will use log2 size normalized data  (`S_norm`)
        metric: str
            "euclidean" or "correlation"
        diag: int, default=1
            before smoothing this value is substituted in the diagonal of the knn contiguity matrix
            Resulting in a reduction of the smoothing effect.
            E.g. if diag=8 and k=10 value of Si = (8 * S_i + sum(S_n, with n in 5nn of i)) / (8+5)
        maximum: bool, default=False
            If True the maximum value of the smoothing and the original matrix entry is taken.
        n_pca_dims: int, default=None
            number of pca to use for the knn distance metric. If None all pcs will be used. (used only if pca_space == True)
        balanced: bool
            whether to use BalancedKNN version
        b_sight: int
            the sight parameter of BalancedKNN (used only if balanced == True)
        b_maxl: int
            the maxl parameter of BalancedKNN (used only if balanced == True)
        group_constraint: str or np.ndarray[int]:
            currently implemented only for balanced = True
            if "clusters" the the clusters will be used as a constraint so that cells of different clusters cannot be neighbors
            if an array of integers of shape vlm.S.shape[1] it will be interpreted as labels of the groups
        n_jobs: int, default 8
            number of parallel jobs in knn calculation

        Returns
        -------
        Nothing but it creates the attributes:
        knn: scipy.sparse.csr_matrix
            knn contiguity matrix
        knn_smoothing_w: scipy.sparse.lil_matrix
            the weights used for the smoothing
        Sx: np.ndarray
            smoothed spliced
        Ux: np.ndarray
            smoothed unspliced
        
        """
        N = self.S.shape[1]
        if k is None:
            k = int(N * 0.025)
        if b_sight is None and balanced:
            b_sight = np.maximum(int(k * 8), N - 1)
        if b_maxl is None and balanced:
            b_maxl = np.maximum(int(k * 4), N - 1)
        if pca_space:
            space = self.pcs[:, :n_pca_dims]
        else:
            space = self.S_norm.T
        if balanced:
            if group_constraint is not None:
                if isinstance(group_constraint, str) and group_constraint == "clusters":
                    constraint = np.array(self.cluster_ix)
                bknn = BalancedKNN(k=k, sight_k=b_sight, maxl=b_maxl, metric=metric, constraint=constraint, mode="distance", n_jobs=n_jobs)
            else:
                bknn = BalancedKNN(k=k, sight_k=b_sight, maxl=b_maxl, metric=metric, mode="distance", n_jobs=n_jobs)
            bknn.fit(space)
            self.knn = bknn.kneighbors_graph(mode="distance")
        else:
            if group_constraint is not None:
                raise ValueError("group_constraint is currently supported only if the argument balanced is set to True")
            self.knn = knn_distance_matrix(space, metric=metric, k=k, mode="distance", n_jobs=n_jobs)
        connectivity = (self.knn > 0).astype(float)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # SparseEfficiencyWarning: Changing the sparsity structure of a csr_matrix is expensive. lil_matrix is more efficient.
            connectivity.setdiag(diag)
        self.knn_smoothing_w = connectivity_to_weights(connectivity)
        if size_norm:
            self.Sx = convolve_by_sparse_weights(self.S_sz, self.knn_smoothing_w)
            self.Ux = convolve_by_sparse_weights(self.U_sz, self.knn_smoothing_w)
        else:
            self.Sx = convolve_by_sparse_weights(self.S, self.knn_smoothing_w)
            self.Ux = convolve_by_sparse_weights(self.U, self.knn_smoothing_w)
        if maximum:
            self.Sx = np.maximum(self.S_sz, self.Sx)
            self.Ux = np.maximum(self.U_sz, self.Ux)

        # Make a differently named varaible for backwards compatibility
        self.Sx_sz = np.copy(self.Sx)
        self.Ux_sz = np.copy(self.Ux)

    def knn_imputation_precomputed(self, knn_smoothing_w: sparse.lil_matrix, maximum: bool=False) -> None:
        """Performs k-nn imputation (like `.knn_imputation()`) but with a precomputed weight matrix
        
        Arguments
        ---------
        knn_smoothing_w: sparse.lil_matrix
            the sparse matrix to be convolved with self.S_sz and self.U_sz
            This should be the result of something like:
            connectivity.setdiag(diagonal_value)
            knn_smoothing_w = connectivity_to_weights(connectivity)
        maximum: bool, default=False
            whether to take the maximum value of the smoothing and the original matrix
        
        Returns
        -------
        Nothing but it creates the attributes:
        Sx: np.ndarray
            smoothed spliced
        Ux: np.ndarray
            smoothed unspliced
        """
        self.Sx = convolve_by_sparse_weights(self.S_sz, knn_smoothing_w)
        self.Ux = convolve_by_sparse_weights(self.U_sz, knn_smoothing_w)
        if maximum:
            self.Sx = np.maximum(self.S_sz, self.Sx)
            self.Ux = np.maximum(self.U_sz, self.Ux)

        self.Sx_sz = np.copy(self.Sx)
        self.Ux_sz = np.copy(self.Ux)

    def gene_knn_imputation(self, k: int=15, pca_space: float=False, metric: str="correlation", diag: float=1,
                            scale_weights: bool=True, balanced: bool=True, b_sight: int=100, b_maxl: int=18,
                            n_jobs: int=8) -> None:
        """Performs genes k-nn smoothing of the genes

        Arguments
        ---------
        k: int, default=15
            number of neighbors
        pca_space: bool, default=False
            if True the knn will be performed in PCA space (`pcs`)
            otherwise it will use log2 size normalized data  (`S_norm`)
        metric: str, default="correlation"
            "euclidean" or "correlation"
        diag: int, default=1
            before smoothing this value is substituted in the diagonal of the knn contiguity matrix
            Resulting in a reduction of the smoothing effect
            E.g. if diag=8 and k=10 value of Si = (8 * S_i + sum(S_n, with n in 5nn of i)) / (8+5)
        scale_weights: bool, default=True
            whether to scale weights by gene total expression/yield
        balanced: bool, default=True
            whether to use BalancedKNN version
        b_sight: int, default=100
            the sight parameter of BalancedKNN (used only if balanced == True)
        b_maxl: int, default=18
            the maxl parameter of BalancedKNN (used only if balanced == True)
        n_jobs: int, default=8
            number of parallel jobs in knn calculation

        Returns
        -------
        Nothing but it creates the attributes:
        gknn: scipy.sparse.csr_matrix
            genes knn contiguity matrix
        gknn_smoothing_w: scipy.sparse.lil_matrix
            the weights used for the smoothing of the genes
        Sx: np.ndarray
            smoothed spliced
        Ux: np.ndarray
            smoothed unspliced
        
        """
        if pca_space:
            assert NotImplementedError
        else:
            space = self.Sx_sz  # imputed size normalized counts
        if balanced:
            bknn = BalancedKNN(k=k, sight_k=b_sight, maxl=b_maxl, mode="distance", metric=metric, n_jobs=n_jobs)
            bknn.fit(space)
            self.gknn = bknn.kneighbors_graph(mode="distance")
        else:
            self.gknn = knn_distance_matrix(space, metric=metric, k=k, mode="distance", n_jobs=n_jobs)
        connectivity = (self.knn > 0).astype(float)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            connectivity.setdiag(diag)
        self.gknn_smoothing_w = connectivity_to_weights(connectivity).tocsr()

        if scale_weights:
            genes_total = space.sum(1)
            self.gknn_smoothing_w = scale_to_match_median(self.gknn_smoothing_w, genes_total)
        # NOTE This might be not computationally efficient after transpose, maybe better to use csc for the genes
        self.Sx_sz = convolve_by_sparse_weights(self.Sx_sz.T, self.gknn_smoothing_w).T
        self.Ux_sz = convolve_by_sparse_weights(self.Ux_sz.T, self.gknn_smoothing_w).T
        
    def fit_gammas(self, steady_state_bool: np.ndarray=None, use_imputed_data: bool=True, use_size_norm: bool=True,
                   fit_offset: bool=True, fixperc_q: bool=False, weighted: bool=True, weights: np.ndarray = "maxmin_diag",
                   limit_gamma: bool=False, maxmin_perc: List[float]=[2, 98], maxmin_weighted_pow: float=15) -> None:
        """Fit gamma using spliced and unspliced data

        Arguments
        ---------
        steady_state_bool: np.ndarray, default=None
            if a boolean array is specified, gamma is fitted using only the corresponding cells
        use_imputed_data: bool, default=True
            use knn smoothed data
        use_size_norm: bool, default=False
            use size normalized data for the fit
        fit_offset: bool, default=True
            Fit with offset
        fixperc_q: bool, default=False
            (when fit_offset==False) Wether to fix the offset to a lower percentile of the unspliced
        weighted: bool, default=True
            use weights for the least squares fit
        weights: string or np.ndarray, default="maxmin_diag"
            the method to determine the weights of the least squares fit.
            "maxmin_diag", "maxmin", "sum", "prod", "maxmin_weighted" are supported
            if a 2d np.ndarray is provided the entry (i,j) is the weight of the cell j when fitting gamma to gene i
        limit_gamma: np.ndarray, default=True
            whether to limit gamma when unspliced is much higher than spliced
        maxmin_perc: List[float], default=[2,98]
            the percentile to use if weights = "maxmin" or "maxmin_diag"

        Returns
        -------
        Nothing it just creates the attributes:
        gammas: np.ndarray
            the vector of the gammas fit to each gene
        q: np.ndarray
            the vector of offsets of the fit
        R2: np.ndarray (optional)
            The vector of squared coefficient of determination

        """
        if steady_state_bool:
            self.steady_state = steady_state_bool
        else:
            self.steady_state = np.ones(self.S.shape[1], dtype=bool)

        if use_imputed_data:
            if use_size_norm:
                tmpS = self.Sx_sz
                tmpU = self.Ux_sz
            else:
                tmpS = self.Sx
                tmpU = self.Ux
        else:
            if use_size_norm:
                tmpS = self.S_sz
                tmpU = self.U_sz
            else:
                tmpS = self.S
                tmpU = self.U

        if weighted:
            if type(weights) is np.ndarray:
                W = weights
            elif weights == "sum":
                W = (tmpS / np.percentile(tmpS, 99, 1)[:, None]) + (tmpU / np.percentile(tmpU, 99, 1)[:, None])
            elif weights == "prod":
                W = (tmpS / np.percentile(tmpS, 99, 1)[:, None]) * (tmpU / np.percentile(tmpU, 99, 1)[:, None])
            elif weights == "maxmin_weighted":
                # Slightly smoother than just takin top and bottom percentile
                down, up = np.percentile(tmpS, maxmin_perc, 1)  # Do this asymmetrically, data is sparse!
                Srange = np.clip(tmpS, down[:, None], up[:, None])
                Srange -= Srange.min(1)[:, None]
                Srange /= Srange.max(1)[:, None]
                W = 0.5 * (Srange**maxmin_weighted_pow + (1 - Srange)**maxmin_weighted_pow)
            elif weights == "maxmin":
                down, up = np.percentile(tmpS, maxmin_perc, 1)  # Do this asymmetrically, data is sparse!
                W = ((tmpS <= down[:, None]) | (tmpS >= up[:, None])).astype(float)
            elif weights == "maxmin_diag":
                denom_Sx = np.percentile(self.Sx, 99.9, 1)
                if np.sum(denom_Sx == 0):
                    denom_Sx[denom_Sx == 0] = np.maximum(np.max(self.Sx[denom_Sx == 0, :], 1), 0.001)
                denom_Ux = np.percentile(self.Ux, 99.9, 1)
                if np.sum(denom_Ux == 0):
                    denom_Ux[denom_Ux == 0] = np.maximum(np.max(self.Ux[denom_Ux == 0, :], 1), 0.001)
                Sx_maxnorm = self.Sx / denom_Sx[:, None]
                Ux_maxnorm = self.Ux / denom_Ux[:, None]
                X = Sx_maxnorm + Ux_maxnorm
                down, up = np.percentile(X, maxmin_perc, axis=1)
                W = ((X <= down[:, None]) | (X >= up[:, None])).astype(float)
            elif weights == "maxmin_double":
                denom_Sx = np.percentile(self.Sx, 99.9, 1)
                denom_Sx[denom_Sx == 0] = np.maximum(np.max(self.Sx[denom_Sx == 0, :], 1), 0.001)
                denom_Ux = np.percentile(self.Ux, 99.9, 1)
                denom_Ux[denom_Ux == 0] = np.maximum(np.max(self.Ux[denom_Ux == 0, :], 1), 0.001)
                Sx_maxnorm = self.Sx / denom_Sx[:, None]
                Ux_maxnorm = self.Ux / denom_Ux[:, None]
                X = Sx_maxnorm + Ux_maxnorm
                down, up = np.percentile(X, maxmin_perc, axis=1)
                W = ((X <= down[:, None]) | (X >= up[:, None])).astype(float)
                down, up = np.percentile(self.Sx, maxmin_perc, 1)
                W += ((self.Sx <= down[:, None]) | (self.Sx >= up[:, None])).astype(float)
                
        if fit_offset:
            if weighted:
                self.gammas, self.q, self.R2 = fit_slope_weighted_offset(tmpU[:, self.steady_state],
                                                                         tmpS[:, self.steady_state],
                                                                         W,
                                                                         return_R2=True,
                                                                         limit_gamma=limit_gamma)
            else:
                if limit_gamma:
                    logging.warning("limit_gamma not implemented with this settings")
                self.gammas, self.q = fit_slope_offset(tmpU[:, self.steady_state],
                                                       tmpS[:, self.steady_state])
        elif fixperc_q:
            if weighted:
                self.gammas, self.q = fit_slope_weighted_offset(tmpU[:, self.steady_state],
                                                                tmpS[:, self.steady_state],
                                                                W, fixperc_q=True, limit_gamma=limit_gamma)
            else:
                if limit_gamma:
                    logging.warning("limit_gamma not implemented with this settings")
                self.gammas, self.q = fit_slope_offset(tmpU[:, self.steady_state],
                                                       tmpS[:, self.steady_state],
                                                       fixperc_q=True)
        else:
            if weighted:
                self.gammas, self.R2 = fit_slope_weighted(tmpU[:, self.steady_state],
                                                          tmpS[:, self.steady_state],
                                                          W, 
                                                          return_R2=True,
                                                          limit_gamma=limit_gamma)
                self.q = np.zeros_like(self.gammas)
            else:
                if limit_gamma:
                    logging.warning("limit_gamma not implemented with this settings")
                self.gammas = fit_slope(tmpU[:, self.steady_state],
                                        tmpS[:, self.steady_state])
                self.q = np.zeros_like(self.gammas)

        # Fix gammas
        self.gammas[~np.isfinite(self.gammas)] = 0

    def filter_genes_good_fit(self, minR: float=0.1, min_gamma: float=0.01) -> None:
        """For backwards compatibility a wrapper around filter_genes_by_phase_portrait
        """
        return self.filter_genes_by_phase_portrait(minR2=minR, min_gamma=min_gamma, minCorr=None)

    def filter_genes_by_phase_portrait(self, minR2: float=0.1, min_gamma: float=0.01, minCorr: float=0.1) -> None:
        """Use the coefficient of determination to filter away genes that have an irregular/complex phase portrait

        Arguments
        ---------
        minR2: float, default=0.1
            Filter away low coefficient of determination fits. If None this filtering will be skipped
        min_gamma: float, default=0.01
            Filter away low gammas. If None this filtering will be skipped
        minCorr: flaot, default=0.2
            Filter away low spliced-usnpliced correlation. If None this filtering will be skipped

        Returns
        -------
        Nothing but modifies it filters out the genes that do not satisfy the conditions
        This affects: "U", "U_sz", "U_norm", "Ux", "Ux_sz", "Ux_norm", "S", "S_sz", "S_norm", "Sx", "Sx_sz", "Sx_norm", "gammas", "q", "R2"
        """

        def paired_correlation_rows(A: np.array, B: np.array) -> np.array:
            A_m = A - A.mean(1)[:, None]
            B_m = B - B.mean(1)[:, None]
            return (A_m * B_m).sum(1) / (np.linalg.norm(A_m, 2, 1) * np.linalg.norm(B_m, 2, 1))

        # @numba.njit()
        # def paired_correlation_rows(A, B):
        #     res = np.zeros(A.shape[0])
        #     for i in range(A.shape[0]):
        #         a = A[i,:] - np.sum(A[i,:]) / A.shape[1]
        #         b = B[i,:] - np.sum(B[i,:]) / B.shape[1]
        #         res[i] = np.sum(a * b) / (np.sqrt(np.sum(a*a)) * np.sqrt(np.sum(b*b)))
        #     return res 

        tmp_filter = np.ones(self.gammas.shape, dtype=bool)
        if minR2 is not None:
            # NOTE Should be: tmp_filter = np.sqrt(self.R2) > minR but since the fit is weighted and constrained R2 can be negative
            R2_corrected = np.sqrt(np.abs(self.R2)) * np.sign(self.R2)
            tmp_filter = tmp_filter & (R2_corrected > minR2)
        if min_gamma is not None:
            tmp_filter = tmp_filter & (self.gammas > min_gamma)
        if minCorr is not None:
            Corr = paired_correlation_rows(self.Sx_sz, self.Ux_sz)
            tmp_filter = tmp_filter & (Corr > minCorr)
        # Perform the filtering
        self.ra = {k: v[tmp_filter] for k, v in self.ra.items()}
        matrixes2filter = ["U", "U_sz", "U_norm", "Ux", "Ux_sz", "Ux_norm",
                           "S", "S_sz", "S_norm", "Sx", "Sx_sz", "Sx_norm"]
        vectors2filter = ["gammas", "q", "R2"]
        for name_attr in matrixes2filter:
            if hasattr(self, name_attr):
                setattr(self, name_attr, getattr(self, name_attr)[tmp_filter, :])
        for name_attr in vectors2filter:
            if hasattr(self, name_attr):
                setattr(self, name_attr, getattr(self, name_attr)[tmp_filter])

    def predict_U(self, which_gamma: str="gammas", which_S: str= "Sx_sz", which_offset: str="q") -> None:
        """Predict U (gamma * S) given the gamma model fit

        Arguments
        ---------
        which_gamma: str, default="gammas"
            name of the attribute to use as gamma
        which_S: str, default="Sx_sz"
            name of the attribute to use as S
        which_offset: str, default="q"
            name of the attribute containing the offset

        Returns
        ------
        Noting but it creates the attribute
        Upred: np.ndarray
           unspliced estimated as `gamma * S`
        """
        self.which_S_for_pred = which_S
        if which_offset is None:
            if hasattr(self, "q_W") or hasattr(self, "q"):
                logging.warn("Predicting U without intercept but intercept was previously fit! Set which_offset='q' or 'q_W' ")
            self.Upred = getattr(self, which_gamma)[:, None] * getattr(self, which_S)
            # self.Upred = self.gammas[:, None] * self.Sx_sz
        else:
            self.Upred = getattr(self, which_gamma)[:, None] * getattr(self, which_S) + getattr(self, which_offset)[:, None]

    def calculate_velocity(self, kind: str="residual", eps: float=None) -> None:
        """Calculate velocity

        Arguments
        ---------
        kind: str, default="residual"
            "residual" calculates the velocity as U_measured - U_predicted
        
        eps: float, default=None
            if specified, velocities with absolute value smaller than eps * range_of_U will be set to 0
            if None this step will be skipped
        
        Results
        -------
        Nothing but it creates the attribute:
        velocity: np.ndarray
            U_measured - U_predicted

        """
        if kind == "residual":
            if self.which_S_for_pred == "Sx_sz":
                self.velocity = self.Ux_sz - self.Upred
            elif self.which_S_for_pred == "Sx":
                self.velocity = self.Ux - self.Upred
            else:
                NotImplementedError(f"Not implemented with which_S = {self.which_S_for_pred}")
        else:
            raise NotImplementedError(f"Velocity calculation kind={kind} is not implemented")

        if eps:
            minimal_signed_res = self.Upred.max(1) * eps
            self.velocity[np.abs(self.velocity) < minimal_signed_res[:, None]] = 0

    def calculate_shift(self, assumption: str="constant_velocity", delta_t: float=1) -> None:
        """Find the change (deltaS) in gene expression for every cell

        Arguments
        ---------
        assumption: str, default="constant_velocity"
            constant_velocity (described in the paper as Model I)
            constant_unspliced (described in the paper as Model II)
        delta_t: float, default=1
            the time step for extrapolation

        Returns
        -------
        Nothing it only creates the following attributes
        delta_S: np.ndarray
            The variation in gene expression
        """
        if assumption == "constant_velocity":
            self.delta_S = delta_t * self.velocity
        elif assumption == "constant_unspliced":
            # Ux_sz = self.Ux_sz - offset; Ux_sz[Ux_sz<0] = 0
            # maybe I should say ratio see below
            Ux_szo = self.Ux_sz - self.q[:, None]
            Ux_szo[Ux_szo < 0] = 0
            egt = np.exp(-self.gammas * delta_t)[:, None]
            self.delta_S = self.Sx_sz * egt + (1 - egt) * Ux_szo / self.gammas[:, None] - self.Sx_sz
        else:
            raise NotImplementedError(f"Assumption {assumption} is not implemented")

    def extrapolate_cell_at_t(self, delta_t: float=1, clip: bool=True) -> None:
        """Extrapolate the gene expression profile for each cell after delta_t
        
        Arguments
        ---------
        delta_t: float, default=1
            the time step considered for the extrapolation
        clip: bool, default=True
            If True negative values are clipped to zero

        Returns
        -------
        Nothing but it creates the attributes:
        Sx_sz_t: np.ndarray
            the extrapolated expression profile
        used_delta_t: float
            stores delta_t for future usage
        """
        if self.which_S_for_pred == "Sx_sz":
            self.Sx_sz_t = self.Sx_sz + delta_t * self.delta_S
            if clip:
                self.Sx_sz_t = np.clip(self.Sx_sz_t, 0, None)
                self.used_delta_t = delta_t
        elif self.which_S_for_pred == "Sx":
            self.Sx_t = self.Sx + delta_t * self.delta_S
            if clip:
                self.Sx_t = np.clip(self.Sx_t, 0, None)
                self.used_delta_t = delta_t
        else:
            NotImplementedError("not implemented for other situations other than Sx or Sx_sz")

    def perform_TSNE(self, n_dims: int=2, perplexity: float=30, initial_pos: np.ndarray=None,
                     theta: float=0.5, n_pca_dim: int=None, max_iter: int=1000) -> None:
        """Perform TSNE on the PCA using barnes hut approximation
        """
        # Perform TSNE
        logging.debug("Running bhtsne")
        if initial_pos is None:
            initial_pos = "random"
        bh_tsne = TSNE(n_components=n_dims, perplexity=perplexity, angle=theta, init=initial_pos, n_iter=max_iter)
        self.ts = bh_tsne.fit_transform(self.pcs[:, :n_pca_dim])

    def estimate_transition_prob(self, hidim: str="Sx_sz", embed: str="ts", transform: str="sqrt",
                                 ndims: int=None, n_sight: int=None, psc: float=None,
                                 knn_random: bool=True, sampled_fraction: float=0.3,
                                 sampling_probs: Tuple[float, float]=(0.5, 0.1), max_dist_embed: float=None,
                                 n_jobs: int=4, threads: int=None, calculate_randomized: bool=True,
                                 random_seed: int=15071990, **kwargs) -> None:
        """Use correlation to estimate transition probabilities for every cells to its embedding neighborhood
        
        Arguments
        ---------
        hidim: str, default="Sx_sz"
            The name of the attribute containing the high dimensional space. It will be retrieved as getattr(self, hidim)
            The updated vector at time t is assumed to be getattr(self, hidim + "_t")
            Appending .T to the string will transpose the matrix (useful in case we want to use S or Sx)
        embed: str, default="ts"
            The name of the attribute containing the embedding. It will be retrieved as getattr(self, embed)
        transform: str, default="sqrt"
            The transformation that is applies on the high dimensional space.
            If None the raw data will be used
        ndims: int, default=None
            The number of dimensions of the high dimensional space to work with. If None all will be considered
            It makes sense only when using principal components
        n_sight: int, default=None (also n_neighbors)
            The number of neighbors to take into account when performing the projection
        psc: float, default=None
            pseudocount added in variance normalizing transform
            If None, 1 would be used for log, 0 otherwise
        knn_random: bool, default=True
            whether to random sample the neighborhoods to speedup calculation
        sampling_probs: Tuple, default=(0.5, 1)
        max_dist_embed: float, default=None
            CURRENTLY NOT USED
            The maximum distance allowed
            If None it will be set to 0.25 * average_distance_two_points_taken_at_random
        n_jobs: int, default=4
            number of jobs to calculate knn
            this only applies to the knn search, for the more time consuming correlation computation see threads
        threads: int, default=None
            The threads will be used for the actual correlation computation by default half of the total.
        calculate_randomized: bool, default=True
            Calculate the transition probabilities with randomized residuals.
            This can be plotted downstream as a negative control and can be used to adjust the visualization scale of the velocity field.
        random_seed: int, default=15071990
            Random seed to make knn_random mode reproducible
        
        Returns
        -------
        """

        numba_random_seed(random_seed)
        self.which_hidim = hidim

        if "n_neighbors" in kwargs:
            n_neighbors = kwargs.pop("n_neighbors")
            if len(kwargs) > 0:
                logging.warning(f"keyword arguments were passed but could not be interpreted {kwargs}")
        else:
            n_neighbors = None

        if n_sight is None and n_neighbors is None:
            n_neighbors = int(self.S.shape[1] / 5)

        if (n_sight is not None) and (n_neighbors is not None) and n_neighbors != n_sight:
            raise ValueError("n_sight and n_neighbors are different names for the same parameter, they cannot be set differently")

        if n_sight is not None and n_neighbors is None:
            n_neighbors = n_sight

        if psc is None:
            if transform == "log" or transform == "logratio":
                psc = 1.
            elif transform == "sqrt":
                psc = 1e-10  # for numerical stablity
            else:  # transform == "linear":
                psc = 0

        if knn_random:
            np.random.seed(random_seed)
            self.corr_calc = "knn_random"
            if "pcs" in hidim:  # sic
                hi_dim = np.array(getattr(self, hidim).T[:, :ndims], order="C")
                hi_dim_t = np.array(getattr(self, hidim + "_t").T[:, :ndims], order="C")
            else:
                if ndims is not None:
                    raise ValueError(f"ndims was set to {ndims} but hidim != 'pcs'. Set ndims = None for hidim='{hidim}'")
                hi_dim = getattr(self, hidim)  # [:, :ndims]
                hi_dim_t = hi_dim + self.used_delta_t * self.delta_S  # [:, :ndims] [:, :ndims]
                if calculate_randomized:
                    self.delta_S_rndm = np.copy(self.delta_S)
                    permute_rows_nsign(self.delta_S_rndm)
                    hi_dim_t_rndm = hi_dim + self.used_delta_t * self.delta_S_rndm
                
            embedding = getattr(self, embed)
            self.embedding = embedding
            logging.debug("Calculate KNN in the embedding space")
            nn = NearestNeighbors(n_neighbors=n_neighbors + 1, n_jobs=n_jobs)
            nn.fit(embedding)  # NOTE should support knn in high dimensions
            self.embedding_knn = nn.kneighbors_graph(mode="connectivity")

            # Pick random neighbours and prune the rest
            neigh_ixs = self.embedding_knn.indices.reshape((-1, n_neighbors + 1))
            p = np.linspace(sampling_probs[0], sampling_probs[1], neigh_ixs.shape[1])
            p = p / p.sum()

            # There was a problem of API consistency because the random.choice can pick the diagonal value (or not)
            # resulting self.corrcoeff with different number of nonzero entry per row.
            # Not updated yet not to break previous analyses
            # Fix is substituting below `neigh_ixs.shape[1]` with `np.arange(1,neigh_ixs.shape[1]-1)`
            # I change it here since I am doing some breaking changes
            sampling_ixs = np.stack((np.random.choice(neigh_ixs.shape[1],
                                                      size=(int(sampled_fraction * (n_neighbors + 1)),),
                                                      replace=False,
                                                      p=p) for i in range(neigh_ixs.shape[0])), 0)
            self.sampling_ixs = sampling_ixs
            neigh_ixs = neigh_ixs[np.arange(neigh_ixs.shape[0])[:, None], sampling_ixs]
            nonzero = neigh_ixs.shape[0] * neigh_ixs.shape[1]
            self.embedding_knn = sparse.csr_matrix((np.ones(nonzero),
                                                    neigh_ixs.ravel(),
                                                    np.arange(0, nonzero + 1, neigh_ixs.shape[1])),
                                                   shape=(neigh_ixs.shape[0],
                                                          neigh_ixs.shape[0]))

            logging.debug(f"Correlation Calculation '{self.corr_calc}'")
            if transform == "log":
                delta_hi_dim = hi_dim_t - hi_dim
                self.corrcoef = colDeltaCorLog10partial(hi_dim, np.log10(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim), neigh_ixs, threads=threads, psc=psc)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                    self.corrcoef_random = colDeltaCorLog10partial(hi_dim, np.log10(np.abs(delta_hi_dim_rndm) + psc) * np.sign(delta_hi_dim_rndm), neigh_ixs, threads=threads, psc=psc)
            elif transform == "logratio":
                log2hidim = np.log2(hi_dim + psc)
                delta_hi_dim = np.log2(np.abs(hi_dim_t) + psc) - log2hidim
                self.corrcoef = colDeltaCorpartial(log2hidim, delta_hi_dim, neigh_ixs, threads=threads)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    delta_hi_dim_rndm = np.log2(np.abs(hi_dim_t_rndm) + psc) - log2hidim
                    self.corrcoef_random = colDeltaCorpartial(log2hidim, delta_hi_dim_rndm, neigh_ixs, threads=threads)
            elif transform == "linear":
                self.corrcoef = colDeltaCorpartial(hi_dim, hi_dim_t - hi_dim, neigh_ixs, threads=threads)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    self.corrcoef_random = colDeltaCorpartial(hi_dim, hi_dim_t_rndm - hi_dim, neigh_ixs, threads=threads)
            elif transform == "sqrt":
                delta_hi_dim = hi_dim_t - hi_dim
                self.corrcoef = colDeltaCorSqrtpartial(hi_dim, np.sqrt(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim), neigh_ixs, threads=threads, psc=psc)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                    self.corrcoef_random = colDeltaCorSqrtpartial(hi_dim, np.sqrt(np.abs(delta_hi_dim_rndm) + psc) * np.sign(delta_hi_dim_rndm), neigh_ixs, threads=threads, psc=psc)
            else:
                raise NotImplementedError(f"transform={transform} is not a valid parameter")
            np.fill_diagonal(self.corrcoef, 0)
            if np.any(np.isnan(self.corrcoef)):
                self.corrcoef[np.isnan(self.corrcoef)] = 1
                logging.warning("Nans encountered in corrcoef and corrected to 1s. If not identical cells were present it is probably a small isolated cluster converging after imputation.")
            if calculate_randomized:
                np.fill_diagonal(self.corrcoef_random, 0)
                if np.any(np.isnan(self.corrcoef_random)):
                    self.corrcoef_random[np.isnan(self.corrcoef_random)] = 1
                    logging.warning("Nans encountered in corrcoef_random and corrected to 1s. If not identical cells were present it is probably a small isolated cluster converging after imputation.")
            logging.debug(f"Done Correlation Calculation")
        else:
            self.corr_calc = "full"
            if "pcs" in hidim:  # sic
                hi_dim = np.array(getattr(self, hidim).T[:, :ndims], order="C")
                hi_dim_t = np.array(getattr(self, hidim + "_t").T[:, :ndims], order="C")
            else:
                if ndims is not None:
                    raise ValueError(f"ndims was set to {ndims} but hidim != 'pcs'. Set ndims = None for hidim='{hidim}'")
                hi_dim = getattr(self, hidim)  # [:, :ndims]
                hi_dim_t = hi_dim + self.used_delta_t * self.delta_S  # [:, :ndims] [:, :ndims]
                if calculate_randomized:
                    self.delta_S_rndm = np.copy(self.delta_S)
                    permute_rows_nsign(self.delta_S_rndm)
                    hi_dim_t_rndm = hi_dim + self.used_delta_t * self.delta_S_rndm
                
            embedding = getattr(self, embed)
            self.embedding = embedding
            logging.debug("Calculate KNN in the embedding space")
            nn = NearestNeighbors(n_neighbors=n_neighbors + 1, n_jobs=n_jobs)
            nn.fit(embedding)
            self.embedding_knn = nn.kneighbors_graph(mode="connectivity")
            
            logging.debug("Correlation Calculation 'full'")
            if transform == "log":
                delta_hi_dim = hi_dim_t - hi_dim
                self.corrcoef = colDeltaCorLog10(hi_dim, np.log10(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim), threads=threads, psc=psc)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                    self.corrcoef_random = colDeltaCorLog10(hi_dim, np.log10(np.abs(delta_hi_dim_rndm) + psc) * np.sign(delta_hi_dim_rndm), threads=threads, psc=psc)
            elif transform == "logratio":
                log2hidim = np.log2(hi_dim + psc)
                delta_hi_dim = np.log2(np.abs(hi_dim_t) + psc) - log2hidim
                self.corrcoef = colDeltaCor(log2hidim, delta_hi_dim, threads=threads)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    delta_hi_dim_rndm = np.log2(np.abs(hi_dim_t_rndm) + 1) - log2hidim
                    self.corrcoef_random = colDeltaCor(log2hidim, delta_hi_dim_rndm, threads=threads)
            elif transform == "linear":
                self.corrcoef = colDeltaCor(hi_dim, hi_dim_t - hi_dim, threads=threads)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    self.corrcoef_random = colDeltaCor(hi_dim, hi_dim_t_rndm - hi_dim, threads=threads, psc=psc)
            elif transform == "sqrt":
                delta_hi_dim = hi_dim_t - hi_dim
                self.corrcoef = colDeltaCorSqrt(hi_dim, np.sqrt(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim), threads=threads, psc=psc)
                if calculate_randomized:
                    logging.debug(f"Correlation Calculation for negative control")
                    delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                    self.corrcoef_random = colDeltaCorSqrt(hi_dim, np.sqrt(np.abs(delta_hi_dim_rndm) + psc) * np.sign(delta_hi_dim_rndm), threads=threads, psc=psc)
            else:
                raise NotImplementedError(f"transform={transform} is not a valid parameter")
            np.fill_diagonal(self.corrcoef, 0)
            if calculate_randomized:
                np.fill_diagonal(self.corrcoef_random, 0)

    def calculate_embedding_shift(self, sigma_corr: float=0.05, expression_scaling: bool=True, scaling_penalty: float=1.) -> None:
        """Use the transition probability to project the velocity direction on the embedding

        Arguments
        ---------
        sigma_corr: float, default=0.05
            the kernel scaling
        expression_scaling: bool, default=True
            rescale arrow intensity penalizing arrows that explain very small amount of expression differences
        scaling_penalty: float, default=1
            Higher values correspond to a stronger penalty
        

        Returns
        -------
        Nothing but it creates the following attributes:
        transition_prob: np.ndarray
            the transition probability calculated using the exponential kernel on the correlation coefficient
        delta_embedding: np.ndarray
            The resulting vector
        """
        # Kernel evaluation
        logging.debug("Calculate transition probability")

        if self.corr_calc == "full" or self.corr_calc == "knn_random":
            # NOTE maybe sparse matrix here are slower than dense
            # NOTE if knn_random this could be made much faster either using sparse matrix or neigh_ixs
            self.transition_prob = np.exp(self.corrcoef / sigma_corr) * self.embedding_knn.A  # naive
            self.transition_prob /= self.transition_prob.sum(1)[:, None]
            if hasattr(self, "corrcoef_random"):
                logging.debug("Calculate transition probability for negative control")
                self.transition_prob_random = np.exp(self.corrcoef_random / sigma_corr) * self.embedding_knn.A  # naive
                self.transition_prob_random /= self.transition_prob_random.sum(1)[:, None]

            unitary_vectors = self.embedding.T[:, None, :] - self.embedding.T[:, :, None]  # shape (2,ncells,ncells)
            with np.errstate(divide='ignore', invalid='ignore'):
                unitary_vectors /= np.linalg.norm(unitary_vectors, ord=2, axis=0)  # divide by L2
                np.fill_diagonal(unitary_vectors[0, ...], 0)  # fix nans
                np.fill_diagonal(unitary_vectors[1, ...], 0)

            self.delta_embedding = (self.transition_prob * unitary_vectors).sum(2)
            self.delta_embedding -= (self.embedding_knn.A * unitary_vectors).sum(2) / self.embedding_knn.sum(1).A.T
            self.delta_embedding = self.delta_embedding.T

            if expression_scaling:
                hi_dim = getattr(self, self.which_hidim)
                estim_delta = hi_dim.dot(self.transition_prob.T) - hi_dim.dot((self.embedding_knn.A / self.embedding_knn.sum(1).A).T)
                cos_proj = (self.delta_S * estim_delta).sum(0) / np.sqrt((estim_delta**2).sum(0))
                self.scaling = np.clip(cos_proj / scaling_penalty, 0, 1)
                self.delta_embedding = self.delta_embedding * self.scaling[:, None]

            if hasattr(self, "corrcoef_random"):
                self.delta_embedding_random = (self.transition_prob_random * unitary_vectors).sum(2)
                self.delta_embedding_random -= (self.embedding_knn.A * unitary_vectors).sum(2) / self.embedding_knn.sum(1).A.T
                self.delta_embedding_random = self.delta_embedding_random.T

                if expression_scaling:
                    estim_delta_rndm = hi_dim.dot(self.transition_prob_random.T) - hi_dim.dot((self.embedding_knn.A / self.embedding_knn.sum(1).A).T)
                    cos_proj_rndm = (self.delta_S_rndm * estim_delta_rndm).sum(0) / np.sqrt((estim_delta_rndm**2).sum(0))
                    self.scaling_rndm = np.clip(cos_proj_rndm / scaling_penalty, 0, 1)
                    self.delta_embedding_random = self.delta_embedding_random * self.scaling_rndm[:, None]
        else:
            # NOTE should implement a version with cython
            raise NotImplementedError(f"Weird value self.corr_calc={self.corr_calc}")

    def calculate_grid_arrows(self, embed: str="embedding", smooth: float=0.5, steps: Tuple=(40, 40),
                              n_neighbors: int=100, n_jobs: int=4) -> None:
        """Calculate the velocity using a points on a regular grid and a gaussian kernel

        Note: the function should work also for n-dimensional grid

        Arguments
        ---------
        embed: str, default=embedding
            The name of the attribute containing the embedding. It will be retrieved as getattr(self, embed)
            The difference vector is getattr(self, 'delta' + '_' + embed)
        smooth: float, smooth=0.5
            Higher value correspond to taking in consideration further points
            the standard deviation of the gaussian kernel is smooth * stepsize
        steps: tuple, default
            the number of steps in the grid for each axis
        n_neighbors:
            number of neighbors to use in the calculation, bigger number should not change too much the results..
            ...as soon as smooth is small
            Higher value correspond to slower execution time
        n_jobs:
            number of processes for parallel computing

        Returns
        -------
        Nothing but it sets the attributes:
        flow_embedding: np.ndarray
            the coordinates of the embedding
        flow_grid: np.ndarray
            the gridpoints
        flow: np.ndarray
            vector field coordinates
        flow_magnitude: np.ndarray
            magnitude of each vector on the grid
        total_p_mass: np.ndarray
            density at each point of the grid

        """
        embedding = getattr(self, embed)
        if hasattr(self, f"delta_{embed}"):
            delta_embedding = getattr(self, f"delta_{embed}")
            if hasattr(self, "corrcoef_random"):
                delta_embedding_random = getattr(self, f"delta_{embed}_random")
        else:
            raise KeyError("This embedding does not have a delta_*")
        # Prepare the grid
        grs = []
        for dim_i in range(embedding.shape[1]):
            m, M = np.min(embedding[:, dim_i]), np.max(embedding[:, dim_i])
            m = m - 0.025 * np.abs(M - m)
            M = M + 0.025 * np.abs(M - m)
            gr = np.linspace(m, M, steps[dim_i])
            grs.append(gr)
            
        meshes_tuple = np.meshgrid(*grs)
        gridpoints_coordinates = np.vstack([i.flat for i in meshes_tuple]).T

        nn = NearestNeighbors(n_neighbors=n_neighbors, n_jobs=n_jobs)
        nn.fit(embedding)
        dists, neighs = nn.kneighbors(gridpoints_coordinates)
        
        std = np.mean([(g[1] - g[0]) for g in grs])
        # isotropic gaussian kernel
        gaussian_w = normal.pdf(loc=0, scale=smooth * std, x=dists)
        self.total_p_mass = gaussian_w.sum(1)

        UZ = (delta_embedding[neighs] * gaussian_w[:, :, None]).sum(1) / np.maximum(1, self.total_p_mass)[:, None]  # weighed average
        magnitude = np.linalg.norm(UZ, axis=1)
        # Assign attributes
        self.flow_embedding = embedding
        self.flow_grid = gridpoints_coordinates
        self.flow = UZ
        self.flow_norm = UZ / np.percentile(magnitude, 99.5)
        self.flow_norm_magnitude = np.linalg.norm(self.flow_norm, axis=1)

        if hasattr(self, "corrcoef_random"):
            UZ_rndm = (delta_embedding_random[neighs] * gaussian_w[:, :, None]).sum(1) / np.maximum(1, self.total_p_mass)[:, None]  # weighed average
            magnitude_rndm = np.linalg.norm(UZ, axis=1)
            # Assign attributes
            self.flow_rndm = UZ_rndm
            self.flow_norm_rndm = UZ_rndm / np.percentile(magnitude_rndm, 99.5)
            self.flow_norm_magnitude_rndm = np.linalg.norm(self.flow_norm_rndm, axis=1)

    def prepare_markov(self, sigma_D: np.ndarray, sigma_W: np.ndarray, direction: str="forward", cells_ixs: np.ndarray=None) -> None:
        """Prepare a transition probability for Markov process

        Arguments
        ---------
        sigma_D: float
            the standard deviation used on the locality-limiting component
        sigma_W: float
            the standard deviation used on the noise component
        direction: str, default="backwards"
            whether to diffuse forward of backwards
        cells_ixs: np.ndarray, default=None
            Cells to use, if None all the cells will be considered.

        Returns
        -------
        Nothing but it creates the following attributes:
        tr: np.ndarray
            the transition probability matrix
        
        """
        if cells_ixs is None:
            cells_ixs = np.arange(self.transition_prob.shape[0])

        # NOTE: This implementation is not speed optimized to improve the speed of the implementation:
        # - the C/Fortran contiguity of the transition matrix should be taken into account
        # - a knn implementation would reduce computation
        # - should avoid transformation to and from dense-sparse formats
        if direction == "forward":
            self.tr = np.array(self.transition_prob[cells_ixs, :][:, cells_ixs])
        elif direction == "backwards":
            self.tr = np.array((self.transition_prob[cells_ixs, :][:, cells_ixs]).T, order="C")
        else:
            raise NotImplementedError(f"{direction} is not an implemented direction")
        dist_matrix = squareform(pdist(self.embedding[cells_ixs, :]))
        K_D = gaussian_kernel(dist_matrix, sigma=sigma_D)
        self.tr = self.tr * K_D
        # Fill diagonal with max or the row and sum=1 normalize
        np.fill_diagonal(self.tr, self.tr.max(1))
        self.tr = self.tr / self.tr.sum(1)[:, None]
       
        K_W = gaussian_kernel(dist_matrix, sigma=sigma_W)
        K_W = K_W / K_W.sum(1)[:, None]
        self.tr = 0.8 * self.tr + 0.2 * K_W
        self.tr = self.tr / self.tr.sum(1)[:, None]
        self.tr = scipy.sparse.csr_matrix(self.tr)

    def run_markov(self, starting_p: np.ndarray=None, n_steps: int=2500, mode: str="time_evolution") -> None:
        """Run a Markov process

        Arguments
        ---------
        starting_p: np.ndarray, default=None
            specifies the starting density
            if None is passed an array of 1/self.tr.shape[0] will be created
        n_steps: np.ndarray, default=2500
            Numbers of steps to be performed
        mode: str, default="time_evolution"
            this argument is passed to the Diffusion.diffuse call

        Returns
        -------
        Nothing but it creates the attribute:
        diffused: np.ndarray
            The probability to be found at any of the states
        """
        if starting_p is None:
            starting_p = np.ones(self.tr.shape[0]) / self.tr.shape[0]
        diffusor = Diffusion()
        self.diffused = diffusor.diffuse(starting_p, self.tr, n_steps=n_steps, mode=mode)[0]

    def default_filter_and_norm(self, min_expr_counts: int=None, min_cells_express: int=None,
                                N: int=None, min_avg_U: float=None, min_avg_S: float=None) -> None:
        """Useful function to get started with velocyto: it performs initial filtering and feature selection, it uses some heuristics to determine the thresholds, results might be suboptimal.

        See `the analysis quick start guide <http://velocyto.org/velocyto.py/tutorial/analysis.html>`_ for further info.

        Arguments
        ---------
        min_expr_counts: int, default=None
            filtering condition: the minimum spliced counts
        min_cells_express: int, default=None
            filtering condition: the minimum number of cells expressing the gene
        N: int, default=None
            number of genes selected by the feature selection procedure
        min_avg_U: float, default=None
            if cluster have been specified beforehand (using the function set_clusters) then this is the minimum average unspliced molecules per cluster
        min_avg_S: float, default=None
            if cluster have been specified beforehand (using the function set_clusters) then this is the minimum average spliced molecules per cluster
        """
        logging.warning("DEPRECATION WARNING - the current function is deprecated. Please refer to documetation for default parameters usage")
        if min_expr_counts is None:
            min_expr_counts = max(20, min(100, self.S.shape[1] * 2.25e-3))
        if min_cells_express is None:
            min_cells_express = max(10, min(50, self.S.shape[1] * 1.5e-3))
        if N is None:
            N = max(1000, min(int((self.S.shape[1] / 1000)**(1 / 3) / 0.0008), 5000))
        if min_avg_U is None:
            min_avg_U = 0.01
        if min_avg_S is None:
            min_avg_S = 0.08

        # This is called just to compute the initial cell size, normalized value will be recalculated
        self.normalize("S", size=True, log=False)
        self.normalize("U", size=True, log=False)

        self.score_detection_levels(min_expr_counts=min_expr_counts, min_cells_express=min_cells_express)
        self.filter_genes(by_detection_levels=True)

        self.score_cv_vs_mean(N=N, max_expr_avg=40)
        self.filter_genes(by_cv_vs_mean=True)

        self.score_detection_levels(min_expr_counts=0, min_cells_express=0,
                                    min_expr_counts_U=int(min_expr_counts / 2) + 1,
                                    min_cells_express_U=int(min_cells_express / 2) + 1)
        
        if hasattr(self, "cluster_labels"):
            self.score_cluster_expression(min_avg_U=min_avg_U, min_avg_S=min_avg_S)
            self.filter_genes(by_detection_levels=True, by_cluster_expression=True)
        else:
            self.filter_genes(by_detection_levels=True)
        self.normalize_by_total()
        self.adjust_totS_totU(normalize_total=True)

    def default_fit_preparation(self, k: int=None, n_comps: int=None) -> None:
        """Useful function to get started with velocyto: it performs PCA and kNN smoothing, it uses some heuristics to determine the parameters, results might be suboptimal.

        See `the analysis quick start guide <http://velocyto.org/velocyto.py/tutorial/analysis.html>`_ for further info.

        Arguments
        ---------
        k: int, default=None
            k in k-NearestNeighbours smoothing
        n_comps: int, default=None
            numbed of components in pca
        """
        logging.warning("DEPRECATION WARNING - the current function is deprecated. Please refer to documetation for default parameters usage")
        self.perform_PCA()
        # Choose the number of components to use for the kNN graph
        if n_comps is None:
            n_comps = int(np.where(np.diff(np.diff(np.cumsum(self.pca.explained_variance_ratio_)) > 0.002))[0][0])
        if k is None:
            k = int(min(1000, max(10, np.ceil(self.S.shape[1] * 0.02))))
        self.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True,
                            b_sight=int(min(k * 8, self.S.shape[1] - 1)),
                            b_maxl=int(min(k * 4, self.S.shape[1] - 1)))
        self.normalize_median()

    def _plot_phase_portrait(self, gene: str, gs_i: Any=None) -> None:
        """Plot spliced-unspliced scatterplot resembling phase portrait
        """
        if gene is None:
            plt.subplot(111)
        else:
            plt.subplot(gs_i)
        ix = np.where(self.ra["Gene"] == gene)[0][0]
        scatter_viz(self.Sx_sz[ix, :], self.Ux_sz[ix, :], c=self.colorandum, s=5, alpha=0.4)
        plt.title(gene)
        xnew = np.linspace(0, self.Sx_sz[ix, :].max())
        plt.plot(xnew, self.gammas[ix] * xnew + self.q[ix], c="k")

    def plot_phase_portraits(self, genes: List[str]) -> None:
        """Plot spliced-unspliced scatterplots resembling phase portraits

        Arguments
        ---------
        genes: List[str]
            A list of gene symbols.
        """
        n = len(genes)
        sqrtn = int(np.ceil(np.sqrt(n)))
        gs = plt.GridSpec(sqrtn, int(np.ceil(n / sqrtn)))
        for i, gn in enumerate(genes):
            self._plot_phase_portrait(gn, gs[i])
    
    def plot_grid_arrows(self, quiver_scale: Union[str, float]="auto", scale_type: str= "relative", min_mass: float=1, min_magnitude: float=None,
                         scatter_kwargs_dict: Dict= None, plot_dots: bool=False, plot_random: bool=False, **quiver_kwargs: Any) -> None:
        """Plots vector field averaging velocity vectors on a grid

        Arguments
        ---------
        quiver_scale: float, default="auto"
            Rescaling factor applied to the arrow field to enhance visibility
            If "auto" the scale is selected using the randomized (negative) control (even if `plot_random`=False)
            If a float is provided the interpretation of the value depends on the parameter `scale_type`, see below.
            NOTE: In the situation where "auto" is resulting in very small or big velocities, pass a float to this parameter
            The float will be interpreted as a scaling, importantly both the data and the control will be scaled
            in this way you can rescale the velocity arbitrarily without the risk of observing just an overfit of the noise
        scale_type: str, default="relative"
            How to interpret `quiver_scale`:
            If "relative" (default) the value will be used as a scaling factor and multiplied by the value from "auto"
            (it follows that quiver_scale="auto" is equivalent to quiver_scale=1)
            If "absolute" the value will be passed to the matplotlib quiver function (not recommended if you are not sure what this implies)
        min_mass: float, default=1
            the minimum density around a grid point for it to be considered and plotted
        min_magnitude: float, default=None
            the minimum magnitude of the velocity for it to be considered and plotted
        scatter_kwargs_dict: dict, default=None
            a dictionary of keyword arguments to pass to scatter
            by default the following are passed: s=20, zorder=-1, alpha=0.2, lw=0, c=self.colorandum. But they can be overridden.
        plot_dots: bool, default=True
            whether to plot dots in correspondence of all low velocity grid points
        plot_random: bool, default=True
            whether to plot the randomized control next to the plot
        **quiver_kwargs: dict
            keyword arguments to pass to quiver
            By default the following are passed angles='xy', scale_units='xy', minlength=1.5. But they can be overridden.
        """
        # plt.figure(figsize=(10, 10))
        _quiver_kwargs = {"angles": 'xy', "scale_units": 'xy', "minlength": 1.5}
        _quiver_kwargs.update(quiver_kwargs)
        
        scatter_dict = {"s": 20, "zorder": -1, "alpha": 0.2, "lw": 0, "c": self.colorandum}
        if scatter_kwargs_dict is not None:
            scatter_dict.update(scatter_kwargs_dict)

        # Determine quiver scale
        if scale_type == "relative":
            if hasattr(self, "flow_rndm"):
                plot_scale = np.linalg.norm(np.max(self.flow_grid, 0) - np.min(self.flow_grid, 0), 2)  # Diagonal of the plot
                arrows_scale = np.percentile(np.linalg.norm(self.flow_rndm[self.total_p_mass >= min_mass, :], 2, 1), 90)  # Tipical lenght of an arrow
                if quiver_scale == "auto":
                    quiver_scale = arrows_scale / (plot_scale * 0.0025)
                else:
                    quiver_scale = quiver_scale * arrows_scale / (plot_scale * 0.0025)
            else:
                raise ValueError(""""`scale_type` was set to 'relative' but the randomized control was not computed when running estimate_transition_prob
                Please run estimate_transition_prob or set `scale_type` to `absolute`""")
        else:
            logging.warning("The arrow scale was set to be 'absolute' make sure you know how to properly interpret the plots")
        
        mass_filter = self.total_p_mass < min_mass
        if min_magnitude is None:
            XY, UV = np.copy(self.flow_grid), np.copy(self.flow)
            if not plot_dots:
                UV = UV[~mass_filter, :]
                XY = XY[~mass_filter, :]
            else:
                UV[mass_filter, :] = 0
        else:
            XY, UV = np.copy(self.flow_grid), np.copy(self.flow_norm)
            if not plot_dots:
                UV = UV[~(mass_filter | (self.flow_norm_magnitude < min_magnitude)), :]
                XY = XY[~(mass_filter | (self.flow_norm_magnitude < min_magnitude)), :]
            else:
                UV[mass_filter | (self.flow_norm_magnitude < min_magnitude), :] = 0

        if plot_random:
            if min_magnitude is None:
                XY, UV_rndm = np.copy(self.flow_grid), np.copy(self.flow_rndm)
                if not plot_dots:
                    UV_rndm = UV_rndm[~mass_filter, :]
                    XY = XY[~mass_filter, :]
                else:
                    UV_rndm[mass_filter, :] = 0
            else:
                XY, UV_rndm = np.copy(self.flow_grid), np.copy(self.flow_norm_rndm)
                if not plot_dots:
                    UV_rndm = UV_rndm[~(mass_filter | (self.flow_norm_magnitude_rndm < min_magnitude)), :]
                    XY = XY[~(mass_filter | (self.flow_norm_magnitude_rndm < min_magnitude)), :]
                else:
                    UV_rndm[mass_filter | (self.flow_norm_magnitude_rndm < min_magnitude), :] = 0

            plt.subplot(122)
            plt.title("Randomized")
            plt.scatter(self.flow_embedding[:, 0], self.flow_embedding[:, 1], **scatter_dict)
            plt.quiver(XY[:, 0], XY[:, 1], UV_rndm[:, 0], UV_rndm[:, 1],
                       scale=quiver_scale, zorder=20000, **_quiver_kwargs)
            plt.axis("off")
            plt.subplot(121)
            plt.title("Data")

        plt.scatter(self.flow_embedding[:, 0], self.flow_embedding[:, 1], **scatter_dict)
        plt.quiver(XY[:, 0], XY[:, 1], UV[:, 0], UV[:, 1],
                   scale=quiver_scale, zorder=20000, **_quiver_kwargs)
        plt.axis("off")

    def plot_arrows_embedding(self, choice: Union[str, int]="auto", quiver_scale: Union[str, float]="auto", scale_type: str="relative",
                              plot_scatter: bool=False, scatter_kwargs: Dict={}, color_arrow: str="cluster",
                              new_fig: bool=False, plot_random: bool=True, **quiver_kwargs: Any) -> None:
        """Plots velocity on the embedding cell-wise
        
        Arguments
        ---------
        choice: int, default = "auto"
            the number of cells to randomly pick to plot the arrows (To avoid overcrowding)
        quiver_scale: float, default="auto"
            Rescaling factor applied to the arrow field to enhance visibility
            If "auto" the scale is selected using the randomized (negative) control (even if `plot_random`=False)
            If a float is provided the interpretation of the value depends on the parameter `scale_type`, see below.
            NOTE: Despite a similar option than plot_grid_arrows, here there is no strong motivation to calculate the scale relative to the randomized control
            This is because the randomized doesn't have to have smaller velocity cell-wise, there might be for example
            scattered cells that will have strong velocity but they will, correctly just average out when calculating the average velocity field.
        scale_type: str, default="relative"
            How to interpret `quiver_scale`:
            If "relative" (default) the value will be used as a scaling factor and multiplied by the value from "auto"
            (it follows that quiver_scale="auto" is equivalent to quiver_scale=1)
            If "absolute" the value will be passed to the matplotlib quiver function
        plot_scatter: bool, default = False
            whether to plot the points
        scatter_kwargs: Dict
            A dictionary containing all the keywords arguments to pass to matplotlib scatter
            by default the following are passed: c="0.8", alpha=0.4, s=10, edgecolor=(0, 0, 0, 1), lw=0.3. But they can be overridden.
        color_arrow: str, default = "cluster"
            the color of the arrows, if "cluster" the arrows are colored the same as the cluster
        new_fig: bool, default=False
            whether to create a new figure
        plot_random: bool, default=True
            whether to plot the randomized control next to the plot
        **quiver_kwargs: dict
            keyword arguments to pass to quiver
            By default the following are passed angles='xy', scale_units='xy', minlength=1.5. But they can be overridden.

        Returns
        -------
        Nothing, just plots the tsne with arrows
        """
        if choice == "auto":
            choice = int(self.S.shape[1] / 3)
            logging.warning(f"Only {choice} arrows will be shown to avoid overcrowding, you can choose the exact number setting the `choice` argument")
        _quiver_kwargs = {"angles": 'xy', "scale_units": 'xy', "minlength": 1.5}
        _scatter_kwargs = dict(c="0.8", alpha=0.4, s=10, edgecolor=(0, 0, 0, 1), lw=0.3)
        _scatter_kwargs.update(scatter_kwargs)
        if new_fig:
            if plot_random and hasattr(self, "delta_embedding_random"):
                plt.figure(figsize=(22, 12))
            else:
                plt.figure(figsize=(14, 14))
        
        ix_choice = np.random.choice(self.embedding.shape[0], size=choice, replace=False)

        # Determine quiver scale
        if scale_type == "relative":
            if hasattr(self, "delta_embedding_random"):
                plot_scale = np.linalg.norm(np.max(self.flow_grid, 0) - np.min(self.flow_grid, 0), 2)  # Diagonal of the plot
                arrows_scale = np.percentile(np.linalg.norm(self.delta_embedding_random, 2, 1), 80)  # Tipical length of an arrow
                if quiver_scale == "auto":
                    quiver_scale = arrows_scale / (plot_scale * 0.005)
                else:
                    quiver_scale = quiver_scale * arrows_scale / (plot_scale * 0.005)
            else:
                raise ValueError("""`scale_type` was set to 'relative' but the randomized control was not computed when running estimate_transition_prob
                Please run estimate_transition_prob or set `scale_type` to `absolute`""")
        else:
            logging.warning("The arrow scale was set to be 'absolute' make sure you know how to properly interpret the plots")

        if color_arrow == "cluster":
            colorandum = self.colorandum[ix_choice, :]
        else:
            colorandum = color_arrow

        _quiver_kwargs.update({"color": colorandum})
        _quiver_kwargs.update(quiver_kwargs)

        if plot_random and hasattr(self, "delta_embedding_random"):
            plt.subplot(122)
            plt.title("Randomized")
            if plot_scatter:
                plt.scatter(self.embedding[:, 0], self.embedding[:, 1], **_scatter_kwargs)
            plt.quiver(self.embedding[ix_choice, 0], self.embedding[ix_choice, 1],
                       self.delta_embedding_random[ix_choice, 0], self.delta_embedding_random[ix_choice, 1],
                       scale=quiver_scale, **_quiver_kwargs)
            plt.axis("off")
            plt.subplot(121)
            plt.title("Data")

        if plot_scatter:
            plt.scatter(self.embedding[:, 0], self.embedding[:, 1], **_scatter_kwargs)

        plt.quiver(self.embedding[ix_choice, 0], self.embedding[ix_choice, 1],
                   self.delta_embedding[ix_choice, 0], self.delta_embedding[ix_choice, 1],
                   scale=quiver_scale, **_quiver_kwargs)
        plt.axis("off")
    
    def plot_cell_transitions(self, cell_ix: int=0, alpha: float=0.1, alpha_neigh: float=0.2,
                              cmap_name: str="RdBu_r", plot_arrow: bool=True,
                              mark_cell: bool=True, head_width: int=3) -> None:
        """Plot the probability of a cell to transition to any other cell

        This function is untested
        """
        cmap = plt.cm.get_cmap(name=cmap_name)
        colorandum = np.ones((self.embedding.shape[0], 4))
        colorandum *= 0.3
        colorandum[:, -1] = alpha
        
        plt.scatter(self.embedding[:, 0], self.embedding[:, 1],
                    c=colorandum, s=50, edgecolor="")
        if mark_cell:
            plt.scatter(self.embedding[cell_ix, 0], self.embedding[cell_ix, 1],
                        facecolor="none", s=100, edgecolor="k")
        if plot_arrow:
            plt.arrow(self.embedding[cell_ix, 0], self.embedding[cell_ix, 1],
                      self.delta_embedding[cell_ix, 0], self.delta_embedding[cell_ix, 1],
                      head_width=head_width, length_includes_head=True)
    
    def plot_velocity_as_color(self, gene_name: str=None, cmap: Any= plt.cm.RdBu_r,
                               gs: Any=None, which_tsne: str="ts", **kwargs: Dict) -> None:
        """Plot velocity as color on the Tsne embedding

        Arguments
        ---------
        gene_name: str
            The name of the gene, should be present in self.S
        cmap: maplotlib.cm.Colormap, default=maplotlib.cm.RdBu_r
            Colormap to use, divergent ones are better, RdBu_r is default
            Notice that 0 will be always set as the center of the colormap. (e.g. white in RdBu_r)
        gs: Gridspec subplot
            Gridspec subplot to plot on.
        which_tsne: str, default="ts"
            the name of the attributed where the desired embedding is stored
        **kwargs: dict
            other keywords arguments will be passed to the plt.scatter call

        Returns
        -------
        Nothing
        """

        ix = np.where(self.ra["Gene"] == gene_name)[0][0]
        kwarg_plot = {"alpha": 0.5, "s": 8, "edgecolor": "0.8", "lw": 0.15}
        kwarg_plot.update(kwargs)
        if gs is None:
            fig = plt.figure(figsize=(10, 10))
            plt.subplot(111)
        else:
            plt.subplot(gs)
    
        tsne = getattr(self, which_tsne)
        if self.which_S_for_pred == "Sx_sz":
            tmp_colorandum = self.Sx_sz_t[ix, :] - self.Sx_sz[ix, :]
        else:
            tmp_colorandum = self.Sx_t[ix, :] - self.Sx[ix, :]
        if (np.abs(tmp_colorandum) > 0.00005).sum() < 10:  # If S vs U scatterplot it is flat
            print("S vs U scatterplot it is flat")
            return
        limit = np.max(np.abs(np.percentile(tmp_colorandum, [1, 99])))  # upper and lowe limit / saturation
        tmp_colorandum = tmp_colorandum + limit  # that is: tmp_colorandum - (-limit)
        tmp_colorandum = tmp_colorandum / (2 * limit)  # that is: tmp_colorandum / (limit - (-limit))
        tmp_colorandum = np.clip(tmp_colorandum, 0, 1)

        scatter_viz(tsne[:, 0], tsne[:, 1],
                    c=cmap(tmp_colorandum), **kwarg_plot)
        plt.axis("off")
        plt.title(f"{gene_name}")

    def plot_expression_as_color(self, gene_name: str=None, imputed: bool= True, cmap: Any= plt.cm.Greens,
                                 gs: Any=None, which_tsne: str="ts", **kwargs: Dict) -> None:
        """Plot expression as color on the Tsne embedding

        Arguments
        ---------
        gene_name: str
            The name of the gene, should be present in self.S
        imputed: bool, default=True
            whether to plot the smoothed or the raw data
        cmap: maplotlib.cm.Colormap, default=maplotlib.cm.Greens
            Colormap to use.
        gs: Gridspec subplot
            Gridspec subplot to plot on.
        which_tsne: str, default="ts"
            the name of the attributed where the desired embedding is stored
        **kwargs: dict
            other keywords arguments will be passed to the plt.scatter call

        Returns
        -------
        Nothing
        """
        ix = np.where(self.ra["Gene"] == gene_name)[0][0]
        kwarg_plot = {"alpha": 0.5, "s": 8, "edgecolor": "0.8", "lw": 0.15}
        kwarg_plot.update(kwargs)
        if gs is None:
            fig = plt.figure(figsize=(10, 10))
            plt.subplot(111)
        else:
            plt.subplot(gs)
    
        tsne = getattr(self, which_tsne)
        if imputed:
            if self.which_S_for_pred == "Sx_sz":
                tmp_colorandum = self.Sx_sz[ix, :]
            else:
                tmp_colorandum = self.Sx[ix, :]
        else:
            tmp_colorandum = self.S_sz[ix, :]
            
        tmp_colorandum = tmp_colorandum / np.percentile(tmp_colorandum, 99)
        # tmp_colorandum = np.log2(tmp_colorandum+1)
        tmp_colorandum = np.clip(tmp_colorandum, 0, 1)

        scatter_viz(tsne[:, 0], tsne[:, 1],
                    c=cmap(tmp_colorandum), **kwarg_plot)
        plt.axis("off")
        plt.title(f"{gene_name}")

    def reload_raw(self, substitute: bool=False) -> None:
        """Reload raw data as it was before filtering steps

        Arguments
        ---------
        substitute: bool=False
            if True `S, U, A, ca, ra` will be all overwritten
            if False `S, U, A, ca, ra` will be loaded separately as `raw_S, raw_U, raw_A, raw_ca, raw_ra`
        """
        if substitute:
            ds = loompy.connect(self.loom_filepath)
            self.S = ds.layer["spliced"][:, :]
            self.U = ds.layer["unspliced"][:, :]
            self.A = ds.layer["ambiguous"][:, :]
            self.initial_cell_size = self.S.sum(0)
            self.initial_Ucell_size = self.U.sum(0)
            self.ca = dict(ds.col_attrs.items())
            self.ra = dict(ds.row_attrs.items())
            ds.close()
        else:
            ds = loompy.connect(self.loom_filepath)
            self.raw_S = ds.layer["spliced"][:, :]
            self.raw_U = ds.layer["unspliced"][:, :]
            self.raw_A = ds.layer["ambiguous"][:, :]
            self.raw_initial_cell_size = self.raw_S.sum(0)
            self.raw_initial_Ucell_size = self.raw_U.sum(0)
            self.raw_ca = dict(ds.col_attrs.items())
            self.raw_ra = dict(ds.row_attrs.items())
            ds.close()


def scatter_viz(x: np.ndarray, y: np.ndarray, *args: Any, **kwargs: Any) -> Any:
    """A wrapper of scatter plot that guarantees that every point is visible in a very crowded scatterplot

    Args
    ----
    x: np.ndarray
        x axis coordinates
    y: np.ndarray
        y axis coordinates
    args and kwargs:
        positional and keyword arguments as in matplotplib.pyplot.scatter

    Returns
    -------
    Plots the graph and returns the axes object
    """
    ix_x_sort = np.argsort(x, kind="mergesort")
    ix_yx_sort = np.argsort(y[ix_x_sort], kind="mergesort")
    args_new = []
    kwargs_new = {}
    for arg in args:
        if type(arg) is np.ndarray:
            args_new.append(arg[ix_x_sort][ix_yx_sort])
        else:
            args_new.append(arg)
    for karg, varg in kwargs.items():
        if type(varg) is np.ndarray:
            kwargs_new[karg] = varg[ix_x_sort][ix_yx_sort]
        else:
            kwargs_new[karg] = varg
    ax = plt.scatter(x[ix_x_sort][ix_yx_sort], y[ix_x_sort][ix_yx_sort], *args_new, **kwargs_new)
    return ax


def ixs_thatsort_a2b(a: np.ndarray, b: np.ndarray, check_content: bool=True) -> np.ndarray:
    "This is super duper magic sauce to make the order of one list to be like another"
    if check_content:
        assert len(np.intersect1d(a, b)) == len(a), f"The two arrays are not matching"
    return np.argsort(a)[np.argsort(np.argsort(b))]

colors20 = np.vstack((plt.cm.tab20b(np.linspace(0., 1, 20))[::2], plt.cm.tab20c(np.linspace(0, 1, 20))[1::2]))


def colormap_fun(x: np.ndarray) -> np.ndarray:
    return colors20[np.mod(x, 20)]


@jit("float64[:](float64[:], int32[:], int32[:], float64[:])", nopython=True)
def _scale_to_match_median(data: np.ndarray, indices: np.ndarray,
                           indptr: np.ndarray, genes_total: np.ndarray) -> np.ndarray:
    # Helper function that operates directly on the .data array of a sparse matrix object
    new_data = np.zeros(data.shape)
    # Loop through the columns
    for i in range(genes_total.shape[0]):
        # Retrieve the values
        non_zero_genes_total = genes_total[indices[indptr[i]:indptr[i + 1]]]
        # Find the normalization factor
        w = np.minimum(1, np.median(non_zero_genes_total) / non_zero_genes_total)
        new_data[indptr[i]:indptr[i + 1]] = w * data[indptr[i]:indptr[i + 1]]
    return new_data


@jit(nopython=True)
def numba_random_seed(value: int) -> None:
    """Same as np.random.seed but for numba"""
    np.random.seed(value)


@jit(nopython=True)
def permute_rows_nsign(A: np.ndarray) -> None:
    """Permute in place the entries and randomly switch the sign for each row of a matrix independently.
    """
    plmi = np.array([+1, -1])
    for i in range(A.shape[0]):
        np.random.shuffle(A[i, :])
        A[i, :] = A[i, :] * np.random.choice(plmi, size=A.shape[1])


def scale_to_match_median(sparse_matrix: sparse.csr_matrix, genes_total: np.ndarray) -> sparse.csr_matrix:
    """Normalize contribution of different neighbor genes to match the median totals
    
    Arguments
    ---------
    sparse_matrix: sparse.csr_matrix
        weights matrix
    
    genes_total: sparse.csr_matrix shape=(sparse_matrix.shape[0])
        array of the total molecules detected for each gene
    
    Returns
    -------
    knn_weights: sparse.csr_matrix
        sparse_matrix after the normalization
    
    # NOTE, since the use I made of this later I could have changed sparse_matrix in place
    """
    newdata = _scale_to_match_median(sparse_matrix.data, sparse_matrix.indices, sparse_matrix.indptr, genes_total)
    return sparse.csc_matrix((newdata,
                              sparse_matrix.indices,
                              sparse_matrix.indptr),
                             shape=sparse_matrix.shape,
                             copy=True)


def gaussian_kernel(X: np.ndarray, mu: float=0, sigma: float=1) -> np.ndarray:
    """Compute gaussian kernel"""
    return np.exp(-(X - mu)**2 / (2 * sigma**2)) / np.sqrt(2 * np.pi * sigma**2)


def load_velocyto_hdf5(filename: str) -> VelocytoLoom:
    """Loads a Velocyto loom object from an hdf5 file

    Arguments
    ---------
    filename: str
        The name of the serialized file

    Returns
    -------
    A VelocytoLoom object

    Note
    ----
    The hdf5 file must have been created using ``VelocytoLoom.to_hdf5`` or the ``dump_hdf5`` function
    """
    return load_hdf5(filename, obj_class=VelocytoLoom)
