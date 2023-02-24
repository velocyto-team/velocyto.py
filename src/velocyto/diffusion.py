from typing import Any

import numpy as np
from numpy.random import default_rng
from scipy import sparse
from scipy.stats import norm
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize


class Diffusion:
    def __init__(self) -> None:
        pass

    def compute_transition_matrix2(
        self, x0: np.ndarray, v: np.ndarray, sigma: float = 0.0, reverse: bool = False
    ) -> sparse.csr_matrix:
        """
        Compute a right-stochastic matrix representing transition probabilities from each node

        Args:
                x0         Embedding positions (n_cells, n_dims)
                v          Velocities on the embedding (n_cells, n_dims
                sigma      The kernel size
                reverse    Compute the reverse transition matrix (for backwards diffusion)

        Remarks:
                Computes a Markov transition matrix for the KNN graph. The probability of transition along an edge
                is proportional to the scalar projection of the velocity vector onto that edge, times the reciprocal
                of the edge length. Edges that get negative scalar projections are clipped to zero and the total
                non-zero outgoing edges are normalized to a sum of 1.0.
        """
        n_cells = x0.shape[0]
        n_neighbors = 20

        # project into the past or future
        x1 = x0 - v if reverse else x0 + v
        # Find nearest neighbors
        nn = NearestNeighbors(n_neighbors=n_neighbors, algorithm="auto", n_jobs=4)
        (dists, nearest) = nn.fit(x0).kneighbors(x1)  # These are shaped (n_cells, n_neighbors), but we flatten them
        dists = dists.reshape(n_cells * n_neighbors)
        nearest = nearest.reshape(n_cells * n_neighbors)

        # Calculate transition probabilities
        probs = norm.pdf(dists, 0, sigma)

        # Make a sparse transition matrix
        cells = np.repeat(np.arange(n_cells), n_neighbors)
        tr = sparse.coo_matrix((probs, (cells, nearest)), shape=(n_cells, n_cells))
        # Make it right stochastic
        tr = normalize(tr.tocsr(), axis=1, norm="l1")
        return tr

    def compute_transition_matrix(
        self,
        knn: sparse.coo_matrix,
        x: np.ndarray,
        v: np.ndarray,
        epsilon: float = 0.0,
        reverse: bool = False,
    ) -> sparse.csr_matrix:
        """
        Compute a right-stochastic matrix representing transition probabilities from each node

        Args:
                knn        KNN graph (n_cells, n_cells)
                x          Embedding positions (n_cells, n_dims)
                v          Velocities on the embedding (n_cells, n_dims)
                reverse    Compute the reverse transition matrix (for backwards diffusion)

        Remarks:
                Computes a Markov transition matrix for the KNN graph. The probability of transition along an edge
                is proportional to the scalar projection of the velocity vector onto that edge, times the reciprocal
                of the edge length. Edges that get negative scalar projections are clipped to zero and the total
                non-zero outgoing edges are normalized to a sum of 1.0.
        """
        # vertices for each edge
        knn = knn.tocoo()
        (v0, v1) = (knn.row, knn.col)

        # calculate edge unit vectors
        uv = x[v1] - x[v0]  # Vector corresponding to an edge from v0 to v1, shape (n_edges, n_dims)
        norms = np.linalg.norm(uv, axis=1)
        uv = uv / norms[:, None]  # Convert to unit vector

        # Project the velocity vectors onto edges, and clip to zero
        scalar_projection = np.array([a.dot(b) for a, b in zip(v[v0], uv, strict=True)])  # Shape: (n_edges)
        if reverse:
            scalar_projection = -scalar_projection
        scalar_projection += epsilon
        # scalar_projection += scalar_projection.min()
        np.clip(scalar_projection, a_min=0, a_max=None, out=scalar_projection)

        # Calculate transition probabilities
        p = scalar_projection * (1 / norms)
        return normalize(sparse.coo_matrix((p, (v0, v1))).tocsr(), axis=1, norm="l1")

    def diffuse(
        self,
        x: np.ndarray,
        tr: sparse.csr_matrix,
        n_steps: int = 10,
        mode: str = "path_integral",
    ) -> Any:
        result = np.zeros(x.shape)
        if mode == "path_integral":
            x = sparse.csr_matrix(x / x.sum())
            for _ix in range(n_steps):
                x = x.dot(tr)
                result = result + x
            return result
        elif mode == "time_evolution":
            x = sparse.csr_matrix(x / x.sum())
            for _ix in range(n_steps):
                x = x.dot(tr)
            return x.toarray()
        elif mode == "map_trajectory":
            x = sparse.csr_matrix(x / x.sum())
            result = [np.argmax(x.toarray())]
            for _ix in range(n_steps):
                x = x.dot(tr)
                result.append(np.argmax(x.toarray()))
            return result
        elif mode == "frontier":
            x = sparse.csr_matrix(x / x.sum())
            result = [np.argmax(x.toarray())]
            for _ix in range(n_steps):
                x_next = x.dot(tr)
                result.append(np.argmax((x_next.toarray() + 1) / (x.toarray() + 1)))
                x = x_next
            return result
        elif mode == "trajectory":
            rng = default_rng()
            node = rng.choice(np.arange(x.shape[0]), p=x)
            trajectories = [node]
            for _ix in range(n_steps):
                x = np.zeros(tr.shape[0])
                x[node] = 1
                x = sparse.csr_matrix(x)
                x_next = x.dot(tr).toarray()
                x_next = normalize(x_next, norm="l1")[0]
                if x_next.sum() == 0:
                    x_next = x.toarray()[0]
                node = rng.choice(np.arange(x_next.shape[0]), p=x_next)
                trajectories.append(node)
                x = x_next
            return trajectories
