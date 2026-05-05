"""PERMANOVA (permutational multivariate ANOVA) on a distance matrix.

The pseudo-F statistic and sum-of-squares decomposition follow Anderson (2001)
and match the single-factor ``vegan::adonis2`` / ``adonis`` test and the
implementation in scikit-bio (BSD-licensed reference implementation).

References
----------
Anderson, M.J. (2001). A new method for non-parametric multivariate analysis of
variance. Austral Ecology 26, 32-46.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def _permanova_s_w(distance_sq: np.ndarray, grouping: np.ndarray) -> float:
    """Within-group SS component (``s_W``) matching skbio/vegan conventions."""
    n = distance_sq.shape[0]
    counts = np.bincount(grouping)
    s_w = 0.0
    for i in range(n):
        g = grouping[i]
        ng = counts[g]
        row_sum = 0.0
        for j in range(i + 1, n):
            if grouping[j] == g:
                row_sum += distance_sq[i, j]
        s_w += row_sum / ng
    return s_w


def permanova_pseudo_f_statistic(distances: np.ndarray, grouping: np.ndarray) -> float:
    """Compute the PERMANOVA pseudo-F for a symmetric distance matrix.

    Parameters
    ----------
    distances
        Square ``(n, n)`` distance matrix (zeros on diagonal).
    grouping
        Length-``n`` integer vector of group codes in ``0 .. n_groups-1``.

    Returns
    -------
    float
        Pseudo-F statistic (larger values indicate stronger separation).
    """
    d = np.asarray(distances, dtype=float)
    if d.ndim != 2 or d.shape[0] != d.shape[1]:
        raise ValueError("distances must be a square matrix.")
    n = d.shape[0]
    grp = np.asarray(grouping)
    if grp.shape[0] != n:
        raise ValueError("grouping length must match the distance matrix size.")
    _, grp = np.unique(grp, return_inverse=True)

    num_groups = int(np.unique(grp).size)
    if num_groups < 2:
        raise ValueError("PERMANOVA requires at least two groups.")
    if n <= num_groups:
        raise ValueError("PERMANOVA requires n_samples > n_groups.")

    d_sq = d * d
    s_t = d_sq.sum() / n / 2.0
    s_w = _permanova_s_w(d_sq, grp)
    s_a = s_t - s_w
    return (s_a / (num_groups - 1)) / (s_w / (n - num_groups))


def permanova_test(
    distances: np.ndarray,
    target: np.ndarray | pd.Series,
    permutations: int = 999,
    random_state: int | np.random.Generator | None = None,
) -> dict:
    """Permutation test for PERMANOVA.

    Parameters
    ----------
    distances
        Square distance matrix between samples.
    target
        Categorical group labels (one per row/column of ``distances``).
    permutations
        Number of random permutations of labels (excluding the observed layout).
        If ``0``, the p-value is ``nan``.
    random_state
        Seed or ``Generator`` for permutations.

    Returns
    -------
    dict
        Keys: ``pseudo_f``, ``p_value``, ``permutations``.
    """
    rng = np.random.default_rng(random_state)
    codes, _ = pd.factorize(np.asarray(target), sort=True)
    f_obs = permanova_pseudo_f_statistic(distances, codes.astype(int, copy=False))
    if permutations == 0:
        return {"pseudo_f": f_obs, "p_value": np.nan, "permutations": 0}

    ge = 1
    for _ in range(permutations):
        perm = rng.permutation(codes)
        if permanova_pseudo_f_statistic(distances, perm) >= f_obs - 1e-12:
            ge += 1
    p_value = ge / (permutations + 1)
    return {"pseudo_f": f_obs, "p_value": p_value, "permutations": permutations}


def validate_permanova_target(target: pd.Series) -> None:
    """Raise if ``target`` looks like a continuous outcome unsuitable for PERMANOVA."""
    t = pd.Series(target)
    n = len(t)
    n_u = t.nunique(dropna=False)
    if n_u < 2:
        raise ValueError("PERMANOVA requires at least two distinct groups.")
    if pd.api.types.is_numeric_dtype(t) and n_u > max(25, n // 2):
        raise ValueError(
            "PERMANOVA expects a categorical grouping; for continuous outcomes use method='knn' "
            "with task='regression'."
        )
