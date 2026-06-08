"""Self-contained port of `StemRobustLTS.fit_circles_rlts` for FORTLS.

Public entry point: ``fit_circles_macif(points, ...)`` — see the function
docstring. Replicates the pipeline of ``main_nurunnabi_slice_stem`` in
``tools/tree_diameters.py`` without any dependency on the custom packages
``pointools``, ``itdtools``, or ``dendromatics``.

Dependencies: numpy, scipy, scikit-learn.
"""

from __future__ import annotations

import math
from typing import Optional

import numpy as np
from scipy.spatial import KDTree
from sklearn.linear_model import LinearRegression


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------

def _shift(points: np.ndarray, vector: np.ndarray) -> np.ndarray:
    return points + vector


def _pca(points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    pts = points.astype(np.float64, copy=False)
    pts = pts - pts.mean(axis=0)
    cov = (pts.T @ pts) / pts.shape[0]
    svd = np.linalg.svd(cov)
    return svd.S, svd.U


def _pca_rotation_matrix(points: np.ndarray) -> np.ndarray:
    _, eig_vec = _pca(points)
    return eig_vec[:, [1, 2, 0]]


def _get_strip(points: np.ndarray, height: float, width: float) -> np.ndarray:
    z_min = max(height - width / 2, 0.0)
    z_max = z_min + width
    z = points[:, 2]
    return points[(z >= z_min) & (z <= z_max)]


def _align_pca(points: np.ndarray, min_h: float, max_h: float) -> np.ndarray:
    mid = (max_h - min_h) / 2 + min_h
    half_w = (max_h - min_h) / 2
    sub = _get_strip(points, mid, 2 * half_w)
    sub = sub - sub.mean(axis=0)
    rot = _pca_rotation_matrix(sub)

    min_anchor = points[np.argmin(points[:, 2])]
    max_anchor = points[np.argmax(points[:, 2])]
    post_min = float(np.dot(min_anchor, rot)[-1])
    post_max = float(np.dot(max_anchor, rot)[-1])
    assert post_min != post_max
    rot = rot * (1 - 2 * (post_min > post_max))
    return points @ rot


def _count_decimals(value: float) -> int:
    s = f"{value:.10g}"
    if "." in s:
        return len(s.split(".")[1])
    return 0


def _get_heights(min_h: float, max_h: float, step: float, width: float) -> np.ndarray:
    dec = max(_count_decimals(v) for v in (min_h, max_h, width, step))
    max_h = max_h - (max_h - min_h) % step
    return np.round(np.arange(min_h, max_h, step), dec)


def _slice_stem(
    points: np.ndarray,
    heights: np.ndarray,
    width: float,
    min_points: int = 10,
) -> tuple[list[np.ndarray], np.ndarray]:
    slices: list[np.ndarray] = []
    kept: list[float] = []
    # Sort once and binary-search per slice instead of full mask per slice.
    order = np.argsort(points[:, 2])
    sorted_pts = points[order]
    z_sorted = sorted_pts[:, 2]
    half = width / 2
    for h in heights:
        z_min = max(float(h) - half, 0.0)
        z_max = z_min + width
        lo = np.searchsorted(z_sorted, z_min, side="left")
        hi = np.searchsorted(z_sorted, z_max, side="right")
        if hi - lo == 0:
            continue
        if min_points is not None and (hi - lo) < min_points:
            continue
        slices.append(sorted_pts[lo:hi])
        kept.append(float(h))
    return slices, np.array(kept)


def _project(slc: np.ndarray) -> np.ndarray:
    # Projected stem points are essentially never bit-exact duplicates; skip
    # the O(N log N) `np.unique` from the reference implementation.
    return np.ascontiguousarray(slc[:, :2])


# ---------------------------------------------------------------------------
# Outlier / robustness helpers
# ---------------------------------------------------------------------------

def _compute_iteration_number(h_0: int = 3, pr: float = 0.9999, eps: float = 0.5) -> int:
    it = np.log(1 - pr) / np.log(1 - (1 - eps) ** h_0)
    return int(np.round(it))


def _median_absolute_deviation(arr: np.ndarray) -> float:
    return 1.4826 * float(np.median(np.abs(arr - np.median(arr))))


def _robust_zscore(arr: np.ndarray) -> np.ndarray:
    return np.abs(arr - np.median(arr)) / (_median_absolute_deviation(arr) + 1e-8)


def _inner_circle(X: np.ndarray, Y: np.ndarray, xc: float, yc: float, r: float, inner_r: float) -> int:
    distance = np.sqrt((X - xc) ** 2 + (Y - yc) ** 2)
    return int(np.sum(distance < r * inner_r))


def _circle_intersection(distance: float, r_a: float, r_b: float, eps: float = 1e-5) -> float:
    if distance >= (r_a + r_b):
        return 0.0
    if (distance <= abs(r_a - r_b)) or ((distance <= eps) and (abs(r_a - r_b) <= eps)):
        return math.pi * min(r_a, r_b) ** 2
    d_a = (r_a ** 2 - r_b ** 2 + distance ** 2) / (2 * distance)
    d_b = distance - d_a
    term_1 = r_a ** 2 * math.acos(d_a / r_a) - d_a * math.sqrt(r_a ** 2 - d_a ** 2)
    term_2 = r_b ** 2 * math.acos(d_b / r_b) - d_b * math.sqrt(r_b ** 2 - d_b ** 2)
    return term_1 + term_2


def _circle_iou(circle_a: np.ndarray, circle_b: np.ndarray) -> float:
    r_a = float(circle_a[2])
    r_b = float(circle_b[2])
    dist = float(np.linalg.norm(np.asarray(circle_a[:2]) - np.asarray(circle_b[:2])))
    inter = _circle_intersection(dist, r_a, r_b)
    union = math.pi * (r_a ** 2 + r_b ** 2) - inter
    return inter / union if union > 0 else 0.0


def _batched_circle_iou(
    centers: np.ndarray, radii: np.ndarray, target: np.ndarray, eps: float = 1e-5,
) -> np.ndarray:
    """Vectorized IoU of (B,) candidate circles against a single ``target``
    circle ``[x, y, r]``. Returns shape ``(B,)``."""
    tx, ty, tr = float(target[0]), float(target[1]), float(target[2])
    dist = np.sqrt((centers[:, 0] - tx) ** 2 + (centers[:, 1] - ty) ** 2)
    inter = np.zeros_like(dist)
    # Branch 1: circles do not overlap.
    no_overlap = dist >= (radii + tr)
    # Branch 2: one fully contains the other.
    contained = (dist <= np.abs(radii - tr)) | (
        (dist <= eps) & (np.abs(radii - tr) <= eps)
    )
    inter = np.where(contained, math.pi * np.minimum(radii, tr) ** 2, inter)
    # Branch 3: proper partial overlap.
    partial = ~(no_overlap | contained)
    if partial.any():
        d = dist[partial]
        ra = radii[partial]
        d_a = (ra * ra - tr * tr + d * d) / (2.0 * d)
        d_b = d - d_a
        # Clamp arguments of acos/sqrt to avoid numerical NaNs.
        ca = np.clip(d_a / ra, -1.0, 1.0)
        cb = np.clip(d_b / tr, -1.0, 1.0)
        ta = ra * ra * np.arccos(ca) - d_a * np.sqrt(np.maximum(ra * ra - d_a * d_a, 0.0))
        tb = tr * tr * np.arccos(cb) - d_b * np.sqrt(np.maximum(tr * tr - d_b * d_b, 0.0))
        inter[partial] = ta + tb
    union = math.pi * (radii * radii + tr * tr) - inter
    return np.where(union > 0, inter / union, 0.0)


# ---------------------------------------------------------------------------
# Vectorized Hyper Least-Squares circle fit (Kanatani & Rangarajan 2011)
# ---------------------------------------------------------------------------

def _hyperLSQ_batch(coords: np.ndarray, iter_max: int = 12) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Vectorized Hyper-LSQ. ``coords`` is ``(B, P, 2)`` → returns
    ``(xc, yc, r)`` each of shape ``(B,)``.

    Newton-Raphson runs for a fixed ``iter_max`` (12 is more than enough in
    float64 — convergence is typically reached in 5-7 steps). Degenerate
    batch elements yield ``nan``/``inf`` here, filtered by the caller.
    """
    x = coords[..., 0]
    y = coords[..., 1]

    mx = x.mean(axis=-1, keepdims=True)
    my = y.mean(axis=-1, keepdims=True)
    Xi = x - mx
    Yi = y - my
    Zi = Xi * Xi + Yi * Yi

    Mxy = (Xi * Yi).mean(axis=-1)
    Mxx = (Xi * Xi).mean(axis=-1)
    Myy = (Yi * Yi).mean(axis=-1)
    Mxz = (Xi * Zi).mean(axis=-1)
    Myz = (Yi * Zi).mean(axis=-1)
    Mzz = (Zi * Zi).mean(axis=-1)

    Mz = Mxx + Myy
    Cov_xy = Mxx * Myy - Mxy * Mxy
    Var_z = Mzz - Mz * Mz

    A2 = 4 * Cov_xy - 3 * Mz * Mz - Mzz
    A1 = Var_z * Mz + 4.0 * Cov_xy * Mz - Mxz * Mxz - Myz * Myz
    A0 = (
        Mxz * (Mxz * Myy - Myz * Mxy)
        + Myz * (Myz * Mxx - Mxz * Mxy)
        - Var_z * Cov_xy
    )
    A22 = A2 + A2

    # Newton-Raphson root finding, lockstep across the batch. We accept that
    # ``Dy == 0`` for some entries produces inf/nan; degenerate fits are
    # detected by ``np.isfinite`` in the caller and dropped.
    with np.errstate(divide="ignore", invalid="ignore"):
        X = np.zeros_like(A0)
        for _ in range(iter_max):
            Dy = A1 + X * (A22 + 16.0 * X * X)
            Y = A0 + X * (A1 + X * (A2 + 4.0 * X * X))
            X = X - Y / Dy

        det = X * X - X * Mz + Cov_xy
        inv2det = 0.5 / det
        Xcenter = (Mxz * (Myy - X) - Myz * Mxy) * inv2det
        Ycenter = (Myz * (Mxx - X) - Mxz * Mxy) * inv2det

    xc = Xcenter + mx[..., 0]
    yc = Ycenter + my[..., 0]
    r = np.sqrt(np.abs(Xcenter * Xcenter + Ycenter * Ycenter + Mz))
    return xc, yc, r


def _three_point_circle_batch(triples: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Closed-form circle through 3 points, vectorized over a batch.

    ``triples`` is ``(B, 3, 2)``. Returns ``(xc, yc, r)`` each ``(B,)``.
    Collinear triples yield ``inf`` — callers should filter beforehand.
    """
    p1 = triples[:, 0, :]
    p2 = triples[:, 1, :]
    p3 = triples[:, 2, :]
    x1, y1 = p1[:, 0], p1[:, 1]
    x2, y2 = p2[:, 0], p2[:, 1]
    x3, y3 = p3[:, 0], p3[:, 1]

    a = 2.0 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    sq1 = x1 * x1 + y1 * y1
    sq2 = x2 * x2 + y2 * y2
    sq3 = x3 * x3 + y3 * y3
    safe = a != 0.0
    inv_a = np.where(safe, 1.0 / a, 0.0)
    xc = (sq1 * (y2 - y3) + sq2 * (y3 - y1) + sq3 * (y1 - y2)) * inv_a
    yc = (sq1 * (x3 - x2) + sq2 * (x1 - x3) + sq3 * (x2 - x1)) * inv_a
    r = np.sqrt((x1 - xc) ** 2 + (y1 - yc) ** 2)
    return xc, yc, r


def _sample_non_collinear_triples(
    n_points: int, n_triples: int, points: np.ndarray, rng: np.random.Generator,
    tol: float = 1e-6, max_attempts: int = 6,
) -> np.ndarray:
    """Draw ``n_triples`` triples of distinct indices from ``range(n_points)``
    such that the three picked points are not collinear (within ``tol``).
    Vectorized rejection sampling — typically 1-2 passes.
    """
    out = np.empty((n_triples, 3), dtype=np.int64)
    needed = np.arange(n_triples)
    for _ in range(max_attempts):
        if needed.size == 0:
            break
        m = needed.size
        idx = rng.integers(0, n_points, size=(m, 3))
        # Disallow repeats within a triple.
        valid = (idx[:, 0] != idx[:, 1]) & (idx[:, 0] != idx[:, 2]) & (idx[:, 1] != idx[:, 2])
        # Non-collinearity: cross-product magnitude in 2D.
        p = points[idx]
        v1 = p[:, 1] - p[:, 0]
        v2 = p[:, 2] - p[:, 0]
        cross = v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0]
        valid &= np.abs(cross) > tol
        accepted = needed[valid]
        out[accepted] = idx[valid]
        needed = needed[~valid]
    if needed.size:
        # Fallback: deterministic enumeration (very rare for stem clouds).
        for k in needed:
            for a in range(n_points - 2):
                for b in range(a + 1, n_points - 1):
                    for c in range(b + 1, n_points):
                        v1 = points[b] - points[a]
                        v2 = points[c] - points[a]
                        if abs(v1[0] * v2[1] - v1[1] * v2[0]) > tol:
                            out[k] = (a, b, c)
                            break
                    else:
                        continue
                    break
                else:
                    continue
                break
    return out


# ---------------------------------------------------------------------------
# Vectorized Robust LTS circle fitters
# ---------------------------------------------------------------------------

def _rlts_fit(
    points: np.ndarray,
    n_iter: int,
    h_size: int,
    rng: np.random.Generator,
) -> tuple[float, float, float, float, np.ndarray]:
    """Vectorized Robust Least-Trimmed-Squares circle fit (replacement for
    ``RobustLTS.compute_rlts``). Returns ``(xc, yc, r, residual, idx)``.
    """
    N = points.shape[0]
    h_size = max(3, min(h_size, N))

    # Mean-center for numerical conditioning, exactly like the reference.
    shift = points.mean(0)
    pts = points - shift

    # Draw n_iter non-collinear triples and fit closed-form circles.
    triples = _sample_non_collinear_triples(N, n_iter, pts, rng)
    sub_pts = pts[triples]  # (n_iter, 3, 2)
    xc0, yc0, r0 = _three_point_circle_batch(sub_pts)

    # Residuals from each initial circle to ALL points: (n_iter, N).
    diff_x = pts[None, :, 0] - xc0[:, None]
    diff_y = pts[None, :, 1] - yc0[:, None]
    dist = np.sqrt(diff_x * diff_x + diff_y * diff_y)
    res0 = (dist - r0[:, None]) ** 2

    # Top h_size indices per iter via argpartition.
    if h_size < N:
        h_idx0 = np.argpartition(res0, h_size - 1, axis=-1)[:, :h_size]
    else:
        h_idx0 = np.broadcast_to(np.arange(N), (n_iter, N))[:, :h_size]

    # Refit on the h_size subset with batched Hyper-LSQ.
    h_pts = pts[h_idx0]  # (n_iter, h_size, 2)
    xc1, yc1, r1 = _hyperLSQ_batch(h_pts)

    # Final residuals + take top h_size per iter again.
    diff_x = pts[None, :, 0] - xc1[:, None]
    diff_y = pts[None, :, 1] - yc1[:, None]
    dist = np.sqrt(diff_x * diff_x + diff_y * diff_y)
    res1 = (dist - r1[:, None]) ** 2
    if h_size < N:
        h_idx1 = np.argpartition(res1, h_size - 1, axis=-1)[:, :h_size]
    else:
        h_idx1 = np.broadcast_to(np.arange(N), (n_iter, N))[:, :h_size]
    h_res = np.take_along_axis(res1, h_idx1, axis=-1).sum(axis=-1)

    # Drop NaN/Inf fits (degenerate iterations).
    finite = np.isfinite(h_res) & np.isfinite(xc1) & np.isfinite(yc1) & np.isfinite(r1)
    if not finite.any():
        # Worst-case fallback: best of whatever we have.
        h_res[~finite] = np.inf
        best = int(np.argmin(h_res))
    else:
        h_res = np.where(finite, h_res, np.inf)
        best = int(np.argmin(h_res))

    best_xc = float(xc1[best] + shift[0])
    best_yc = float(yc1[best] + shift[1])
    best_r = float(r1[best])
    best_residual = float(h_res[best] / h_size)
    return best_xc, best_yc, best_r, best_residual, h_idx1[best]


def _iou_rlts_fit(
    points: np.ndarray,
    prev_circle: np.ndarray,
    n_iter: int,
    h_size: int,
    top_k: int,
    rng: np.random.Generator,
) -> tuple[float, float, float, float, np.ndarray]:
    """Vectorized IoU-scored RLTS (port of ``IoURobustLTS.compute_rlts``)."""
    N = points.shape[0]
    h_size = max(3, min(h_size, N))

    shift = points.mean(0)
    pts = points - shift

    triples = _sample_non_collinear_triples(N, n_iter, pts, rng)
    sub_pts = pts[triples]
    xc0, yc0, r0 = _three_point_circle_batch(sub_pts)

    diff_x = pts[None, :, 0] - xc0[:, None]
    diff_y = pts[None, :, 1] - yc0[:, None]
    dist = np.sqrt(diff_x * diff_x + diff_y * diff_y)
    res0 = (dist - r0[:, None]) ** 2
    if h_size < N:
        h_idx0 = np.argpartition(res0, h_size - 1, axis=-1)[:, :h_size]
    else:
        h_idx0 = np.broadcast_to(np.arange(N), (n_iter, N))[:, :h_size]

    h_pts = pts[h_idx0]
    xc1, yc1, r1 = _hyperLSQ_batch(h_pts)

    diff_x = pts[None, :, 0] - xc1[:, None]
    diff_y = pts[None, :, 1] - yc1[:, None]
    dist = np.sqrt(diff_x * diff_x + diff_y * diff_y)
    res1 = (dist - r1[:, None]) ** 2
    if h_size < N:
        h_idx1 = np.argpartition(res1, h_size - 1, axis=-1)[:, :h_size]
    else:
        h_idx1 = np.broadcast_to(np.arange(N), (n_iter, N))[:, :h_size]
    h_res = np.take_along_axis(res1, h_idx1, axis=-1).sum(axis=-1)
    finite = np.isfinite(h_res) & np.isfinite(xc1) & np.isfinite(yc1) & np.isfinite(r1)
    h_res = np.where(finite, h_res, np.inf)

    # Top-k by residual, then pick the one with the lowest IoU (matches the
    # original `IoURobustLTS.compute_rlts` semantics — `np.argmin` on the IoU
    # column of the top-k slice).
    k = min(top_k, n_iter)
    top_idx = np.argpartition(h_res, k - 1)[:k] if k < n_iter else np.arange(n_iter)
    cand_centers = np.column_stack([xc1[top_idx] + shift[0], yc1[top_idx] + shift[1]])
    cand_radii = r1[top_idx]
    ious = _batched_circle_iou(cand_centers, cand_radii, prev_circle[:3])
    best = int(top_idx[int(np.argmin(ious))])

    best_xc = float(xc1[best] + shift[0])
    best_yc = float(yc1[best] + shift[1])
    best_r = float(r1[best])
    best_residual = float(h_res[best] / h_size)
    return best_xc, best_yc, best_r, best_residual, h_idx1[best]


def _momentum_rlts_fit(
    points: np.ndarray,
    prev_circle: np.ndarray,
    n_iter: int,
    crop_n_iter: int,
    h_size: int,
    rng: np.random.Generator,
    min_points: int = 10,
    min_crop: float = 1.1,
    max_crop: float = 1.5,
    early_exit_iou: float = 0.85,
) -> tuple[float, float, float, float, np.ndarray]:
    """Momentum-corrected RLTS (port of ``MomentumRobustLTS.compute_rlts``).

    Adds an *early exit*: when the initial IoU-scored fit already meets
    ``early_exit_iou``, we skip the crop iterations.
    """
    N = points.shape[0]
    h_size = max(3, min(h_size, N))
    h_size = max(h_size, min(min_points, N))

    # Initial fit constrained by prev_circle via IoU scoring.
    xc, yc, r, res, idx = _iou_rlts_fit(points, prev_circle, n_iter, h_size, 10, rng)
    best_iou = _circle_iou(np.array([xc, yc, r]), np.asarray(prev_circle[:3], dtype=np.float64))
    best_xc, best_yc, best_r, best_residual, best_idx = xc, yc, r, res, idx

    if best_iou >= early_exit_iou:
        return best_xc, best_yc, best_r, best_residual, best_idx

    prev_center = np.asarray(prev_circle[:2], dtype=np.float64)
    prev_radius = float(prev_circle[2])
    kdtree = KDTree(points)

    max_rad = prev_radius * max_crop
    min_rad = prev_radius * min_crop
    if crop_n_iter <= 0:
        return best_xc, best_yc, best_r, best_residual, best_idx
    steps = np.arange(max_rad, min_rad, (min_rad - max_rad) / crop_n_iter)
    for r_crop in steps:
        cropped = kdtree.query_ball_point(prev_center, float(r_crop))
        if len(cropped) < min_points:
            break
        crop_pts = points[cropped]
        crop_h_size = max(3, int(crop_pts.shape[0] * h_size / N))
        crop_h_size = max(crop_h_size, min(min_points, crop_pts.shape[0]))
        xc2, yc2, r2, res2, idx2 = _rlts_fit(crop_pts, n_iter, crop_h_size, rng)
        iou = _circle_iou(np.array([xc2, yc2, r2]), np.asarray(prev_circle[:3], dtype=np.float64))
        if iou > best_iou:
            best_iou = iou
            best_xc, best_yc, best_r, best_residual = xc2, yc2, r2, res2
            best_idx = np.asarray(cropped)[idx2]
            if best_iou >= early_exit_iou:
                break
    return best_xc, best_yc, best_r, best_residual, best_idx


# ---------------------------------------------------------------------------
# Linear regression on (x, y, r) vs z — used in `find_best_section_global`
# and `correct_best_section`. Hot path is small (run a handful of times), so
# the simple sklearn-backed wrapper is fast enough.
# ---------------------------------------------------------------------------

class _StemLinearRegression:
    def __init__(self) -> None:
        self._x_model = LinearRegression()
        self._y_model = LinearRegression()
        self._r_model = LinearRegression()

    def fit(self, x: np.ndarray, y: np.ndarray) -> None:
        x_set, y_set, r_set = np.hsplit(y, 3)
        self._x_model.fit(x, x_set)
        self._y_model.fit(x, y_set)
        self._r_model.fit(x, r_set)

    def predict(self, x: np.ndarray) -> np.ndarray:
        return np.hstack([
            self._x_model.predict(x),
            self._y_model.predict(x),
            self._r_model.predict(x),
        ])

    def fit_predict(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        self.fit(x, y)
        return self.predict(x)

    @staticmethod
    def relative_squared_error(y_pred: np.ndarray, y_true: np.ndarray) -> np.ndarray:
        mean = y_true.mean(axis=0)
        mean_diff = y_pred - mean
        return (y_pred - y_true) ** 2 / ((mean_diff ** 2).sum(axis=0) + 1e-12)


def _robust_slr_best_model(
    points: np.ndarray, n_iter: int, h_prop: float, rng: np.random.Generator,
) -> _StemLinearRegression:
    """RANSAC-like robust 3-target linear regression (port of `RobustSLR`)."""
    size = points.shape[0]
    h_size = max(2, int(size * h_prop))
    best_res = np.inf
    best_model: Optional[_StemLinearRegression] = None
    z_col = points[:, -1].reshape(-1, 1)
    xyr = points[:, :3]
    for _ in range(n_iter):
        sub = rng.integers(0, size, size=2)
        if sub[0] == sub[1]:
            continue
        m = _StemLinearRegression()
        m.fit(z_col[sub], xyr[sub])
        preds = m.predict(z_col)
        res = _StemLinearRegression.relative_squared_error(preds, xyr)
        idx = np.argsort(res.mean(1))[:h_size]
        h_res = float(res[idx].sum())
        if h_res < best_res:
            best_res = h_res
            best_model = m
    return best_model


# ---------------------------------------------------------------------------
# Stem-level orchestration (port of StemRobustLTS)
# ---------------------------------------------------------------------------

def _slice_array_to_section(arr: np.ndarray, n: int, sample_size: int) -> np.ndarray:
    return arr[n:n + sample_size].reshape(sample_size, -1)


def _slice_list_to_section(lst: list[np.ndarray], n: int, sample_size: int) -> list[np.ndarray]:
    return lst[n:n + sample_size]


def _update_section(iterable, new_value, n: int, sample_size: int) -> None:
    iterable[n:n + sample_size] = new_value


def _correct_slice(slc: np.ndarray, center: np.ndarray, radius: float) -> np.ndarray:
    kdtree = KDTree(slc)
    outer = kdtree.query_ball_point(center, radius)
    return slc[np.array(outer)]


def _check_best_circles(
    circles: np.ndarray, preds: np.ndarray, max_inner: int = 5, inner_r: float = 0.5,
) -> np.ndarray:
    diff = np.abs(circles[:, :3] - preds)
    xy_out = diff[:, :2] >= np.median(circles[:, 2]) * inner_r
    r_out = _robust_zscore(circles[:, 2]) >= 2.5
    inner_out = circles[:, -2] > max_inner
    return xy_out.any(1) * r_out * inner_out


def _compute_slice_iou(buffer: np.ndarray) -> np.ndarray:
    return np.array([
        _circle_iou(buffer[i, ...], buffer[i + 1, ...])
        for i in range(buffer.shape[0] - 1)
    ])


def _find_best_section_global(
    circles: np.ndarray, sample_size: int, rng: np.random.Generator,
) -> int:
    it = _compute_iteration_number(h_0=2, pr=0.99999, eps=0.7)
    model = _robust_slr_best_model(circles, it, 0.5, rng)
    preds = model.predict(circles[:, -1].reshape(-1, 1))
    res = _StemLinearRegression.relative_squared_error(preds, circles[:, :3])

    mean_res = res.mean(1)
    # Sliding-window mean via cumulative sum (O(N) vs O(N·sample_size)).
    cs = np.concatenate(([0.0], np.cumsum(mean_res)))
    window_means = (cs[sample_size:] - cs[:-sample_size]) / sample_size
    return int(np.argmin(window_means))


def _fit_single_circle(
    pts: np.ndarray, n_iter: int, min_points: int, inner_r: float, rng: np.random.Generator,
) -> tuple:
    h_size = max(pts.shape[0] // 2, min(min_points, pts.shape[0]))
    xc, yc, r, res, _ = _rlts_fit(pts, n_iter, h_size, rng)
    inner = _inner_circle(pts[:, 0], pts[:, 1], xc, yc, r, inner_r)
    return xc, yc, r, res, inner


def _fit_circles_initial(
    slices: list[np.ndarray], heights: np.ndarray, n_iter: int, inner_r: float,
    min_points: int, rng: np.random.Generator,
) -> np.ndarray:
    rows = [_fit_single_circle(pts, n_iter, min_points, inner_r, rng) for pts in slices]
    arr = np.array(rows, dtype=np.float64)
    return np.hstack([arr, heights.reshape(-1, 1)])


def _correct_best_section(
    slices: list[np.ndarray],
    circles: np.ndarray,
    n_iter: int,
    rng: np.random.Generator,
    inner_r: float = 0.5,
    max_inner: int = 5,
    outer_r: float = 1.5,
) -> None:
    heights = circles[:, -1].reshape(-1, 1)
    model = _StemLinearRegression()
    model.fit(heights, circles[:, :3])
    preds = model.predict(heights)

    mask = _check_best_circles(circles, preds, max_inner, inner_r)
    for i in np.where(mask)[0]:
        new_slc = _correct_slice(slices[i], preds[i, :2], preds[i, 2] * outer_r)
        if new_slc.shape[0] < 3:
            continue
        new_circle = _fit_single_circle(new_slc, n_iter, 10, inner_r, rng)
        circles[i, :-1] = np.array(new_circle, dtype=circles.dtype)
        slices[i] = new_slc
        if i == (mask.size - 1):
            return
        model = _StemLinearRegression()
        preds = model.fit_predict(circles[:, -1].reshape(-1, 1), circles[:, :3])


def _correct_section_downwards(
    circles: np.ndarray,
    buffer_slices: list[np.ndarray],
    rng: np.random.Generator,
    n_iter: int = 3,
    inner_n_iter: Optional[int] = None,
    min_crop: float = 1.1,
    max_crop: float = 1.5,
    inner_r: float = 0.5,
) -> None:
    if inner_n_iter is None:
        inner_n_iter = _compute_iteration_number()

    buffer = circles[1:len(buffer_slices), :]
    if buffer.shape[0] < 2:
        return
    slc_iou = float(_compute_slice_iou(buffer).mean())

    circle_0 = circles[0, :]
    circle_1 = circles[1, :]
    init_iou = _circle_iou(circle_0, circle_1)

    best_circle = circle_0[:-1].copy()
    slice_0 = buffer_slices[0]
    best_inn = _inner_circle(
        slice_0[:, 0], slice_0[:, 1],
        float(circle_0[0]), float(circle_0[1]), float(circle_0[2]), inner_r,
    )

    if not ((init_iou < 0.75) or (init_iou < slc_iou * 0.8) or (circle_0[-2] > 0)):
        return

    best_iou = init_iou
    # Incremental buffer accumulation: list-append + single vstack at use time
    # avoids the O(N^2) `np.vstack(buffer_slices[0:n])` + `np.unique` per step.
    acc = [buffer_slices[0]]
    slice_tmp = buffer_slices[0]

    for n in range(2, len(buffer_slices) + 1):
        if slice_tmp.shape[0] < 3:
            break
        h_size = slice_tmp.shape[0] // 2
        xc, yc, r, res, _ = _momentum_rlts_fit(
            slice_tmp, circle_1, inner_n_iter, n_iter, h_size, rng,
            min_crop=min_crop, max_crop=max_crop,
        )
        inn_check = _inner_circle(slice_0[:, 0], slice_0[:, 1], xc, yc, r, inner_r)
        new_circle = np.array([xc, yc, r, res, inn_check])
        temp_iou = _circle_iou(new_circle, circle_1)
        if temp_iou >= best_iou:
            if best_iou >= 0.75:
                if inn_check <= best_inn:
                    best_iou = temp_iou
                    best_circle = new_circle
                    best_inn = inn_check
            else:
                best_iou = temp_iou
                best_circle = new_circle
                best_inn = inn_check
        # Early exit once we have a clearly good correction — saves the rest
        # of the buffer-accumulation loop on easy slices.
        if best_iou >= 0.9:
            break
        acc.append(buffer_slices[n - 1])
        slice_tmp = np.vstack(acc)

    circles[0, :-1] = best_circle


def _correct_sections(
    slices: list[np.ndarray],
    circles: np.ndarray,
    n_iter: int,
    sample_size: int,
    min_crop: float,
    max_crop: float,
    inner_r: float,
    rng: np.random.Generator,
) -> None:
    if circles.shape[0] < sample_size:
        return

    best_section = _find_best_section_global(circles, sample_size, rng)

    sub_slices = _slice_list_to_section(slices, best_section, sample_size)
    sub_circles = _slice_array_to_section(circles, best_section, sample_size)
    _correct_best_section(sub_slices, sub_circles, n_iter, rng, inner_r)
    _update_section(circles, sub_circles, best_section, sample_size)

    inner_n_iter = _compute_iteration_number()
    for new_section in range(best_section - 1, -1, -1):
        if new_section + sample_size + 1 > len(slices):
            continue
        sub_slices = _slice_list_to_section(slices, new_section, sample_size + 1)
        sub_circles = _slice_array_to_section(circles, new_section, sample_size + 1)
        _correct_section_downwards(
            sub_circles, sub_slices, rng,
            inner_n_iter=inner_n_iter, min_crop=min_crop, max_crop=max_crop, inner_r=inner_r,
        )
        _update_section(circles, sub_circles, new_section, sample_size + 1)

    for new_section in range(best_section, circles.shape[0] - sample_size):
        sub_slices = _slice_list_to_section(slices, new_section, sample_size + 1)
        sub_circles = _slice_array_to_section(circles, new_section, sample_size + 1)
        inv_sub_slices = sub_slices[::-1]
        inv_sub_circles = sub_circles[::-1]
        _correct_section_downwards(
            inv_sub_circles, inv_sub_slices, rng,
            inner_n_iter=inner_n_iter, min_crop=min_crop, max_crop=max_crop, inner_r=inner_r,
        )
        _update_section(circles, inv_sub_circles[::-1], new_section, sample_size + 1)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def fit_circles_macif(
    points: np.ndarray,
    *,
    step: float = 0.035,
    width: float = 0.07,
    sample_size: int = 10,
    min_crop: float = 1.2,
    max_crop: float = 1.75,
    inner_r: float = 0.5,
    min_height: float = 0.5,
    max_height: Optional[float] = None,
    min_points: int = 10,
    align: bool = True,
    h_0: int = 3,
    pr: float = 0.9999,
    eps: float = 0.5,
    seed: Optional[int] = None,
) -> np.ndarray:
    """Fit circles along a tree stem using the Nurunnabi RLTS pipeline.

    Parameters
    ----------
    points : (N, 3) float array
        Raw stem point cloud (columns x, y, z, any units consistent with `step`,
        `width`, `min_height`, `max_height`).
    step : float
        Vertical spacing between consecutive slice centres.
    width : float
        Vertical thickness of each slice.
    sample_size : int
        Number of consecutive slices that form a section in the global anchor
        search and correction loops.
    min_crop, max_crop : float
        Radius-multiplier bounds used by the momentum correction.
    inner_r : float
        Ratio of the inner radius to the fitted radius for the inner-point count.
    min_height, max_height : float
        Vertical extent of the analysis. `max_height=None` uses the cloud's z
        maximum.
    min_points : int
        Minimum number of points needed to keep a slice and to fit a circle.
    align : bool
        Apply PCA-based verticalisation to the stem before slicing.
    h_0, pr, eps : int, float, float
        Sampling parameters for the iteration count (see Nurunnabi et al.).
    seed : int | None
        Seed for the random sampler used during RLTS fitting.

    Returns
    -------
    (M, 6) float array with columns ``[x, y, radius, residual, inner, z]``.
    """
    rng = np.random.default_rng(seed)

    points = np.asarray(points, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(f"`points` must be of shape (N, 3); got {points.shape}")

    pts = _shift(points, -points.min(axis=0))
    if align:
        pts = _align_pca(pts, float(pts[:, 2].min()), float(pts[:, 2].max()))
        pts = _shift(pts, -pts.min(axis=0))

    z_max_cloud = float(pts[:, 2].max())
    z_max = z_max_cloud if max_height is None else min(float(max_height), z_max_cloud)

    heights = _get_heights(min_height, z_max, step, width)
    slices_3d, heights = _slice_stem(pts, heights, width, min_points)
    if len(slices_3d) == 0:
        return np.empty((0, 6), dtype=np.float64)

    slices_2d = [_project(s) for s in slices_3d]

    n_iter = _compute_iteration_number(h_0=h_0, pr=pr, eps=eps)
    circles = _fit_circles_initial(slices_2d, heights, n_iter, inner_r, min_points, rng)
    _correct_sections(
        slices_2d, circles, n_iter, sample_size, min_crop, max_crop, inner_r, rng,
    )
    return circles
