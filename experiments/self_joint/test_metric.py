"""Algebraic checks for the learned metric; no biological labels."""
import numpy as np
from threadpoolctl import threadpool_limits
from learn import transform, cross_validate


def main():
    rng = np.random.default_rng(42)
    dim = 140
    a = rng.normal(size=(dim, dim))
    within = a @ a.T + np.eye(dim)
    matrix = transform(within, "joint")
    assert np.allclose(matrix.T @ within @ matrix, np.eye(dim), atol=1e-10)
    delta = rng.normal(size=(100, dim))
    direct = np.einsum("ij,ij->i", delta @ np.linalg.inv(within), delta)
    projected = np.sum((delta @ matrix)**2, axis=1)
    assert np.allclose(direct, projected)
    no_cross = within.copy()
    no_cross[:136, 136:] = 0
    no_cross[136:, :136] = 0
    assert np.allclose(transform(no_cross, "joint"), transform(no_cross, "block"), atol=1e-10)
    only_cov = transform(within, "coverage")
    assert not np.any(only_cov[:136])
    assert np.allclose(only_cov[136:].T @ within[136:, 136:] @ only_cov[136:], np.eye(4))
    midpoint = 0.25*within + np.eye(dim)
    signal = transform(within, "signal_joint", midpoint)
    expected_metric = np.linalg.inv(within) - np.linalg.inv(2*midpoint + 0.5*within)
    assert np.allclose(signal @ signal.T, expected_metric, atol=1e-10)
    assert np.allclose(transform(no_cross, "signal_joint", .25*no_cross + np.eye(dim)) @
                       transform(no_cross, "signal_joint", .25*no_cross + np.eye(dim)).T,
                       transform(no_cross, "signal_block", .25*no_cross + np.eye(dim)) @
                       transform(no_cross, "signal_block", .25*no_cross + np.eye(dim)).T, atol=1e-10)
    latent = rng.normal(size=(300, dim))
    split = latent[:, None, :] + rng.normal(size=(300, 2, dim))*0.1
    metrics = cross_validate(split, [f"parent{i}" for i in range(300)], ("joint", "block", "coverage", "signal_joint", "signal_block", "signal_coverage"))
    assert len(metrics) == 14
    assert all(m["same_closer_fraction"] > 0.95 for m in metrics)
    print("PASS metric whitening, Mahalanobis equivalence, block ablation, coverage control, held-parent proxy")


if __name__ == "__main__":
    with threadpool_limits(limits=1): main()
