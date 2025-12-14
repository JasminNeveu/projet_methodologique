import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from numpy.typing import NDArray
from typing import List

images_path = Path("images/")
images_path.mkdir(parents=True, exist_ok=True)


def generate_data_lasso(n: int):
    p = 8
    i = np.arange(p)
    j = np.arange(p)
    R = 0.5 ** np.abs(i[:, None] - j[None, :])

    X = np.random.multivariate_normal(np.zeros(p), R, size=n)
    beta = np.array([3, 1.5, 0, 0, 2, 0, 0, 0])
    y = X @ beta + np.random.randn(n)

    return y, X


def soft_threshold(u: NDArray, tau: float) -> NDArray:
    return np.sign(u) * np.maximum(np.abs(u) - tau, 0.0)


def proximal_gradient_lasso(
    X: NDArray,
    y: NDArray,
    lam: float,
    w_init: NDArray,
    learning_rate: float,
    tol: float,
    max_iter: int,
) -> NDArray:
    n = X.shape[0]
    w = w_init.copy()

    for _ in range(max_iter):
        w_old = w.copy()
        grad = -(1.0 / n) * X.T @ (y - X @ w)
        w = soft_threshold(w - learning_rate * grad, learning_rate * lam)
        if np.linalg.norm(w - w_old) < tol:
            break

    return w


def lasso_path_prox(
    X: NDArray,
    y: NDArray,
    lambdas_desc: List[float],
    learning_rate: float,
    tol: float,
    max_iter: int,
):
    p = X.shape[1]
    W = np.zeros((p, len(lambdas_desc)))

    w = np.zeros(p)

    for i, lam in enumerate(lambdas_desc):
        w = proximal_gradient_lasso(
            X=X,
            y=y,
            lam=lam,
            w_init=w,
            learning_rate=learning_rate,
            tol=tol,
            max_iter=max_iter,
        )
        W[:, i] = w

    return W


if __name__ == "__main__":
    np.random.seed(0)

    y, X = generate_data_lasso(n=100)

    n = X.shape[0]

    lambda_max = np.linalg.norm(X.T @ y, ord=np.inf) / n

    lambdas_desc = np.logspace(
        np.log10(lambda_max),
        np.log10(lambda_max * 1e-3),
        200,
    )

    W = lasso_path_prox(
        X=X,
        y=y,
        lambdas_desc=lambdas_desc,
        learning_rate=0.05,
        tol=1e-6,
        max_iter=5000,
    )

    lambdas_plot = lambdas_desc[::-1]
    W_plot = W[:, ::-1]

    for j in range(W_plot.shape[0]):
        plt.semilogx(lambdas_plot, W_plot[j])

    plt.xlabel("λ")
    plt.ylabel("coefficients")
    plt.savefig(images_path / "lasso_path.png", format="png")
