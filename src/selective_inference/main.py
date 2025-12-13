import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering, KMeans
from statsmodels.stats.weightstats import ztest
from pathlib import Path

images_path = Path("images/")
images_path.mkdir(parents=True, exist_ok=True)


M = 2000
n = 100
p = 5


def generate_data():
    M = np.zeros((n, p))
    mu = M.ravel(order="F")
    i = np.arange(n * p) + 1
    j = np.arange(n * p) + 1
    C = 0.98 ** (np.abs(i[:, None] - j[None, :]))
    X = np.random.multivariate_normal(
        mean=mu,
        cov=C,
    )
    return X


def pvalues_cah():
    pvalues = np.empty(M)
    for i in range(M):
        if i % 100 == 0:
            print(f"Simulation {i}")
        X = generate_data()
        X = X.reshape((n, p), order="F")
        hac = AgglomerativeClustering(linkage="average")
        labels = hac.fit_predict(X)
        unique_labels = np.unique(labels)
        if unique_labels.size < 2:
            raise ValueError("La CAH a défini moins de 2 clusters.")
        clusters = np.random.choice(unique_labels.size, 2, replace=False)
        print(f"Nombre de clusters {unique_labels.size}")
        print(f"Clusters choisis: {clusters[0]}-{clusters[1]}")

        c1, c2 = unique_labels[clusters[0]], unique_labels[clusters[1]]
        g1 = X[labels == c1].ravel()
        g2 = X[labels == c2].ravel()

        stat_overall, pvalue = ztest(g1, g2, alternative="two-sided")
        pvalues[i] = pvalue

    y = np.arange(1, len(pvalues) + 1) / len(pvalues)
    plt.plot(np.sort(pvalues), y)
    plt.plot(y, y, "--", color="gray")
    plt.title("CAH - average linkage", loc="left")
    plt.xlabel("p-value")
    plt.grid()
    plt.savefig(images_path / "cah_pvalues.png", format="png")


def pvalues_kmeans(n_clusters=3):
    pvalues = np.empty(M)
    for i in range(M):
        if i % 100 == 0:
            print(f"Simulation {i}")
        X = generate_data()
        X = X.reshape((n, p), order="F")
        model = KMeans(n_clusters=n_clusters, n_init="auto").fit(X)
        labels = model.labels_
        unique_labels = np.unique(labels)
        clusters = np.random.choice(unique_labels.size, 2, replace=False)
        print(f"Nombre de clusters {unique_labels.size}")
        print(f"Clusters choisis: {clusters[0]}-{clusters[1]}")

        c1, c2 = unique_labels[clusters[0]], unique_labels[clusters[1]]
        g1 = X[labels == c1].ravel()
        g2 = X[labels == c2].ravel()

        stat_overall, pvalue = ztest(g1, g2, alternative="two-sided")
        pvalues[i] = pvalue

    y = np.arange(1, len(pvalues) + 1) / len(pvalues)
    plt.plot(np.sort(pvalues), y)
    plt.plot(y, y, "--", color="gray")
    plt.title(f"kmeans - n_clusters = {n_clusters}", loc="left")
    plt.xlabel("p-value")
    plt.grid()
    plt.savefig(images_path / "kmeans_df.png", format="png")


# pvalues_kmeans(n_clusters=3)
pvalues_kmeans()
# TODO: faire aussi avec la regression Lasso
