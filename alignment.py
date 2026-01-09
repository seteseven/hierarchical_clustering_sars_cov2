"""
Clustering Hierárquico de Sequências do SARS-CoV-2

Este projeto realiza a análise de similaridade entre sequências genômicas do
SARS-CoV-2 utilizando vetorização por k-mers com TF-IDF e distância cosseno.
O agrupamento é feito por clustering hierárquico (UPGMA), com visualização por
heatmap e dendrograma, e avaliação do número de clusters via Silhouette Score.

Bibliotecas: Biopython, Scikit-learn, SciPy, NumPy, Matplotlib e Seaborn.
Dataset: Genetic Sequences for the SARS-CoV-2 Coronavirus (Kaggle).
"""
# Link do dataset usado: https://www.kaggle.com/datasets/priteshshrivastava/genetic-sequences-for-the-sarscov2-coronavirus
# =========================
# Imports
# =========================

import csv
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from Bio import SeqIO
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics import silhouette_score

from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster



# Converte o CSV contendo sequências genéticas em um arquivo FASTA.
def csv_fasta(file:str,seq_col:int,name_col:int,file_name_output:str):
    with open(file, 'r') as f:
        if f is not None:
            leitor = csv.reader(f)
            linhas = []
            next(leitor, None)
            for linha in leitor:
                if linha is not None:
                    linhas.append(linha)
        with open(f'{file_name_output}.fasta', 'w') as seq:
            for linha in linhas:
                seq.write('>' + linha[name_col] + '\n' + linha[seq_col] + '\n')





# carrega o arquivo FASTA e retorna as sequências e os IDs
def load_fasta(fasta_file):

    sequences = []
    ids = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq).upper())
        ids.append(record.id)

    return sequences, ids

# Converte as sequências em um vetor de k-mers
def tfidf_kmers(sequences, k=6):

    vectorizer = TfidfVectorizer(
        analyzer="char",
        ngram_range=(k, k),
        lowercase=False,
        norm="l2"
    )

    X = vectorizer.fit_transform(sequences)
    return X

# clustering hierárquico
def hierarchical_clustering(X, method="average", metric="cosine"):
    """
    Calcula a matriz de distância e executa clustering hierárquico.
    """
    dist_vector = pdist(X.toarray(), metric=metric)
    dist_matrix = squareform(dist_vector)
    Z = linkage(dist_vector, method=method)

    return dist_vector, dist_matrix, Z


# Plot dos clusters
def plot_heatmap(dist_matrix, output="heatmap.png"):
    plt.figure(figsize=(8, 6))
    sns.heatmap(dist_matrix, cmap="viridis", xticklabels=False, yticklabels=False)
    plt.title("Distance Matrix Heatmap")
    plt.tight_layout()
    plt.savefig(output)
    plt.close()


# Plot do dendrograma
def plot_dendrogram(Z, labels, output="dendrogram.png", constrained_layout=True):
    plt.figure(figsize=(16, 8))
    dendrogram(Z, labels=labels, leaf_rotation=90, leaf_font_size=4)
    plt.ylabel("Cosine distance")
    plt.title("Hierarchical Clustering (UPGMA)")
    plt.savefig(output)
    plt.close()


# Avalisação de clusters
def silhouette_analysis(X, Z, k_min=2, k_max=10, metric="cosine"):
    scores = []
    ks = range(k_min, k_max)

    for k in ks:
        clusters = fcluster(Z, k, criterion="maxclust")
        score = silhouette_score(X, clusters, metric=metric)
        scores.append(score)

    return ks, scores

# plot dos scores
def plot_silhouette(ks, scores, output="silhouette_score.png"):
    plt.figure(figsize=(6, 4))
    plt.plot(ks, scores, marker="o")
    plt.xlabel("Number of clusters")
    plt.ylabel("Silhouette score")
    plt.title("Cluster Validation")
    plt.tight_layout()
    plt.savefig(output)
    plt.close()


if __name__ == "__main__":

    # 1. Converter CSV → FASTA
    csv_fasta(file='Genetic-Sequences-for-the-SARS-CoV-2-Coronavirus.csv',seq_col=6,name_col=0,file_name_output='sequencias')

    # 2. Ler FASTA
    sequences, ids = load_fasta("sequencias.fasta")

    # 3. Vetorização
    X = tfidf_kmers(sequences, k=6)

    # 4. Clustering
    dist_vector, dist_matrix, Z = hierarchical_clustering(X)

    # 5. Visualizações
    plot_heatmap(dist_matrix)
    plot_dendrogram(Z, ids)

    # 6. Validação
    ks, scores = silhouette_analysis(X, Z)
    plot_silhouette(ks, scores)
