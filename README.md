## ⚠️ Aviso Importante

Este repositório foi desenvolvido **exclusivamente para fins de estudo, aprendizado e prática** em bioinformática e ciência de dados. Os métodos, análises e resultados apresentados **não devem ser utilizados para aplicações clínicas, diagnósticas ou decisões médicas**.

O projeto tem caráter **exploratório e educacional**, servindo como base para experimentação com técnicas de aprendizado de máquina aplicadas a sequências biológicas. As interpretações biológicas não substituem análises filogenéticas formais ou validações experimentais.

Utilize este material como referência acadêmica e para fins de aprendizado.

# 🧬 Clustering Hierárquico de Sequências Genômicas do SARS-CoV-2

##  Visão Geral
Este projeto realiza a análise de similaridade e o agrupamento hierárquico de sequências genômicas do SARS-CoV-2 utilizando técnicas de aprendizado de máquina aplicadas a dados biológicos. A abordagem combina a representação das sequências por k-mers com vetorização TF-IDF e o uso da distância cosseno para identificar padrões de similaridade genética.

Os resultados são explorados por meio de visualizações gráficas, como mapas de calor e dendrogramas, além de uma avaliação quantitativa da qualidade dos agrupamentos utilizando o Silhouette Score.

---

## 📂 Dataset
- **Fonte:** Kaggle  
- **Nome:** *Genetic Sequences for the SARS-CoV-2 Coronavirus*  
- **Formato original:** CSV  
- **Formato de trabalho:** FASTA (gerado durante o pré-processamento)

Cada sequência representa um fragmento genômico do SARS-CoV-2 e é tratada como uma sequência de caracteres para fins de análise computacional.

---

##  Metodologia

### 1. Pré-processamento dos Dados
- Conversão das sequências do formato CSV para FASTA
- Validação de identificadores e integridade das sequências
- Leitura e manipulação das sequências utilizando Biopython

### 2. Extração de Características
- Representação das sequências por meio de **k-mers** de tamanho fixo
- Vetorização utilizando **TF-IDF** com n-gramas de caracteres
- Normalização dos vetores para padronização das magnitudes

### 3. Cálculo de Distâncias
- Cálculo da dissimilaridade entre sequências utilizando **distância cosseno**
- Construção de uma matriz de distâncias para análise exploratória

### 4. Clustering Hierárquico
- Agrupamento hierárquico aglomerativo utilizando o método **UPGMA (average linkage)**
- Visualização da estrutura hierárquica por meio de **dendrogramas**

### 5. Avaliação dos Clusters
- Avaliação da qualidade dos agrupamentos utilizando o **Silhouette Score**
- Comparação entre diferentes números de clusters para auxiliar na escolha do particionamento ideal

---

## 📊 Visualizações Geradas
O projeto gera automaticamente as seguintes visualizações:
- **Mapa de calor (heatmap)** da matriz de distâncias
- **Dendrograma** do clustering hierárquico
- **Gráfico do Silhouette Score** em função do número de clusters

As figuras são salvas localmente, permitindo reprodutibilidade e uso em relatórios ou apresentações.

---

## 🛠️ Tecnologias Utilizadas
- **Biopython** – manipulação de sequências biológicas  
- **Scikit-learn** – vetorização TF-IDF e métricas de avaliação  
- **SciPy** – cálculo de distâncias e clustering hierárquico  
- **NumPy** – operações numéricas  
- **Matplotlib e Seaborn** – visualização de dados  

---

## 🎯 Aplicações
- Análise exploratória de variabilidade genética viral  
- Identificação de padrões de similaridade entre genomas  
- Base metodológica para estudos em bioinformática e filogenia computacional  
- Projeto demonstrativo de aplicação de Machine Learning em Biologia  
