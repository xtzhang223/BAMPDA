# BAMPDA: A Matrix Refactoring and Heterogeneous Inference Framework for Predicting PTM-Disease Associations

BAMPDA is a matrix refactoring and heterogeneous inference framework for predicting potential post-translational modification (PTM)-disease associations.

## Overview

 BAMPDA integrates **Bounded Nuclear Norm Regularization (BNNR)** and **Matrix Decomposition and Heterogeneous Graph Inference (MDHGI)** into a unified computational framework.

By combining disease semantic similarity, protein sequence similarity, and Gaussian interaction profile similarity, BAMPDA captures latent topological structures in heterogeneous biological networks and infers potential PTM-disease associations.

## Datasets

### Dataset 1

- **Source:** [PhosphoSitePlus]
- **Scale:** 1,751 known PTM-disease associations
- **Entities:** 1,036 proteins and 391 diseases

### Dataset 2

- **Source:** [PTMD v1.0]
- **Preprocessing:** Protein-disease associations without specific PTM-site information were removed to construct a high-confidence benchmark dataset.
- **Scale:** 905 refined PTM-disease associations
- **Entities:** 749 proteins and 275 diseases

## Data Formulation

Although the original datasets contain PTM-related information, BAMPDA formulates the prediction task as a binary association matrix completion problem. Known protein-disease associations are represented by an adjacency matrix \(A\), where rows correspond to proteins and columns correspond to diseases.

\[
A(i,j)=
\begin{cases}
1, & \text{if protein } p_i \text{ is associated with disease } d_j \\
0, & \text{otherwise}
\end{cases}
\]

This binary formulation allows BAMPDA to infer potential PTM-disease associations based on latent topological and similarity information rather than directly relying on PTM-type labels.

Each dataset folder contains the following files:

### Biological features

- `Disease_Name.xlsx`: Disease name information
- `Protein_Numbers.xlsx`: Protein index information
- `Disease_Semantic_Similarity_Matrix.xlsx`: Disease semantic similarity matrix
- `Protein_Sequence_Similarity_Matrix.xlsx`: Protein sequence similarity matrix

### Gaussian interaction profile similarities

- `Disease_Gaussian_Similarity_Matrix.xlsx`
- `Protein_Gaussian_Similarity_Matrix.xlsx`

### Association matrices

- `Protein_Disease_adj.xlsx`: Known protein-disease association pairs
- `Protein_Disease_adj_name.xlsx`: Known protein-disease association 
- `Protein_Disease_Associations.xlsx`: Binary protein-disease association matrix



### Core MATLAB files

- `main.m`: Main script for running BAMPDA
- `pre_recall.m`: Calculates Precision, Recall, and F1-score based on the predicted association matrix and the ground-truth association matrix.
- `BNNR.m`: Implementation of Bounded Nuclear Norm Regularization
- `MDHGIMDA.m`: Implementation of matrix decomposition and heterogeneous graph inference
- `ensemble_average.m`: Ensemble integration of prediction scores
- `computing_ensemble_ranking.m`: Ranking calculation for cross-validation evaluation
- `integration_protein_disease_similarity.m`: Integration of protein and disease similarity matrices
- `similarity_protein.m`: Calculation of protein Gaussian similarity
- `similarity_disease.m`: Calculation of disease Gaussian similarity

### Supporting MATLAB files

- `logistic.m`
- `lrr.m`
- `softmax_sample.m`
- `solve_l12.m`
- `solve_l2.m`
- `svt.m`
- `t_2.m`
- `test2_generate_array.m`
- `pre_recall.m`



### PTM-stratified robustness analysis

The following Python scripts are used for the PTM-stratified robustness analysis described in the manuscript. Specifically, they split the benchmark associations according to five representative PTM categories and support the evaluation of BAMPDA under uneven PTM coverage.

- `split_Phosphorylation.py`: Extracts phosphorylation-related associations
- `split_Ubiquitination.py`: Extracts ubiquitination-related associations
- `split_Methylation.py`: Extracts methylation-related associations
- `split_Acetylation.py`: Extracts acetylation-related associations
- `split_Glycosylation.py`: Extracts glycosylation-related associations

These scripts were used for the stratified evaluation reported in the robustness experiment, where BAMPDA was evaluated separately on each PTM category using AUC, Recall, Precision, and F1-score.
## Requirements

The code requires the following environment:

- MATLAB R2023b or later
- Python 3.9

