# BAMPDA: A matrix refactoring and heterogeneous inference framework for predicting PTM-disease associations.

BAMPDA: A matrix refactoring and heterogeneous inference framework for predicting PTM-disease associations.
Prepared by Necla Nisa Soylu and Emre Sefer
 

## 📖 Overview
Predicting PTM-disease associations is crucial for understanding pathogenic mechanisms. BAMPDA integrates **Bounded Nuclear Norm Regularization (BNNR)** and **Matrix Decomposition and Heterogeneous Graph Inference (MDHGI)** into a unified computational pipeline. By combining disease semantic similarities, protein sequence similarities, and Gaussian interaction profile similarities, BAMPDA effectively captures the complex latent topological structures to infer potential PTM-disease pairs.

---

### Dataset 1
* **Source:** [PhosphoSitePlus®](https://www.phosphosite.org/)
* **Scale:** 1,751 known PTM-disease associations.
* **Entities:** 1,036 proteins and 391 diseases.

### Dataset 2
* **Source:** [PTMD v1.0](http://ptmd.biocuckoo.org/)
* **Preprocessing:** We strictly filtered the original PTMD database to exclude any protein-disease associations that did not involve specific PTM sites, ensuring a highly confident benchmark.
* **Scale:** 905 refined PTM-disease associations.
* **Entities:** 749 proteins and 275 diseases.

---

## 🧮 Data Formulation (Adjacency Matrix)

Although our datasets inherently contain specific PTM classifications, **BAMPDA formulates the prediction task purely as a binary matrix refactoring problem**. Crucially, the specific types of PTMs are *not* provided as explicit input labels to the model during training.

All known protein-disease associations from the datasets are computationally structured into an adjacency matrix $A$ (size $np \times nd$, where $np$ is the number of proteins and $nd$ is the number of diseases):

$$
A(i, j) = 
\begin{cases} 
1, & \text{if protein } p_{i} \text{ is associated with disease } d_{j} \\ 
0, & \text{otherwise} 
\end{cases}
$$

*Note: This purely binary formulation ensures that the diverse PTM states recovered in our subsequent case studies are genuine computational discoveries of latent biological mechanisms, rather than mere memorization of initial dataset labels.*

Each dataset folder (`Dataset1` and `Dataset2`) contains the following preprocessed components in Excel format:

**1. Biological Features:**
- `Disease_Name.xlsx`: The list of targeted diseases.
- `Protein_Numbers.xlsx`: The list of PTM-related proteins.
- `Disease_Semantic_Similarity_Matrix.xlsx`: Calculated based on the MeSH descriptors of diseases.
- `Protein_Sequence_Similarity_Matrix.xlsx`: Calculated using standard sequence alignment tools.

**2. Gaussian Interaction Profile (GIP) Similarities:**
- `Disease_Gaussian_Similarity_Matrix.xlsx`
- `Protein_Gaussian_Similarity_Matrix.xlsx`

**3. Association Matrices:**
- `Protein_Disease_adj.xlsx`: The adjacency matrix representing the topological network.
- `Protein_Disease_Associations.xlsx`: Ground-truth known PTM-disease associations used for training/validation.

---

## ⚙️ Model Architecture & Code Structure
*(借鉴点：将代码文件与处理步骤对应起来)*
The framework is highly modularized. You can find the data processing codes and the core algorithms in the following files:

* **Step 1: Feature Weighting**
  - `protein_disease_weight_matrix.m`: Processes the `Disease_Semantic_Similarity_Matrix` and `Protein_Sequence_Similarity_Matrix` to generate the initial weight matrix.
* **Step 2: Similarity Integration**
  - `integration_protein_disease_similarity.m`: Incorporates the Gaussian similarities (`Disease_Gaussian_Similarity_Matrix.xlsx` and `Protein_Gaussian_Similarity_Matrix.xlsx`) to compute the final, integrated `Disease_Similarity` and `Protein_Similarity` matrices.
* **Step 3: Core Prediction (BAMPDA Integration)**
  - `BNNR.m`: Executes the Bounded Nuclear Norm Regularization algorithm.
  - `MDHGI.m`: Executes the Matrix Decomposition and Heterogeneous Graph Inference algorithm.
  - `main.m`: The primary orchestration script. It loads the datasets, calls the preprocessing modules, executes the BAMPDA core, and outputs the final prediction scores and performance metrics (e.g., AUC).

---

## 🚀 Step-by-Step Illustrated Protocol (Quick Start)
*(针对审稿人意见的核心回应：图文教程)*

### 1. Requirements
- **Environment:** MATLAB R2023b or later.

### 2. Loading the Data
Clone this repository and add the folder to your MATLAB path. Open `main.m`. By default, the script is configured to load `Dataset1`. 

*(Note for authors: Paste a screenshot here showing the MATLAB workspace after loading the Excel files)*

### 3. Running the Framework
Simply click the **"Run"** button in MATLAB or type `main` in the Command Window. The script will automatically trigger the integration modules and start the matrix refactoring process.

*(Note for authors: Paste a screenshot here showing the command window printing out the progress or parameters)*

### 4. Interpreting the Results
Upon completion, the script will output the predicted association score matrix. Higher scores represent a higher probability of association between the corresponding PTM and disease. 

*(Note for authors: Paste a screenshot here showing the final ROC curve figure or the output variable containing the results)*

---
## ✉️ Contact
For any questions regarding the code or datasets, please refer to the manuscript or open an issue in this repository.