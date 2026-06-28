function [ protein_similarity,disease_similarity ] = integration_protein_disease_similarity(protein_gauss, disease_gauss)

ss = xlsread('C:\Users\zxt\Desktop\run_MDHGIBNNR\demo\Protein_Sequence_Similarity_Matrix.xlsx'); 

vs = xlsread('C:\Users\zxt\Desktop\run_MDHGIBNNR\demo\Disease_Semantic_Similarity_Matrix.xlsx');          

protein_similarity = ss+protein_gauss;     

disease_similarity= vs+disease_gauss;  

