function [protein_similarity, disease_similarity] = integration_protein_disease_similarity(protein_gauss, disease_gauss)

persistent ss vs

data_dir = 'C:\Users\zxt\Desktop\run_MDHGIBNNR\demo';

if isempty(ss)
    ss = readmatrix(fullfile(data_dir, 'Protein_Sequence_Similarity_Matrix.xlsx'));
end

if isempty(vs)
    vs = readmatrix(fullfile(data_dir, 'Disease_Semantic_Similarity_Matrix.xlsx'));
end

protein_similarity = ss + protein_gauss;

disease_similarity = vs + disease_gauss;

end