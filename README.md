## Reproducibility


This GitHub repository contains the Python code for reproducing the analysis and figures in the following manuscript:

Van Hoan Do, Stefan Canzar, A generalization of tSNE and UMAP to single cell multimodal
data, xxx.


### Installing of Jvis

The instruction of how to install the package:
https://pypi.org/project/Jvis-learn/

### Reproducibility

1/ Section "Proofs of principle". Here we analyzed a SNARE-seq data set and performed a number of experiments to validate the ability of Jvis to learn the correct representation under the noise. The data and python scripts to reproduce this section are provided in the directory "proof_of_principle"

2/ Section "Jvis utilizes multi-modal data to resolve subtle transcriptomic difference". We ran Jvis on two CITE-seq datasets and compared it to the visualizations based on each modality alone. The data and python scripts to reproduce this section are provided in the directory "cell_heterogeneity"

3/ Section "Jvis improves the visualization of joint velocity landscapes of protein and RNA". In this section we visualized the single cell acceleration on the joint visualization in constrast to the transcriptomic space as in the original publication [1]. The datasets are download from the paper:
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1945-3#availability-of-data-and-materials
* Gorin G, Svensson V, Pachter L. CITE-seq protein and mRNA counts. figshare; 2019 [cited 2020 Jan 20]. Available from: https://figshare.com/articles/CITE-seq_protein_and_mRNA_counts/8309696
* Gorin G, Svensson V, Pachter L. REAP-seq protein and mRNA counts. figshare; 2019 [cited 2020 Jan 20]. Available from: https://figshare.com/articles/REAP-seq_protein_and_mRNA_counts/8309708
* Gorin G, Svensson V, Pachter L. ECCITE-seq protein and mRNA counts. figshare; 2019 [cited 2020 Jan 20]. Available from: https://figshare.com/articles/ECCITE-seq_protein_and_mRNA_counts/8309714
* Gorin G, Pachter L, Svensson V. 10X protein and mRNA counts. figshare; 2019 [cited 2020 Jan 20]. Available from: https://figshare.com/articles/10X_protein_and_mRNA_counts/9912734/1

The Python scripts for reproducing this section are given in the directory "sc_velocity"

#### Reference
[1] Gorin, G., Svensson, V. & Pachter, L. Protein velocity and acceleration from single-cell multiomics experiments. Genome Biol 21, 39 (2020). https://doi.org/10.1186/s13059-020-1945-3


