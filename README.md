# Bayesian-Relevance-Networks-in-Matlab

This Matlab source code is in support of the publication "Uncovering Robust Patterns of MicroRNA Co-Expression across Cancers using Bayesian Relevance Networks" by Ramachandran, Sanchez-Taltavull, and Perkins, PLoS ONE, Vol. 12, No. 8, Art. e0183103 (also appeared at GLBIO 2017, where it won "Outstanding Presentation" prize).

See the tops of the .m files for descriptions of functionality. The file TestScript.m demonstrates how the code would typically be used to perform an analysis of a matrix of count data. The basic steps are to 1) Compute Bayesian (or Pearson) correlations between entities (genes, microRNAs, etc.), 2) Compute a null distribution of those correlations by permuting the data, and 3) threshold the correlations based on the permuted distribution and a desired false discovery rate threshold.
