# mRank
# method description
mRank is a computational method, which can discover reliable module
biomarkers by driving from phenotypic states of gene expression data, and it consists of two major
procedures: (a) an iterative supervised module detection guided by phenotypic states from a given
network; (b) a block-based module ranking via network topological centrality.
# code description
(1) download, preprocess (including filtering and normalizing) TCGA HCC RNA-seq data with TCGAbiolinks r package

(2) download RegNetwork and map it with TCGA HCC gene expression by calculating differential mutual information

(3) implement phenotype-driven module detection algorithm to detect modules in the weighted RegNetwork 

(4) create hyper-graph by combining RegNetwork and detected modules, then implement block-based module ranking algorithm

(5) rank those modules by their PR values

(6) select the top-10 ranked modules as our HCC candidate network biomarkers, implement NOA analysis on those candidate 
network biomarkers, extract edges with selected GO terms such as cell death, cell proliferation and so on

(7) download GSE14520, GSE22058, GSE25097, GSE45436, and GSE64041 from GEO database, map them with their corresponding GPL data
and preprocess them then use them for validation

(8) use igraph package to split RegNetwork into small subnetworks

(9)  input small subnetworks to GSEA, GSA, SAFE( r codes), Crank and Conductance(http://snap.stanford.edu/crank/,  C++, .exe)
