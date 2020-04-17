In module detecting,we need input preprocessed,normalized,filtered TCGA HCC RNA-seq, its normal and tumor label,
DEGs, Regnetwork, then by implementing our phenotype module detection algorithm, we can get lots of modules.
Then, we deleted repeated modules
In module ranking,we input detected and unique modules and weighted Regnetwork, we implement our block-based module ranking algorithm.
Then,we rank modules by their global PageRank values.
We select the top-10 ranked modules for fruther declaration.

