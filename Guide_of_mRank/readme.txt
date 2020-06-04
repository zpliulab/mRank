(1) Module detection 
In module detecting, we firstly input preprocessed (including filtering and normalizing )TCGA HCC RNA-seq, its normal and tumor label,
DEGs(differetial expressed genes), RegNetwork. Then by implementing our phenotype module detection algorithm, we can get lots of modules.
At last, we need to delete those repeated modules if the same gene sets exist in more than one module .
(2) Module ranking
In module ranking, we firstly input detected unique modules and weighted RegNetwork, we implement our block-based module ranking algorithm. Then, we rank modules by their global PageRank values.
Finally, We select the top-10 ranked modules for fruther declaration.

