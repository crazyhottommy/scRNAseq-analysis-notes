# scRNAseq-analysis-notes
my scRNAseq analysis notes 

## The reason
Single cell RNAseq is becoming more and more popular, and as a technique, it might become as common as PCR. I just got some 10x genomics single cell RNAseq data to play with, it is a good time for me to take down notes here. I hope it is useful for other people as well.


## readings before doing anything

### single cell tutorials
* [Course material in notebook format for learning about single cell bioinformatics methods](https://github.com/YeoLab/single-cell-bioinformatics)
* [Analysis of single cell RNA-seq data course, Cambridge University](https://github.com/hemberg-lab/scRNA.seq.course) Great tutorial!
* [f1000 workflow paper A step-by-step workflow for low-level analysis of single-cell RNA-seq data](http://f1000research.com/articles/5-2122/v1) by Aaron Lun, the athour of diffHiC, GenomicInteractions and csaw.
* [2016 Bioconductor workshop: Analysis of single-cell RNA-seq data with R and Bioconductor](https://github.com/drisso/bioc2016singlecell)
* [paper: Single-Cell Transcriptomics Bioinformatics and Computational Challenges](http://journal.frontiersin.org/article/10.3389/fgene.2016.00163/full)


### single cell RNA-seq normalization
* [paper: Assessment of single cell RNA-seq normalization methods](http://biorxiv.org/content/early/2016/07/17/064329)
* [paper: A practical guide to single-cell RNA-sequencing for biomedical research and clinical applications](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0467-4)
* [Normalizing single-cell RNA sequencing data: challenges and opportunities](https://www.nature.com/nmeth/journal/v14/n6/full/nmeth.4292.html) Nature Methods
* [SinQC: A Method and Tool to Control Single-cell RNA-seq Data Quality](http://www.morgridge.net/SinQC.html).
* [Scone](https://github.com/YosefLab/scone) Single-Cell Overview of Normalized Expression data

### single cell batch effect
* [Overcoming confounding plate effects in differential expression analyses of single-cell RNA-seq data](http://biorxiv.org/content/early/2016/09/08/073973)

* [Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors](https://www.nature.com/articles/nbt.4091)

### Single cell RNA-seq

* [a collection of single RNA-seq tools by Sean Davis ](https://github.com/seandavi/awesome-single-cell)
* [paper: Design and computational analysis of single-cell RNA-sequencing experiments](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0927-y)
* [paper by Mark Robinson: Bias, Robustness And Scalability In Differential Expression Analysis Of Single-Cell RNA-Seq Data](http://biorxiv.org/content/early/2017/05/28/143289)

> Considerable differences are found between the methods in terms of the number and characteristics of the genes that are called differentially expressed. Pre-filtering of lowly expressed genes can have important effects on the results, particularly for some of the methods originally developed for analysis of bulk RNA-seq data. Generally, however, **methods developed for bulk RNA-seq analysis do not perform notably worse than those developed specifically for scRNA-seq.**

* [paper: Power Analysis of Single Cell RNA‐Sequencing Experiments](http://biorxiv.org/content/early/2016/09/08/073692)
* [paper: The contribution of cell cycle to heterogeneity in single-cell RNA-seq data](http://www.nature.com/nbt/journal/v34/n6/full/nbt.3498.html)
* [paper: Batch effects and the effective design of single-cell gene expression studies](http://biorxiv.org/content/early/2016/07/08/062919)
* [On the widespread and critical impact of systematic bias and batch effects in single-cell RNA-Seq data](http://biorxiv.org/content/early/2015/08/25/025528) 
* [paper: Comparison of methods to detect differentially expressed genes between single-cell populations](http://www.ncbi.nlm.nih.gov/pubmed/27373736)
* [review: Single-cell genome sequencing: current state of the science](http://www.nature.com/nrg/journal/vaop/ncurrent/full/nrg.2015.16.html)  
* [Ginkgo](http://qb.cshl.edu/ginkgo/?q=/ESjKTTeZIdnoGwEB4WTu) A web tool for analyzing single-cell sequencing data.
* [SingleCellExperiment bioc package](http://bioconductor.org/packages/devel/bioc/html/SingleCellExperiment.html) Defines a S4 class for storing data from single-cell experiments. This includes specialized methods to store and retrieve spike-in information, dimensionality reduction coordinates and size factors for each cell, along with the usual metadata for genes and libraries.
* [ASAP](http://biorxiv.org/content/early/2016/12/22/096222): a Web-based platform for the analysis and inter-active visualization of single-cell RNA-seq data
* [Seurat](http://www.satijalab.org/seurat.html) is an R package designed for the analysis and visualization of single cell RNA-seq data. It contains easy-to-use implementations of commonly used analytical techniques, including the identification of highly variable genes, dimensionality reduction (PCA, ICA, t-SNE), standard unsupervised clustering algorithms (density clustering, hierarchical clustering, k-means), and the discovery of differentially expressed genes and markers.
* [R package for the statistical assessment of cell state hierarchies from single-cell RNA-seq data](http://bioconductor.org/packages/devel/bioc/html/sincell.html)  
* [Monocle](http://cole-trapnell-lab.github.io/monocle-release/) Differential expression and time-series analysis for single-cell RNA-Seq and qPCR experiments.
* Single Cell Differential Expression: bioconductor package [scde](http://bioconductor.org/packages/devel/bioc/html/scde.html)
* [Sincera](https://research.cchmc.org/pbge/sincera.html):A Computational Pipeline for Single Cell RNA-Seq Profiling Analysis. Bioconductor package will be available soon. 
* [MAST](https://github.com/RGLab/MAST): a flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data
* [scDD](http://biorxiv.org/content/early/2015/11/13/031021): A statistical approach for identifying differential distributions in single-cell RNA-seq experiments
* [Inference and visualisation of Single-Cell RNA-seq Data data as a hierarchical tree structure: bioconductor CellTree](http://bioconductor.org/packages/devel/bioc/html/cellTree.html)  
* [Fast and accurate single-cell RNA-Seq analysis by clustering of transcript-compatibility counts](http://biorxiv.org/content/early/2016/01/15/036863) by Lior Pachter et.al
* [cellity](https://github.com/ti243/cellity): Classification of low quality cells in scRNA-seq data using R.
* [bioconductor: using scran to perform basic analyses of single-cell RNA-seq data](http://bioconductor.org/packages/devel/bioc/vignettes/scran/inst/doc/scran.html)
* [scater](https://github.com/davismcc/scater): single-cell analysis toolkit for expression with R
* [Monovar](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3835.html): single-nucleotide variant detection in single cells
* [paper: Comparison of methods to detect differentially expressed genes between single-cell populations](http://m.bib.oxfordjournals.org/content/early/2016/07/02/bib.bbw057.full)
* [Single-cell mRNA quantification and differential analysis with Census](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.4150.html)
* [CIDR](https://github.com/VCCRI/CIDR): Ultrafast and accurate clustering through imputation for single-cell RNA-seq data
* [CellView](http://biorxiv.org/content/early/2017/04/04/123810): Interactive Exploration Of High Dimensional Single Cell RNA-Seq Data
* [Scanpy](https://github.com/theislab/scanpy) is a scalable toolkit for analyzing single-cell gene expression data. It includes preprocessing, visualization, clustering, pseudotime and trajectory inference, differential expression testing and simulation of gene regulatory networks. The Python-based implementation efficiently deals with datasets of more than one million cells.

### single cell RNA-seq clustering

* [Geometry of the Gene Expression Space of Individual Cells](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004224)
* [pcaReduce](http://biorxiv.org/content/early/2015/09/08/026385): Hierarchical Clustering of Single Cell Transcriptional Profiles.
* [CountClust](https://www.bioconductor.org/packages/3.3/bioc/html/CountClust.html): Clustering and Visualizing RNA-Seq Expression Data using Grade of Membership Models. Fits grade of membership models (GoM, also known as admixture models) to cluster RNA-seq gene expression count data, identifies characteristic genes driving cluster memberships, and provides a visual summary of the cluster memberships
* [FastProject](http://biorxiv.org/content/early/2016/03/12/043463): A Tool for Low-Dimensional Analysis of Single-Cell RNA-Seq Data
* [SNN-Cliq](http://bioinformatics.oxfordjournals.org/content/31/12/1974.full) Identification of cell types from single-cell transcriptomes using a novel clustering method
* [Compare clusterings for single-cell sequencing](http://bioconductor.org/packages/devel/bioc/html/clusterExperiment.html) bioconductor package.The goal of this package is to encourage the user to try many different clustering algorithms in one package structure. We give tools for running many different clusterings and choices of parameters. We also provide visualization to compare many different clusterings and algorithm tools to find common shared clustering patterns.
* [CIDR: Ultrafast and accurate clustering through imputation for single cell RNA-Seq data](http://biorxiv.org/content/early/2016/08/10/068775)
* [SC3](http://bioconductor.org/packages/release/bioc/html/SC3.html)- consensus clustering of single-cell RNA-Seq data.  SC3 achieves high accuracy and robustness by consistently integrating different clustering solutions through a consensus approach. Tests on twelve published datasets show that SC3 outperforms five existing methods while remaining scalable, as shown by the analysis of a large dataset containing 44,808 cells. Moreover, an interactive graphical implementation makes SC3 accessible to a wide audience of users, and SC3 aids biological interpretation by identifying marker genes, differentially expressed genes and outlier cells.
* [GiniClust2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1431-3): a cluster-aware, weighted ensemble clustering method for cell-type detection
* [FateID infers cell fate bias in multipotent progenitors from single-cell RNA-seq data](https://www.nature.com/articles/nmeth.4662)
* [matchSCore: Matching Single-Cell Phenotypes Across Tools and Experiments](https://www.biorxiv.org/content/early/2018/05/07/314831) In this work we introduce matchSCore (https://github.com/elimereu/matchSCore), an approach to match cell populations fast across tools, experiments and technologies. We compared 14 computational methods and evaluated their accuracy in clustering and gene marker identification in simulated data sets.
* [Cluster Headache: Comparing Clustering Tools for 10X Single Cell Sequencing Data](https://www.biorxiv.org/content/early/2017/10/19/203752)
* The [celaref](https://github.com/MonashBioinformaticsPlatform/celaref) (cell labelling by reference) package aims to streamline the cell-type identification step, by suggesting cluster labels on the basis of similarity to an already-characterised reference dataset - wheather that's from a similar experiment performed previously in the same lab, or from a public dataset from a similar sample.

### dimention reduction and visualization of clusters

* [Principal Component Analysis Explained Visually](http://setosa.io/ev/principal-component-analysis/)  
* [PCA, MDS, k-means, Hierarchical clustering and heatmap](https://rpubs.com/crazyhottommy/PCA_MDS). I wrote it.
* [A tale of two heatmaps](https://rpubs.com/crazyhottommy/a-tale-of-two-heatmap-functions). I wrote it.
* [Heatmap demystified](https://rpubs.com/crazyhottommy/heatmap_demystified). I wrote it.
* [Cluster Analysis in R - Unsupervised machine learning](http://www.sthda.com/english/wiki/cluster-analysis-in-r-unsupervised-machine-learning#at_pco=smlre-1.0&at_si=58765a95fcb21379&at_ab=per-2&at_pos=3&at_tot=4) very practical intro on STHDA website.
* [I wrote on PCA, and heatmaps on Rpub](https://rpubs.com/crazyhottommy)
* A most read for clustering analysis for high-dimentional biological data:[Avoiding common pitfalls when clustering
biological data](http://stke.sciencemag.org/content/9/432/re6)
* [How does gene expression clustering work?](http://www.nature.com/nbt/journal/v23/n12/full/nbt1205-1499.html) A must read for 
clustering.
* [How to read PCA plots for scRNAseq](http://www.nxn.se/valent/2017/6/12/how-to-read-pca-plots) by VALENTINE SVENSSON.

<blockquote class="twitter-tweet" data-lang="en"><p lang="en" dir="ltr">See <a href="https://t.co/yxCb85ctL1">https://t.co/yxCb85ctL1</a>: &quot;MDS best choice for preserving outliers,  PCA for variance, &amp; T-SNE for clusters&quot; <a href="https://twitter.com/mikelove">@mikelove</a> <a href="https://twitter.com/AndrewLBeam">@AndrewLBeam</a></p>&mdash; Rileen Sinha (@RileenSinha) <a href="https://twitter.com/RileenSinha/status/768873620521250816">August 25, 2016</a></blockquote>
<script async src="//platform.twitter.com/widgets.js" charset="utf-8"></script>

[paper: Outlier Preservation by Dimensionality Reduction Techniques](http://oai.cwi.nl/oai/asset/22628/22628B.pdf)
>"MDS best choice for preserving outliers, PCA for variance, & T-SNE for clusters"

* [How to Use t-SNE Effectively](http://distill.pub/2016/misread-tsne/)
* [t-SNE explained in plain javascript](https://beta.observablehq.com/@nstrayer/t-sne-explained-in-plain-javascript)
* [Rtsne](https://github.com/jkrijthe/Rtsne) R package for T-SNE
* [rtsne](https://github.com/jdonaldson/rtsne) An R package for t-SNE (t-Distributed Stochastic Neighbor Embedding)
a bug was in `rtsne`: https://gist.github.com/mikelove/74bbf5c41010ae1dc94281cface90d32
* [t-SNE-Heatmaps](https://github.com/KlugerLab/t-SNE-Heatmaps) Beta version of 1D t-SNE heatmaps to visualize expression patterns of hundreds of genes simultaneously in scRNA-seq.
* [PHATE dimensionality reduction method](https://github.com/KrishnaswamyLab/PHATE) paper: http://biorxiv.org/content/early/2017/03/24/120378  PHATE also uncovers and emphasizes progression and transitions (when they exist) in the data, which are often missed in other visualization-capable methods. Such patterns are especially important in biological data that contain, for example, single-cell phenotypes at different phases of differentiation, patients at different stages of disease progression, and gut microbial compositions that vary gradually between individuals, even of the same enterotype.

* [Uniform Manifold Approximation and Projection (UMAP)](https://github.com/lmcinnes/umap) is a dimension reduction technique that can be used for visualisation similarly to t-SNE, but also for general non-linear dimension reduction. The algorithm is founded on three assumptions about the data. Run from R: https://gist.github.com/crazyhottommy/caa5a4a4b07ee7f08f7d0649780832ef
* [umapr](https://github.com/ropenscilabs/umapr) UMAP dimensionality reduction in R

### interesting papers to read 

* [Comparative analysis of droplet-based ultra-high-throughput single-cell RNA-seq systems](https://www.biorxiv.org/content/early/2018/05/02/313130)

* [A Single-Cell Transcriptome Atlas of the Aging Drosophila Brain](https://www.cell.com/cell/fulltext/S0092-8674(18)30720-7)
* [Bias, robustness and scalability in single-cell differential expression analysis](https://www.nature.com/articles/nmeth.4612) by Mark Robinson
 
### database
* [single cell expression atlas](https://www.ebi.ac.uk/gxa/sc/home) 

### advance of scRNA-seq tech
* [Single-cell profiling of the developing mouse brain and spinal cord with split-pool barcoding](http://science.sciencemag.org/content/360/6385/176) no isolation of single cells needed!

* [Dynamics and Spatial Genomics of the Nascent Transcriptome by Intron seqFISH](https://www.cell.com/cell/fulltext/S0092-8674(18)30647-0)

* [Highly Multiplexed Single-Cell RNA-seq for Defining Cell Population and Transcriptional Spaces](https://www.biorxiv.org/content/early/2018/05/05/315333) blog post by Lior Patcher [The benefits of multiplexing](https://liorpachter.wordpress.com/2018/05/30/the-benefits-of-multiplexing/). Need to re-read carefully.

* [Three-dimensional intact-tissue sequencing of single-cell transcriptional states](http://science.sciencemag.org/content/early/2018/06/20/science.aat5691)


### pseudotemporal modelling

* [Uncovering pseudotemporal trajectories with covariates from single cell and bulk expression data](https://www.nature.com/articles/s41467-018-04696-6)

* all different algorithms https://github.com/agitter/single-cell-pseudotime
* [Genomic trajectories with heterogeneous genetic and environmental backgrounds](https://bioconductor.org/packages/release/bioc/html/phenopath.html)

* [A descriptive marker gene approach to single-cell pseudotime inference](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty498/5043298)

* [A collection of 50 trajectory inference methods within a common interface](https://github.com/dynverse/dynmethods) take a look of this!

### large scale single cell analysis
* [bigSCale](https://genome.cshlp.org/content/early/2018/05/03/gr.230771.117.abstract): an analytical framework for big-scale single-cell data. [github link](https://github.com/iaconogi/bigSCale) for millions of cells (starts with a count matrix)

* [Alevin: An integrated method for dscRNA-seq quantification](https://www.biorxiv.org/content/early/2018/06/01/335000) based on Salmon.

## The field is advancing so fast!!

check this website for the tools being added:  
https://www.scrna-tools.org/

paper published:  
[Exploring the single-cell RNA-seq analysis landscape with the scRNA-tools database](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006245)


### contamination of 10x data
https://twitter.com/constantamateur/status/994832241107849216?s=11

> Did you know that droplet based single cell RNA-seq data (like 10X) is contaminated by ambient mRNAs?  Good news though, we've written a paper (https://www.biorxiv.org/content/early/2018/04/20/303727 …) and created an R package called SoupX (https://github.com/constantAmateur/SoupX) to fix this problem.

>Is this really a problem?  It depends on your experiment.  Contamination ranges from 2% - 50%.  10% seems common; it's 8% for 10X PBMC data.  Solid tissues are typically worse, but there's no way to know in advance.  Wouldn't you like to know how contaminated your data are?

>These mRNAs come from the single cell suspension fed into the droplet creation system.  They mostly get their from lysed cells and so resemble the cells being studied.  This means the profile of the contamination is experiment specific and creates a batch effect.

[cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) is the toolkit developed by the 10x genomics company to deal with the data.


