# scRNAseq-analysis-notes
my scRNAseq analysis notes 

## The reason
Single cell RNAseq is becoming more and more popular, and as a technique, it might become as common as PCR. I just got some 10x genomics single cell RNAseq data to play with, it is a good time for me to take down notes here. I hope it is useful for other people as well.


## readings before doing anything

### single cell tutorials
* [Course material in notebook format for learning about single cell bioinformatics methods](https://github.com/YeoLab/single-cell-bioinformatics)
* [Analysis of single cell RNA-seq data course, Cambridge University](https://github.com/hemberg-lab/scRNA.seq.course) Great tutorial!
* [Current best practices in single‐cell RNA‐seq analysis: a tutorial](https://www.embopress.org/doi/10.15252/msb.20188746) github repo https://github.com/theislab/single-cell-tutorial
* [f1000 workflow paper A step-by-step workflow for low-level analysis of single-cell RNA-seq data](http://f1000research.com/articles/5-2122/v1) by Aaron Lun, the athour of diffHiC, GenomicInteractions and csaw.
* [SimplesingleCell](https://bioconductor.org/packages/release/workflows/html/simpleSingleCell.html) A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor by Aaron Lun as well!
* [2016 Bioconductor workshop: Analysis of single-cell RNA-seq data with R and Bioconductor](https://github.com/drisso/bioc2016singlecell)
* [Orchestrating Single-Cell Analysis with Bioconductor](https://bioconductor.github.io/OrchestratingSingleCellAnalysis/)
* [paper: Single-Cell Transcriptomics Bioinformatics and Computational Challenges](http://journal.frontiersin.org/article/10.3389/fgene.2016.00163/full)
* [Single Cell RNA-seq (scRNA-seq) Library Structure](https://teichlab.github.io/scg_lib_structs/)
* [Variance stabilizing scRNA-seq counts](http://www.nxn.se/valent/2017/10/15/variance-stabilizing-scrna-seq-counts) is log2(x+1) reasonable?

### scRNAseq experimental design
* [paper: Design and Analysis of Single-Cell Sequencing Experiments](https://www.cell.com/cell/fulltext/S0092-8674(15)01353-7)
* [paper: Experimental design for single-cell RNA sequencing](https://academic.oup.com/bfg/article/17/4/233/4604806)
* [paper: How to design a single-cell RNA-sequencing experiment: pitfalls, challenges and perspectives](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bby007/4831233)
* [GT-TS: Experimental design for maximizing cell type discovery in single-cell data](https://www.biorxiv.org/content/early/2018/08/07/386540)
* [Tutorial: guidelines for the experimental design of single-cell RNA sequencing studies](https://www.nature.com/articles/s41596-018-0073-y)

### single cell RNA-seq normalization
* [paper: Performance Assessment and Selection of Normalization Procedures for Single-Cell RNA-Seq](https://www.biorxiv.org/content/10.1101/235382v2)
* [paper: Assessment of single cell RNA-seq normalization methods](http://biorxiv.org/content/early/2016/07/17/064329)
* [paper: A practical guide to single-cell RNA-sequencing for biomedical research and clinical applications](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0467-4)
* [Normalizing single-cell RNA sequencing data: challenges and opportunities](https://www.nature.com/nmeth/journal/v14/n6/full/nmeth.4292.html) Nature Methods
* [SinQC: A Method and Tool to Control Single-cell RNA-seq Data Quality](http://www.morgridge.net/SinQC.html).
* [Scone](https://github.com/YosefLab/scone) Single-Cell Overview of Normalized Expression data

### single cell impute 
* [SAVER: gene expression recovery for single-cell RNA sequencing](https://www.nature.com/articles/s41592-018-0033-z)  an expression recovery method for unique molecule index (UMI)-based scRNA-seq data that borrows information across genes and cells to provide accurate expression estimates for all genes.

* [DeepImpute: an accurate, fast and scalable deep neural network method to impute single-cell RNA-Seq data](https://github.com/lanagarmire/DeepImpute) https://www.biorxiv.org/content/early/2018/06/22/353607
* [MAGIC](https://github.com/krishnaswamylab/MAGIC) (Markov Affinity-based Graph Imputation of Cells), is a method for imputing missing values restoring structure of large biological datasets.
* [bayNorm: Bayesian gene expression recovery, imputation and normalisation for single cell RNA-sequencing data](https://www.biorxiv.org/content/early/2018/08/03/384586?) github [page](https://github.com/WT215/bayNorm)
* [Zero-preserving imputation of scRNA-seq data using low-rank approximation](https://www.biorxiv.org/content/early/2018/08/22/397588?rss=1)

### single cell batch effect
* [Overcoming confounding plate effects in differential expression analyses of single-cell RNA-seq data](http://biorxiv.org/content/early/2016/09/08/073973)

* [Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors](https://www.nature.com/articles/nbt.4091)

* [Panoramic stitching of heterogeneous single-cell transcriptomic data](https://www.biorxiv.org/content/early/2018/07/17/371179) Here we present [Scanorama](http://cb.csail.mit.edu/cb/scanorama/), inspired by algorithms for panorama stitching, that overcomes the limitations of existing methods to enable accurate, heterogeneous scRNA-seq data set integration.

* [Fast Batch Alignment of Single Cell Transcriptomes Unifies Multiple Mouse Cell Atlases into an Integrated Landscape](https://www.biorxiv.org/content/early/2018/08/22/397042) [github link](https://github.com/Teichlab/bbknn)

* [Scalable integration of single cell RNAseq data for batch correction and meta analysis](https://github.com/immunogenomics/harmony)

* [liger](https://macoskolab.github.io/liger/)R package for integrating and analyzing multiple single-cell datasets


### Benchmark single cell pipeline

* [MatchScore2](https://github.com/elimereu/matchSCore2) paper: [Benchmarking Single-Cell RNA Sequencing Protocols for Cell Atlas Projects](https://www.biorxiv.org/content/10.1101/630087v1)

* [Benchmarking single cell RNA-sequencing analysis pipelines using mixture control experiments](https://www.nature.com/articles/s41592-019-0425-8) [code and data](https://github.com/LuyiTian/CellBench_data) [cellbench](https://bioconductor.org/packages/release/bioc/html/CellBench.html) bioconductor package.

### Differential expression 

* [A discriminative learning approach to differential expression analysis for single-cell RNA-seq](https://www.nature.com/articles/s41592-018-0303-9) by Lior Patcher group.
* [scde](http://bioconductor.org/packages/release/bioc/html/scde.html) bioconductor package maintained by Jean Fan in Xiaowei Zhuang's lab at Harvard. Need to talk to her once I get a chance.
* [Presto](https://github.com/immunogenomics/presto)Fast Wilcoxon and auROC for single cell RNAseq and scATACseq data. take a look!
* How to compare clusters with multiple samples? https://twitter.com/RoryKirchner/status/1082752967806210048 . work in progess https://github.com/HelenaLC/muscat by Helena from Mark Robinson lab. bioc2019 workshop http://biocworkshops2019.bioconductor.org.s3-website-us-east-1.amazonaws.com/page/muscWorkshop__vignette/
and a blog post by VALENTINE SVENSSON from Lior Patcher's group http://www.nxn.se/valent/2019/2/15/handling-confounded-samples-for-differential-expression-in-scrna-seq-experiments

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

### predict cell type by reference

* [scMatch: a single-cell gene expression profile annotation tool using reference datasets](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz292/5480299)
* [Celaref](https://github.com/MonashBioinformaticsPlatform/celaref)
* [MetaNeighbour](https://github.com/maggiecrow/MetaNeighbor)
* [scID](https://github.com/BatadaLab/scID)
* [scMCA](https://github.com/ggjlab/scMCA)
* [scPred](https://github.com/IMB-Computational-Genomics-Lab/scPred)
* [SingleR](https://github.com/dviraran/SingleR)
* [singleCellNet](https://www.biorxiv.org/content/10.1101/508085v1)

### single cell RNA-seq clustering

* [Single Cell Clustering Comparison](https://jef.works/blog/2018/06/28/single-cell-clustering-comparison/) A blog post.
* [A systematic performance evaluation of clustering methods for single-cell RNA-seq data](https://f1000research.com/articles/7-1141/v1) F1000 paper by Mark Robinson. tl;dr version: "SC3 and Seurat show the most favorable results".
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
* [horseshoe effect from PCA](http://www.huber.embl.de/users/whuber/pub/horseshoe.html) Spurious structures in latent space decomposition and low-dimensional embedding methods
* also read chapter 9 of http://web.stanford.edu/class/bios221/book/Chap-MultivaHetero.html
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
* [A tutorial on t-SNE](http://blog.thegrandlocus.com/2018/08/a-tutorial-on-t-sne-1)
* when to use PCA instead of t-SNE
https://stats.stackexchange.com/questions/238538/are-there-cases-where-pca-is-more-suitable-than-t-sne/249520#249520

>t-SNE is a great piece of Machine Learning but one can find many reasons to use PCA instead of it. Of the top of my head, I will mention five. As most other computational methodologies in use, t-SNE is no silver bullet and there are quite a few reasons that make it a suboptimal choice in some cases. Let me mention some points in brief:

>Stochasticity of final solution. PCA is deterministic; t-SNE is not. One gets a nice visualisation and then her colleague gets another visualisation and then they get artistic which looks better and if a difference of 0.03% in the KL(P||Q) divergence is meaningful... In PCA the correct answer to the question posed is guaranteed. t-SNE might have multiple minima that might lead to different solutions. This necessitates multiple runs as well as raises questions about the reproducibility of the results.

>Interpretability of mapping. This relates to the above point but let's assume that a team has agreed in a particular random seed/run. Now the question becomes what this shows... t-SNE tries to map only local / neighbours correctly so our insights from that embedding should be very cautious; global trends are not accurately represented (and that can be potentially a great thing for visualisation). On the other hand, PCA is just a diagonal rotation of our initial covariance matrix and the eigenvectors represent a new axial system in the space spanned by our original data. We can directly explain what a particular PCA does.

>Application to new/unseen data. t-SNE is not learning a function from the original space to the new (lower) dimensional one and that's a problem. On that matter, t-SNE is a non-parametric learning algorithm so approximating with parametric algorithm is an ill-posed problem. The embedding is learned by directly moving the data across the low dimensional space. That means one does not get an eigenvector or a similar construct to use in new data. In contrast, using PCA the eigenvectors offer a new axes system what can be directly used to project new data. [Apparently one could try training a deep-network to learn the t-SNE mapping (you can hear Dr. van der Maaten at ~46' of this video suggesting something along this lines) but clearly no easy solution exists.]

>Incomplete data. Natively t-SNE does not deal with incomplete data. In fairness, PCA does not deal with them either but numerous extensions of PCA for incomplete data (eg. probabilistic PCA) are out there and are almost standard modelling routines. t-SNE currently cannot handle incomplete data (aside obviously training a probabilistic PCA first and passing the PC scores to t-SNE as inputs).

>The k is not (too) small case. t-SNE solves a problem known as the crowding problem, effectively that somewhat similar points in higher dimension collapsing on top of each other in lower dimensions (more here). Now as you increase the dimensions used the crowding problem gets less severe ie. the problem you are trying to solve through the use of t-SNE gets attenuated. You can work around this issue but it is not trivial. Therefore if you need a k dimensional vector as the reduced set and k is not quite small the optimality of the produce solution is in question. PCA on the other hand offer always the k best linear combination in terms of variance explained. (Thanks to @amoeba for noticing I made a mess when first trying to outline this point.)

>I do not mention issues about computational requirements (eg. speed or memory size) nor issues about selecting relevant hyperparameters (eg. perplexity). I think these are internal issues of the t-SNE methodology and are irrelevant when comparing it to another algorithm.

>To summarise, t-SNE is great but as all algorithms has its limitations when it comes to its applicability. I use t-SNE almost on any new dataset I get my hands on as an explanatory data analysis tool. I think though it has certain limitations that do not make it nearly as applicable as PCA. Let me stress that PCA is not perfect either; for example, the PCA-based visualisations are often inferior to those of t-SNE.

* projection to new data https://twitter.com/EduEyras/status/1032215352623747072
>You can’t add samples to an existing tSNE plot because there is no function outputed by the initial tSNE that maps from the higher dimensional space to the lower dimensions

* [Interpretable dimensionality reduction of single cell transcriptome data with deep generative models](https://www.nature.com/articles/s41467-018-04368-5) 

![](https://github.com/crazyhottommy/scRNAseq-analysis-notes/blob/master/pics/tnse-limit.jpg)
>UMAP is faster, the embeddings are often ++better, and you can use the result to project new data.

* PCA loadings can be used to project new data

e.g. from this paper [Multi-stage Differentiation Defines Melanoma Subtypes with Differential Vulnerability to Drug-Induced Iron-Dependent Oxidative Stress](https://www.cell.com/cancer-cell/abstract/S1535-6108(18)30122-3) Fig 1D.

```{r}
diffStagePCA = prcomp(t(diffStageDataCentered))

# Diff stage PCA (scores for top panel)
diffStagePCA_scores = diffStagePCA$x

# Cell line projected to diff stage PCA (scores for bottom panel)
diffStagePCA_rotation = diffStagePCA$rotation
cellLineProjected_scores <- as.matrix(t(cellLineDataCentered)) %*% as.matrix(diffStagePCA_rotation)
```

* [Sleepwalk: Walk through your embedding](https://anders-biostat.github.io/sleepwalk/) So, can you be sure that the visualisation you get by using t-SNE, UMAP, MDS or the like really give you a faithful representation of your data? Are the points that lie almost on top of each other really all similar? Does the large distance on your 2D representation always mean lots of dissimilarities? Our sleepwalk package for the R statistical programming environment can help you answer these questions.

* [Generalizable and Scalable Visualization of Single-Cell Data Using Neural Networks](https://www.sciencedirect.com/science/article/pii/S2405471218302357?via%3Dihub) standard methods, such as t-stochastic neighbor embedding (t-SNE), are not scalable to datasets with millions of cells and the resulting visualizations cannot be generalized to analyze new datasets. Here we introduce **net-SNE**, a generalizable visualization approach that trains a neural network to learn a mapping function from high-dimensional single-cell gene-expression profiles to a low-dimensional visualization.

* [PHATE dimensionality reduction method](https://github.com/KrishnaswamyLab/PHATE) paper: http://biorxiv.org/content/early/2017/03/24/120378  PHATE also uncovers and emphasizes progression and transitions (when they exist) in the data, which are often missed in other visualization-capable methods. Such patterns are especially important in biological data that contain, for example, single-cell phenotypes at different phases of differentiation, patients at different stages of disease progression, and gut microbial compositions that vary gradually between individuals, even of the same enterotype.

* [Uniform Manifold Approximation and Projection (UMAP)](https://github.com/lmcinnes/umap) is a dimension reduction technique that can be used for visualisation similarly to t-SNE, but also for general non-linear dimension reduction. The algorithm is founded on three assumptions about the data. Run from R: https://gist.github.com/crazyhottommy/caa5a4a4b07ee7f08f7d0649780832ef
* [umapr](https://github.com/ropenscilabs/umapr) UMAP dimensionality reduction in R
* [uwot](https://github.com/jlmelville/uwot) An R package implementing the UMAP dimensionality reduction method. UMAP multi-threaded.
* [Fast Fourier Transform-accelerated Interpolation-based t-SNE (FIt-SNE)](https://github.com/KlugerLab/FIt-SNE) The FIt-SNE implementation is generally faster than UMAP when you have more than 3,000 cells. In the realm of 10,000's of cells FIt-SNE scales at the same rate as UMAP. However, note that this is a log-log scale. Even if FI-tSNE starts scaling at the rate of UMAP, it is still consistently about 4 times faster. In other words, a dataset that takes an hour for UMAP will take 15 minutes for FIt-SNE. see the benchmark here https://nbviewer.jupyter.org/gist/vals/a138b6b13ae566403687a241712e693b by Valentine Svensson.
* [Parallel opt-SNE implementation with Python wrapper](https://github.com/omiq-ai/Multicore-opt-SNE) [preprint:Automated optimal parameters for T-distributed stochastic neighbor embedding improve visualization and allow analysis of large datasets](https://www.biorxiv.org/content/early/2018/11/23/451690.article-metrics)

### regulatory network

* [Scribe](https://github.com/cole-trapnell-lab/Scribe): Towards inferring causal regulations with single cell dynamics-coupled measurements
* single cell gene regulatory network analysis https://github.com/aertslab/SCENIC
* [Single-Cell Transcriptomics Unveils Gene Regulatory Network Plasticity](https://www.biorxiv.org/content/early/2018/10/17/446104)
### useful databases

* [CellMarker: a manually curated resource of cell markers in human and mouse](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gky900/5115823)
* [scRNAseq bioc package](https://bioconductor.org/packages/devel/data/experiment/html/scRNAseq.html) Gene-level counts for a collection of public scRNA-seq datasets, provided as SingleCellExperiment objects with cell- and gene-level metadata.
* [human cell atlas database](https://staging.data.humancellatlas.org/)
* [EMBL-EBI atlas](https://www.ebi.ac.uk/gxa/sc/home)
* (PanglaoDB)[https://panglaodb.se/) is a database for the scientific community interested in exploration of single cell RNA sequencing experiments from mouse and human. We collect and integrate data from multiple studies and present them through a unified framework.
* [scRNASeqDB])https://bioinfo.uth.edu/scrnaseqdb/)database, which contains 36 human single cell gene expression data sets collected from Gene Expression Omnibus (GEO)
* [JingleBell](http://jinglebells.bgu.ac.il/)A repository of standardized single cell RNA-Seq datasets for analysis and visualization at the single cell level.
* [Broad single cell portal](https://portals.broadinstitute.org/single_cell)
* The [conquer](http://imlspenticton.uzh.ch:3838/conquer/) (consistent quantification of external rna-seq data) repository is developed by Charlotte Soneson and Mark D Robinson at the University of Zurich, Switzerland. It is implemented in shiny and provides access to consistently processed public single-cell RNA-seq data sets. 

### interesting papers to read 

* [A single-cell molecular map of mouse gastrulation and early organogenesis](https://www.nature.com/articles/s41586-019-0933-9)
* [The single-cell transcriptional landscape of mammalian organogenesis](https://www.nature.com/articles/s41586-019-0969-x)
* [Comparative analysis of droplet-based ultra-high-throughput single-cell RNA-seq systems](https://www.biorxiv.org/content/early/2018/05/02/313130)
* [scRNA-seq mixology: towards better benchmarking of single cell RNA-seq protocols and analysis methods](https://www.biorxiv.org/content/early/2018/10/08/433102) github [repo](https://github.com/LuyiTian/CellBench_data)
* [A Single-Cell Transcriptome Atlas of the Aging Drosophila Brain](https://www.cell.com/cell/fulltext/S0092-8674(18)30720-7)
* [Bias, robustness and scalability in single-cell differential expression analysis](https://www.nature.com/articles/nmeth.4612) by Mark Robinson
* [Cell type transcriptome atlas for the planarian Schmidtea mediterranea](http://science.sciencemag.org/content/360/6391/eaaq1736)
* [Cell type atlas and lineage tree of a whole complex animal by single-cell transcriptomics](http://science.sciencemag.org/content/360/6391/eaaq1723)
* [Single-cell reconstruction of developmental trajectories during zebrafish embryogenesis](http://science.sciencemag.org/content/early/2018/04/25/science.aar3131)
* [The dynamics of gene expression in vertebrate embryogenesis at single-cell resolution](http://science.sciencemag.org/content/early/2018/04/25/science.aar5780)
* [Single-cell mapping of gene expression landscapes and lineage in the zebrafish embryo](http://science.sciencemag.org/content/early/2018/04/25/science.aar4362)
* [Decomposing cell identity for transfer learning across cellular measurements, platforms, tissues, and species]( https://www.biorxiv.org/content/early/2018/08/20/395004.1)
* [The contribution of cell cycle to heterogeneity in single-cell RNA-seq data](https://www.nature.com/articles/nbt.3498)
### merge different scRNAseq data sets

* [scMerge](https://www.biorxiv.org/content/early/2018/08/16/393280)

### single cell RNAseq copy-number variation
* [Linking transcriptional and genetic tumor heterogeneity through allele analysis of single-cell RNA-seq data.](https://www.ncbi.nlm.nih.gov/pubmed/29898899) tool [HoneyBADGER](https://github.com/JEFworks/HoneyBADGER)

### advance of scRNA-seq tech
* [Single-cell profiling of the developing mouse brain and spinal cord with split-pool barcoding](http://science.sciencemag.org/content/360/6385/176) no isolation of single cells needed!

* [Dynamics and Spatial Genomics of the Nascent Transcriptome by Intron seqFISH](https://www.cell.com/cell/fulltext/S0092-8674(18)30647-0)

* [Highly Multiplexed Single-Cell RNA-seq for Defining Cell Population and Transcriptional Spaces](https://www.biorxiv.org/content/early/2018/05/05/315333) blog post by Lior Patcher [The benefits of multiplexing](https://liorpachter.wordpress.com/2018/05/30/the-benefits-of-multiplexing/). Need to re-read carefully.

* [Three-dimensional intact-tissue sequencing of single-cell transcriptional states](http://science.sciencemag.org/content/early/2018/06/20/science.aat5691)

* [Multi-omic profiling of transcriptome and DNA methylome in single nuclei with molecular partitioning](https://www.biorxiv.org/content/early/2018/10/04/434845)

### Allele specific scRNAseq
* [scBASE](https://github.com/churchill-lab/scBASE) A set of tools for quantitation of allele-specific expression from scRNA-Seq data
* [paper: Genomic encoding of transcriptional burst kinetics](https://www.nature.com/articles/s41586-018-0836-1)

### pseudotemporal modelling

* [Uncovering pseudotemporal trajectories with covariates from single cell and bulk expression data](https://www.nature.com/articles/s41467-018-04696-6)

* all different algorithms https://github.com/agitter/single-cell-pseudotime
* [Genomic trajectories with heterogeneous genetic and environmental backgrounds](https://bioconductor.org/packages/release/bioc/html/phenopath.html)

* [A descriptive marker gene approach to single-cell pseudotime inference](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty498/5043298)

* [A collection of 50 trajectory inference methods within a common interface](https://github.com/dynverse/dynmethods) take a look of this!

* [velocyto](http://velocyto.org/) RNA abundance is a powerful indicator of the state of individual cells, but does not directly reveal dynamic processes such as cellular differentiation. Here we show that RNA velocity - the time derivative of RNA abundance - can be estimated by distinguishing unspliced and spliced mRNAs in standard single-cell RNA sequencing protocols. [paper](https://www.biorxiv.org/content/early/2017/10/19/206052) [comment](https://www.nature.com/articles/d41586-018-05882-8)

* [STREAM](http://stream.pinellolab.org/) is an interactive computational pipeline for reconstructing complex celluar developmental trajectories from sc-qPCR, scRNA-seq or scATAC-seq data from Luca Pinello Lab.

### large scale single cell analysis
* [bigSCale](https://genome.cshlp.org/content/early/2018/05/03/gr.230771.117.abstract): an analytical framework for big-scale single-cell data. [github link](https://github.com/iaconogi/bigSCale) for millions of cells (starts with a count matrix)
 [bigScale2](https://github.com/iaconogi/bigSCale2)
* [Alevin: An integrated method for dscRNA-seq quantification](https://www.biorxiv.org/content/early/2018/06/01/335000) based on Salmon.
* [How to Use Alevin with Seurat Alevin-Seurat Connection](https://combine-lab.github.io/alevin-tutorial/2018/alevin-seurat/) blog post
* [Kallisto BUStools](https://liorpachter.wordpress.com/2019/06/21/near-optimal-single-cell-rna-seq-pre-processing/) paper https://www.biorxiv.org/content/10.1101/673285v1
* [SCope: Visualization of large-scale and high dimensional single cell data](https://github.com/aertslab/SCope)
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

### some tools for 10x
[DropletUtils](https://www.bioconductor.org/packages/release/bioc/html/DropletUtils.html)
Provides a number of utility functions for handling single-cell (RNA-seq) data from droplet technologies such as 10X Genomics. This includes data loading, identification of cells from empty droplets, removal of barcode-swapped pseudo-cells, and downsampling of the count matrix.

