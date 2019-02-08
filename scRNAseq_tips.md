
### markers for immune cells
https://twitter.com/lcymcdnld/status/1093650424303824897

* CD4 is not well captured (drop out) by 10x 3' technology. 5' is better.
* NKG7 is typically well captured and present in CD8 and NK cells but not in CD4
* CD3D, IL7R (CD4 T cells). CD8A. GNLY, NKG7 (NK). CD14, LYZ (Monocytes).
* This paper also presents the markers for each immune cell type https://www.ncbi.nlm.nih.gov/pubmed/29227470

### imputation caveat

https://twitter.com/jamesrhowe6/status/1073644734323683328

>James: Probably not a good idea to impute single cell RNA-seq data... neuron-specific marker VGlut2/Slc17a6 is expressed at low levels in every single CNS cell cluster after SAVER imputation

>George: This was actually our motivation for the development of [ALRA](https://www.biorxiv.org/content/10.1101/397588v1): to create a method for imputation of scRNA-seq data that preserves biological zeros.

>Lana: Please try our tool [deepImpute](https://www.biorxiv.org/content/10.1101/353607v1). Wispering: we don't think it makes much sense to impute genes with very low and extremely sparse expression to begin with; and deepImpute is orders faster than tools like SAVER.
