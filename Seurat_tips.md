

### Merge a lot of Seurat objects

load libraries 

```r
#install.packages("tidyverse")
#install.packages("Seurat")
# here package is to get rid of absolute path
# see details https://github.com/jennybc/here_here
#install.packages("here")
# interact with file systems https://github.com/r-lib/fs
#install.packages("fs")
library(tidyverse)
library(Seurat)
library(here)
library(fs)
library(purrr)
here()
```

read in the data as seruat object
and then merge https://github.com/satijalab/seurat/issues/480

```{r}
input_folders<- dir_ls( path = here("data"), recursive = T) %>% path_dir() %>%
        unique() %>% str_subset(""mm10-1.2.0_premrna")

samples<- str_extract(input_folders, pattern = yourpattern)
names(input_folders)<- samples
seurat_data<- purrr::map(input_folders, Read10X)

#prefix the sample name to the cell name, otherwise merge seurat objects gives error
add_sample_name_to_cell<- function(x, y){
        colnames(x)<- paste(y, colnames(x), sep = "_")
        return(x)
}

seurat_data<- map2(seurat_data, samples, add_sample_name_to_cell)
seurat_objects<- map2(seurat_data, samples, function(x,y) CreateSeuratObject(raw.data = x, project = y))
```

merge

```{r}

merged_seurat<- purrr::reduce(seurat_objects, function(x,y) {MergeSeurat(x,y, do.normalize = F)})
```

### SubsetData

```
When we run SubsetData, we have (by default) not subsetted the raw.data slot as well, as this can be slow and usually unnecessary. This can in some cases cause problems downstream, but setting do.clean=T does a full subset.
```
