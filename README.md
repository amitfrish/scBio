# scBio
scBio is a R pacakge containing a repository of methods aiming to understand complex processes related to the cellular variability within tissues.

### New in version 0.1.5:
1. Added automatic alteration of the neighbourhood size parameter in the cases of cell types with limited number of cells.
2. Fixed a problem when using a single cell or bulk data with small amount of genes.

## CPM
A method based on computational deconvolution for identifying a cell population map from bulk gene expression data of a heterogeneous sample. The CPM method provides an advantageous alternative to existing deconvolution approaches, particularly in providing a fine-resolution mapping. Specifically, CPM is focused on cell alterations within each cell type and not changes in the total number of cells in the cell type. Therefore, it can model cell changes across trajectories and specific cell subtypes.

### Inputs
#### CPM has four mandetory inputs:

SCData - A matrix containing the single-cell RNA-seq data. Each row corresponds to a certain gene and each column to a certain cell. Importantly, CPM relies on many iterative processes and therefore might take a long running time. For extremely large single cell datasets, we suggest to use only a portion of the data, using random sampling of the cells.

SCLabels - A vector containing the labels of each of the cells.

BulkData - A matrix containing heterogenous RNA-seq data for one or more samples. Each row corresponds to a certain gene and each column to a certain sample.

cellSpace - The cell state space corresponding to the single-cell RNA-seq data. It can be a vector for a 1-dim space or a matrix for a multidimensional space where each column represents a different dimension.

```
data(SCLabels)
data(SCFlu)
data(BulkFlu)
data(SCCellSpace)
```

#### CPM has also two important optional inputs:

neighborhoodSize - Cell neighborhood size which will be used for the analysis. This should be lower than the number of cells in the smallest cell type. The defalt is 10. Generally, neighborhoodSize should be the highest number which represents well the number of neighbouring cells that are similar to the center cell. Large number can increase the analysis robustness in dense cell spaces but also can reduce it if the space is too sparse. 

modelSize - The reference subset size in each iteration of CPM. This should be lower than the total number of cells. The defalt is 50. The selected cells within each model will be gathered from all cell types. If you have a low number of cell types or low totel number of cells across all cell types, this should be lower than 50. Generally, up to five or six cells from each cell type is still okay.

#### Other optional inputs:

no_cores - A number for the amount of cores which will be used for the analysis. The defalt (NULL) is total number of cores minus 1.

minSelection - The minimum number of times in which each reference cell is selected. Increasing this value might make CPM more robust but also have a large effect on the algorithm's running time. The defalt is 5.

quantifyTypes - A boolean parameter indicating whether the prediction of cell type quantities is needed. This is recommended only in the case of homogenicity within cell types. Cell types with high inner cellular variability will recieve less reliabe values. The default is FALSE.

typeTransformation - This parameter will have an effect only if quantifyTypes = TRUE. A boolean parameter indicating whether cell type deconvolution should be provided in fractions. This is done by substracting all cell types by values of the minimal cell type and dividing in their sum. This is not recommended, since it reduces the comparability between sample. The default is FALSE.

calculateCI - A boolean parameter indicating whether the calculation of confidence itervals is needed. The default is FALSE.

### Relative prediction
If you have a test group and a control group, for example: disease and healthy, different time points and time points zero, you should use CPM in a relative manner.

To do so, make sure your input bulk and single cell data sets are in log-scale (already logged in the example data) and substract the mean of the control group for all samples in the test group.

```
BulkFluReduced = BulkFlu - rowMeans(BulkFlu[,grep("pbs",colnames(BulkFlu))])
BulkFluReduced = BulkFluReduced[,grep("flu",colnames(BulkFluReduced))]
```

Then, run CPM using the new bulk data set:

```
res = CPM(SCFlu, SCLabels, BulkFluReduced, SCCellSpace, no_cores = 6)
```

### Absolute prediction
Use input bulk and single cell data sets are in a linear-scale.

```
BulkFluAbs = exp(BulkFlu)-1
SCFluAbs = exp(SCFlu)-1
```

Then, run CPM using the new bulk data set:

```
resAbs = CPM(SCFluAbs, SCLabels, BulkFluAbs, SCCellSpace, no_cores = 6)
```

### Cell type prediction
Cell type predition is only a side result in this algorithm but it can be done as well.

Using relative data:
```
res = CPM(SCFlu, SCLabels, BulkFluReduced, SCCellSpace, quantifyTypes = T, no_cores = 6)
```

Using absolute data:
```
resAbs = CPM(SCFluAbs, SCLabels, BulkFluAbs, SCCellSpace, quantifyTypes = T, typeTransformation = T, no_cores = 6)
```

### Outputs
predicted	- CPM predicted cell abundance matrix (The main result of CPM). Each row represents a sample and each column a single cell.

cellTypePredictions	 - CPM predicted cell-type abundance matrix. Each row represnts a sample and each column a single cell-type. This is calculated if quantifyTypes = TRUE.

confIntervals	- A matrix containing the confidence iterval for each cell and sample. Each row represnts a sample and each column a single cell. This is calculated if calculateCI = TRUE.

numOfRuns	- The number of deconvolution repeats preformed by CPM.
