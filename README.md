# scBio
scBio is a R pacakge containing a repository of methods aiming to understand complex processes related to the cellular variability within tissues.

### CPM
A method based on computational deconvolution for identifying a cell population map from bulk gene expression data of a heterogeneous sample. The CPM method provides an advantageous alternative to existing deconvolution approaches, particularly in providing a fine-resolution mapping. Specifically, CPM is focused on cell alterations within each cell type and not changes in the total number of cells in the cell type. Therefore, it can model cell changes across trajectories and specific cell subtypes.

#### Inputs
CPM has four mandetory inputs:

SCData - A matrix containing the single-cell RNA-seq data. Each row corresponds to a certain gene and each column to a certain cell.

SCLabels - A vector containing the labels of each of the cells.

BulkData - A matrix containing heterogenous RNA-seq data for one or more samples. Each row corresponds to a certain gene and each column to a certain sample.

cellSpace - The cell state space corresponding to the single-cell RNA-seq data. It can be a vector for a 1-dim space or a matrix for a multidimensional space where each column represents a different dimension.

```
data(SCLabels)
data(SCFlu)
data(BulkFlu)
data(SCCellSpace)
```

#### Relative prediction
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

#### Absolute prediction
Use input bulk and single cell data sets are in a linear-scale.

```
BulkFluAbs = exp(BulkFlu)-1
SCFluAbs = exp(SCFlu)-1
```

Then, run CPM using the new bulk data set:

```
resAbs = CPM(SCFluAbs, SCLabels, BulkFluAbs, SCCellSpace, no_cores = 6)
```

#### Cell type prediction
Cell type predition is only a side result in this algorithm but it can be done as well.

Using relative data:
```
res = CPM(SCFlu, SCLabels, BulkFluReduced, SCCellSpace, quantifyTypes = T, no_cores = 6)
```

Using absolute data:
```
resAbs = CPM(SCFluAbs, SCLabels, BulkFluAbs, SCCellSpace, quantifyTypes = T, typeTransformation = T, no_cores = 6)
```

#### Outputs
predicted	- CPM predicted cell abundance matrix (The main result of CPM). Each row represents a sample and each column a single cell.

cellTypePredictions	 - CPM predicted cell-type abundance matrix. Each row represnts a sample and each column a single cell-type. This is calculated if quantifyTypes = TRUE.

confIntervals	- A matrix containing the confidence iterval for each cell and sample. Each row represnts a sample and each column a single cell. This is calculated if calculateCI = TRUE.

numOfRuns	- The number of deconvolution repeats preformed by CPM.
