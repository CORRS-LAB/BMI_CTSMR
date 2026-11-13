# BMI-CTS-MR

Code to identify immune cell type-specific DNA methylation regions asscociated with BMI.

## Files

- `code/ewas.R`: discovery and validation models, discovery BH correction within each cell type, and association tables.
- `data/`: discovery cell proportions, plus validation methylation data for the 66 selected regions, phenotype, and cell proportions.
- `results/`: model results. The complete discovery statistics are supplied as `discovery_all.csv.gz`.

The discovery phenotype and methylation matrix should be obtained from GEO using the link provided in the paper; these two files are not included in this repository.

## Run the analysis

Install the required R packages:

```r
install.packages(c("data.table", "plyr"))
```

Set the R working directory to this repository's root. Obtain the discovery phenotype and methylation matrix from GEO using the link in the paper, and save them at the following paths with sample IDs matching the supplied discovery cell-proportion table:

```text
data/discovery_model_phenotype.csv
data/discovery.pre_QC.beta.tsv
```
Open `code/ewas.R` in RStudio and run it, or run:

```r
source("code/ewas.R")
```

The script writes five CSV files under `results/`.
## License

[Apache 2.0](LICENSE).
