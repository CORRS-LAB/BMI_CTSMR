# BMI-CTS-MR

EWAS analysis code for cell-type-specific associations between DNA methylation regions and body mass index.

## Files

- `scripts/ewas.R`: discovery and validation models, discovery BH correction within each cell type, and association tables.
- `data/release/`: validation methylation data for the 66 selected hg19 regions, phenotype, and cell proportions.
- `results/ewas/`: model results. The complete discovery statistics are supplied as `discovery_all.csv.gz`.
- `FILE_MANIFEST.csv`: file checksums.

Discovery individual-level methylation, phenotype, and cell-proportion data are not distributed.

## Run the analysis

Install the required R packages:

```r
install.packages(c("data.table", "plyr"))
```

Set the R working directory to this repository's root. Provide the discovery inputs locally:

```text
data/release/discovery_model_phenotype.csv
data/release/discovery_cell_proportions.csv
data/prepared/discovery_pre_QC/discovery.pre_QC.beta.tsv
```
(The file `discovery.pre_QC.beta.tsv` can be downloaded from GEO using the link provided in the paper.)

Open `scripts/ewas.R` in RStudio and run it, or run:

```r
source("scripts/ewas.R")
```

The script writes five CSV files under `results/ewas/`.
## License

[Apache 2.0](LICENSE).
