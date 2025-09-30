# data

The folders inside this folder should contain all data at various stages.

This data is being loaded/manipulated/changed/saved with code from the `code` folders.

## `data/raw`

The subdirectory `data/raw` contains the original input files which are never
edited manually, they are just used as inputs to our code.

- `clean-data.Rds`: the UGAFluVac data
- `joined-data.csv`: the cohort data with distances from the BQ repo, see
script 01 for details.

## `data/processed`

The subdirectory `data/processed` contains the output files from the data
cleaning operations. Only the script `01-Data-Processing.R` should output
files into this directory.

- `model-data.Rds` contains all of the data that will be used for creating our
statistical models in this repo.
- `reporting-data.Rds` contains all of the data that will need to be reported
to NIH CIVICs for submission as metadata on the CIVICs portal.

We recommend using the `.Rds` files to ensure proper data formatting, and those
are the only files that are used in our downstream pipeline. However, we also
provide `.csv` files of these datasets for ease of use.

<!-- end of file -->
