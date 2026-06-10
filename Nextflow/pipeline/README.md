# Nextflow pipeline

## Run the pipeline

Define the local path for `input` in the `conf/local.config` file.
If using singularity, in the same file, uncomment and update the local path to where the singularity containers are stored.

The pipeline should then run using `docker` or `singularity` containers with this command
from within the `Nextflow/pipeline` folder:
```
nextflow run main.nf -c ../conf/local.config -profile singularity,local --outdir output
```
or, directly passing the input to the CLI command:
```
nextflow run main.nf -profile singularity --input /path/to/data/test_eset.rds --gene_info /path/to/data/gene_info.xlsx --outdir output
```
or, with the pre-made test profile:
```
nextflow run main.nf -profile singularity,test --outdir output -resume


nextflow run main.nf -profile docker --input ../Nextflow_demo/data/test_eset.rds --outdir output --gene_info ../Nextflow_demo/data/gene_info.xlsx --group0 "M" --group1 "U" --comparison "IGHV" -resume

```
