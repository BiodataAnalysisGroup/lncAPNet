# Nextflow pipeline

## Run the pipeline

Define the local path for `input` in the `conf/local.config` file.
In the same file, update the local path to where the singularity (or docker) containers are stored.

The pipeline should then run using `docker` or `singularity` containers with this command
from within the `Nextflow/pipeline` folder:
```
nextflow run main.nf -c ../conf/local.config -profile singularity,local --outdir output
```
