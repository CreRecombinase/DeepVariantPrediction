# Prepare virtual environment

Before getting prepared, please make sure that you have install `conda` on your computer.

```
$ cd DeepVariantPrediction
$ conda env create --name deepvarpred --file environment.yaml
$ source activate deepvarpred
```

If you have already set up the virtual environment, run

```
$ source activate deepvarpred
```

# Run a test

To see the schedule:
```
$ snakemake -np input/pranav_test1_allele1.hdf5 input/pranav_test2_allele1.hdf5
```

To run through the pipeline:
```
$ snakemake -r input/pranav_test1_allele1.hdf5 input/pranav_test2_allele1.hdf5
```
