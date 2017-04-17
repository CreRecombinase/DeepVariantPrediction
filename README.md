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

To run through the variant2input step:
```
$ snakemake -r input/pranav_test1_allele1.hdf5 input/pranav_test2_allele1.hdf5
```

To run through both variant2input and input2score steps:
```
$ snakemake -r score/keras_deepsea_copy/pranav_test1_allele1.hdf5 score/keras_deepsea_copy/pranav_test1_allele2.hdf5
```

You will get:
```
$ h5dump -A -H score/keras_deepsea_copy/pranav_test1_allele1.hdf5 score/keras_deepsea_copy/pranav_test1_allele2.hdf5
HDF5 "score/keras_deepsea_copy/pranav_test1_allele1.hdf5" {
GROUP "/" {
   DATASET "y_pred" {
      DATATYPE  H5T_IEEE_F32LE
      DATASPACE  SIMPLE { ( 9, 919 ) / ( 9, 919 ) }
   }
}
}
HDF5 "score/keras_deepsea_copy/pranav_test1_allele2.hdf5" {
GROUP "/" {
   DATASET "y_pred" {
      DATATYPE  H5T_IEEE_F32LE
      DATASPACE  SIMPLE { ( 9, 919 ) / ( 9, 919 ) }
   }
}
}
```

# Structure of the pipeline

The pipeline contains three modules:

1. `modules/variant2input.snakemake`: takes a variant file and output a hdf5 file that is ready for being used by `keras` model.
    - input: specifies at `config.yaml`
    - output: `input/{data}_allele1.hdf5` and `input/{data}_allele2.hdf5`
2. `modules/input2score.snakemake`: takes the two outputs above and the model specified at `config.yaml`, output the predicted scores.
    - input: see above and model
    - output: `score/[model_name]/[data]_allele1.hdf5` and `score/[model_name]/[data]_allele2.hdf5`
3. `modules/score2performance.snakemake`: takes the predictions along with the group of labels of interest, output data table ready for downstream analysis in RDS format (combining all input variant sets with selected labels) and a quick summary per group of labels and summary method (specified in `config.yaml`).
    - input: see above and groups of labels
    - output: RDS file per group `score/[model]/[group]/result.rds` and `performance/[model]/[group]/[method]/report.html`

# TODO

1. Add more formatting method at `modules/submodules/variant2input/formatting`. Now only Pranav data format has been implemented (from GWAS study, see format at `test/`).
2. ~Debug `modules/input2score.snakemake`~
3. Debug `modules/score2performance.snakemake` first part and implement at least one performance report script
4. Add training snakemake script

# Attention

1. Now we only implemented a very ugly "do nothing" train model sub workflow at `modules/submodules/input2score/_train_model_do_nothing.snakemake`. It takes a pretrained keras model and copy it as `{model}_copy.hdf5`.
2. BED.formatted requirement (strand)
