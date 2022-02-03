# DeWolf, S. et. al.; _Tissue-specific landscape of the human T cell repertoire in graft-versus-host-disease_

# Cloning this repository

As this repo contains submodule, please clone as follows:

    git clone --recurse-submodules https://github.com/vdblab/Dewolf-TCR-landscape-2022

## Accessing Data

Raw TCR sequencing data can be obtained \<adaptivebiotech.com\> under project ID XXXX. Each sample data can be exported by clicking the checkbox and then export-\> Export. Save data to `./raw_TCR_data/`

### Patient ID's

For certain analyses (such as GLIPH), id's of samples and patients needed to be changed. To avoid confusion, a key of patient names is included here: `data/pt_key_forGLIPH.csv`.

# GLIPH Analysis

## Data cleaning

Sample id's were changed to conform to the format needed by gliph, and processed with the script `04_gliph2.Rmd`. The output of that script with the relevant columns is written out to `all_tcr_data_withcomp.tsv`.

## Running GLIPH

The parameters for gliph can be found `GLIPH/test_paramter.txt`. GLIPH version 0.01 was used.

# Calculating Overlap Score of GLIPH results

The script `overlap_score_calculations.py` reads in `test_cluster_allpts_withcomparators` to create the overlap figure.

    # if needed 
    # conda create -n tcr --file conda_env.txt
    # conda activate tcr

    python ./scripts/overlap_score_calculations.py ./data/test_cluster_allpts_withcomparators.csv

## Halo Analysis

See soccin's repository for a description of this analysis. That repo is included as a submodule.

## GTEx comparison

The comparison between the TCR profiles and the GTEx transcription data are in `GTEx/`.  The documents there describe downloading, preprocessing the data, and generating the figures.

