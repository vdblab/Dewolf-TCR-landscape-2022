## README

### prerequisites
- a docker installation

### Running Jupyter
A jupyterlab server can be started with the following command in the terminal to pull the datascience notebook base image from dockerhub and spin up a local server.
```
docker run -it --rm -p 10000:8888 -v "${PWD}":/home/jovyan/work  jupyter/datascience-notebook
```

## `sc-tcr-analysis.ipynb`
Running the notebook will install the required packages, download the raw data, and generate several files summarizing the TCR repertoire of 3 samples.


## Nuc-Seq analysis
 
The Nuc Seq Jupyter notebook contains the code used to run the Nuc Seq processing and downstream analysis. This notebook has been adapted from the public shunPykeR GitHub [repo](https://github.com/kousaa/shunPykeR). Please use instructions [here](https://github.com/kousaa/shunPykeR) to install the shunPykeR.yml environment to fully reproduce the results of this pipeline.
