## conda environment
conda create --name scte r=3.6
conda activate scte
conda isntall python=3.8

## just for jupyter
conda install -y jupyter


conda install matplotlib
#conda install -c bioconda anndata2ri
pip install anndata2ri
conda install -c conda-forge leidenalg


##=== imortmant ===##
pip install scanpy
conda install -c conda-forge rpy2


## R package scran
Rscript install_scran.R