## conda environment
conda create --name scte r=3.6
conda activate scte
conda isntall python=3.8

## just for jupyter
conda install -y jupyter
python -m ipykernel install –name scte –display-name "scte" # for using conda env in jupyter



conda install matplotlib
#conda install -c bioconda anndata2ri
pip install anndata2ri
conda install -c conda-forge leidenalg


##=== imortmant ===##
pip install scanpy
conda install -c conda-forge rpy2


## R package scran
Rscript install_scran.R


## conda for R package
conda install bioconductor-rhdf5=2.30.0
conda install bioconductor-rhdf5lib=1.8.0
