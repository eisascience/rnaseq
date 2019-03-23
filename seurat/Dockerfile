from bioconductor/release_core2
RUN apt-get update -y && apt-get install -y libhdf5-dev python-pip
RUN pip install numba==0.42.0
RUN pip install umap-learn
ADD install.R /
RUN Rscript install.R