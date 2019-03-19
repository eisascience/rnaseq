from bioconductor/release_core2
RUN apt-get update -y && apt-get install -y libhdf5-dev pip
RUN pip install umap-learn
ADD install.R /
RUN Rscript install.R