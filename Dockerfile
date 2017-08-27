from bioconductor/release_core2
RUN apt-get install ed && apt-get clean
ADD install.R /
RUN Rscript install.R