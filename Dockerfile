from bioconductor/release_core2
RUN apt-get -q -y install ed libcairo2-dev libxt-dev && apt-get clean
ADD install.R /
RUN Rscript install.R