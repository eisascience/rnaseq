from rocker/verse:3.5.2
RUN apt-get update -y && apt-get install -y libhdf5-dev
ADD install.R /
RUN Rscript install.R