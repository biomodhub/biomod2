# syntax=docker/dockerfile:1
# This works:
# FROM ghcr.io/osgeo/gdal:ubuntu-small-3.9.3
FROM ghcr.io/osgeo/gdal:ubuntu-small-3.12.1

RUN export DEBIAN_FRONTEND=noninteractive && \
apt-get -y update && apt-get -y upgrade && \
apt-get -y install \
  ca-certificates \
  dirmngr \
  gpg \
  software-properties-common \
  dos2unix \
  git \
  build-essential \
  cmake \
  libcurl4-openssl-dev \
  libgsl0-dev \
  libjansson-dev \
  libssl-dev \
  libtiff-dev \
  libudunits2-dev \
  libproj-dev \
  libsqlite3-dev \
  libgeos-dev \
  lockfile-progs \
  rename \
  libcurl4-openssl-dev \
  python3-pip \
  python-is-python3 \
  parallel \
  r-base \
  wget \
  aria2

RUN pip3 install --break-system-packages --no-cache-dir \
    numpy \
    gsutil \
    scipy && \
rm /usr/local/lib/pkgconfig/proj.pc && \
Rscript -e 'install.packages(c("s2", "proxy", "e1071", "classInt", "DBI", "magrittr", "units"), repos="https://cloud.r-project.org"); if (!library(s2, logical.return=T)) quit(save="no", status=10)' && \
Rscript -e 'install.packages("sf", verbose=T, INSTALL_opts=c("--no-test-load"), repos="https://cloud.r-project.org")'
