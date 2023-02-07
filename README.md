<p align="center">
  <a href="https://odelaneau.github.io/GLIMPSE/">
    <img height="150" src="docs/assets/images/branding/glimpse_logo_400x171.png">
  </a>
  <h1 align="center">GLIMPSE</h1>
</p>

[![GitHub license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/odelaneau/GLIMPSE/blob/master/LICENSE)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/glimpse-bio/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/glimpse-bio/badges/downloads.svg)](https://anaconda.org/bioconda/glimpse-bio)
[![Docker hub](https://img.shields.io/docker/pulls/simrub/glimpse)](https://hub.docker.com/r/simrub/glimpse/)
[![GitHub Releases](https://img.shields.io/github/v/release/odelaneau/GLIMPSE)](https://github.com/odelaneau/GLIMPSE/releases)


### Versions

**Current release: v2.0.0**. Release date: Dec 07, 2022

For details of past changes please see [CHANGELOG](docs/CHANGELOG.md).

### License

GLIMPSE is available under a MIT license. For more information please see the [LICENSE](LICENSE).
 
### Features
GLIMPSE2 is a set of tools for phasing and imputation for low-coverage sequencing datasets:

- **GLIMPSE2_chunk** splits the genome into chunks for imputation and phasing
- **GLIMPSE2_split_reference** creates the reference panel representation used by GLIMPSE2_phase. It allows major speedups for large reference panels.
- **GLIMPSE2_phase** program to impute and phase low coverage sequencing data
- **GLIMPSE2_ligate** concatenates phased chunks of data into chromosome-wide phased files

### Documentation and tutorials

Visit our website for tutorials, documentation and installation instructions:

https://odelaneau.github.io/GLIMPSE/

### GLIMPSE1

At the moment, GLIMPSE2 performs imputation only from a reference panel of samples. To use the joint-model, particularly useful for many samples at higher coverages (>0.5x) and small reference panels, please visit the GLIMPSE1 website and checkout the [GLIMPSE1 branch](https://github.com/odelaneau/GLIMPSE/tree/glimpse1).


### Installation

To build the source code, please refer to the [step-by-step guide on the website](https://odelaneau.github.io/GLIMPSE/docs/installation).

To run and test the pre-made docker images, you can run the following:



#### GLIMPSE v2.0.0:
```
wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.0/glimpse_v2.0.0-27-g0919952_20221207.tar.gz
docker load < glimpse_v2.0.0-27-g0919952_20221207.tar.gz
# test docker run using -it
docker run -it glimpse:v2.0.0-27-g0919952_20221207
```
To run the tools you can use 'GLIMPSE2_[tool name]' (e.g. GLIMPSE2_phase) inside the docker run (-it) session 

#### GLIMPSE v1.1.1:
```
wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/glimpse_v1.1.1-c27e90d_20210521.tar.gz
docker load < glimpse_v1.1.1-c27e90d_20210521.tar.gz
# test docker run using -it
docker run -it glimpse:v1.1.1-c27e90d_20210521
```
To run the tools you can use 'GLIMPSE_[tool name]_v1.1.1' (e.g. GLIMPSE_phase_v1.1.1) inside the docker run (-it) session 
