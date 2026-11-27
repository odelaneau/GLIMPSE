<p align="center">
  <a href="https://odelaneau.github.io/GLIMPSE/">
    <img height="150" src="docs/assets/images/branding/glimpse_logo_400x171.png">
  </a>
  <h1 align="center">GLIMPSE2</h1>
  <p align="center">Phasing and imputation of low-coverage whole-genome sequencing data, designed for reference panels with hundreds of thousands of samples.</p>
</p>

[![GitHub Releases](https://img.shields.io/github/v/release/odelaneau/GLIMPSE?label=release)](https://github.com/odelaneau/GLIMPSE/releases)
[![Release date](https://img.shields.io/github/release-date/odelaneau/GLIMPSE?label=released)](https://github.com/odelaneau/GLIMPSE/releases)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/glimpse-bio/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/glimpse-bio/badges/downloads.svg)](https://anaconda.org/bioconda/glimpse-bio)
[![Docker hub](https://img.shields.io/docker/pulls/simrub/glimpse)](https://hub.docker.com/r/simrub/glimpse/)
[![GitHub license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/odelaneau/GLIMPSE/blob/master/LICENSE)

**[Website](https://odelaneau.github.io/GLIMPSE/)** · **[Install](https://odelaneau.github.io/GLIMPSE/install/)** · **[Tutorial](https://odelaneau.github.io/GLIMPSE/tutorial/getting-started/)** · **[Documentation](https://odelaneau.github.io/GLIMPSE/tools/)** · **[Changelog](https://odelaneau.github.io/GLIMPSE/changelog/)**

### Tools

GLIMPSE2 takes a reference panel and low-coverage BAM/CRAM files, or genotype likelihoods, and returns phased and imputed genotypes. It is made of these tools:

- **GLIMPSE2_chunk** splits the genome into chunks for imputation and phasing
- **GLIMPSE2_split_reference** creates the reference panel representation used by GLIMPSE2_phase. It allows major speedups for large reference panels.
- **GLIMPSE2_phase** imputes and phases low-coverage sequencing data
- **GLIMPSE2_ligate** concatenates phased chunks of data into chromosome-wide phased files
- **GLIMPSE2_concordance** measures imputation accuracy against high-coverage genomes
- **GLIMPSE2_inspect** summarises a binary reference panel

Tutorials, documentation and installation instructions are on the [website](https://odelaneau.github.io/GLIMPSE/). The release history is in the [changelog](https://odelaneau.github.io/GLIMPSE/changelog/).

### GLIMPSE1

GLIMPSE2 imputes from a reference panel of samples only. For the joint model, useful with many samples at higher coverage (>0.5x) and small reference panels, see the [GLIMPSE1 website](https://odelaneau.github.io/GLIMPSE/glimpse1/index.html) and the [GLIMPSE1 branch](https://github.com/odelaneau/GLIMPSE/tree/glimpse1).

### Building from source

GLIMPSE2 supports **Linux x86_64**, **Linux aarch64**, and **macOS ARM64 (Apple Silicon)**.

<details>
<summary><b>macOS (Apple Silicon)</b></summary>

```bash
brew install htslib boost libdeflate
git clone --recursive https://github.com/odelaneau/GLIMPSE.git
cd GLIMPSE
make system
```
</details>

<details>
<summary><b>Linux (Debian/Ubuntu)</b></summary>

```bash
sudo apt install build-essential libboost-iostreams-dev libboost-program-options-dev libboost-serialization-dev libhts-dev libdeflate-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev zlib1g-dev
git clone --recursive https://github.com/odelaneau/GLIMPSE.git
cd GLIMPSE
make system
```

On aarch64, also install `clang`: the build uses it by default there because it makes GLIMPSE2_phase much faster.
</details>

For detailed instructions including building dependencies from source, see the [step-by-step guide on the website](https://odelaneau.github.io/GLIMPSE/install/build/).

#### Docker images

##### Building the image
The `Dockerfile` builds an image with all the GLIMPSE2 tools plus bcftools, bgzip and tabix, for linux/amd64 and linux/arm64. `docker-build.sh` builds it from a commit or tag:
```
./docker-build.sh v2.1.0                             # glimpse:2.1.0 for this machine
./docker-build.sh v2.1.0 --push ghcr.io/me/glimpse   # both architectures, pushed to a registry
```
See the [Docker page on the website](https://odelaneau.github.io/GLIMPSE/install/docker/) for building by hand and the build options.

##### Dockerhub
You can find GLIMPSE on docker hub at this link:
https://hub.docker.com/r/simrub/glimpse/

##### Local download from Github releases

###### GLIMPSE v2.0.0:
```
wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.0/glimpse_v2.0.0-27-g0919952_20221207.tar.gz
docker load < glimpse_v2.0.0-27-g0919952_20221207.tar.gz
# test docker run using -it
docker run -it glimpse:v2.0.0-27-g0919952_20221207
```
To run the tools you can use 'GLIMPSE2_[tool name]' (e.g. GLIMPSE2_phase) inside the docker run (-it) session 

###### GLIMPSE v1.1.1:
```
wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/glimpse_v1.1.1-c27e90d_20210521.tar.gz
docker load < glimpse_v1.1.1-c27e90d_20210521.tar.gz
# test docker run using -it
docker run -it glimpse:v1.1.1-c27e90d_20210521
```
To run the tools you can use 'GLIMPSE_[tool name]_v1.1.1' (e.g. GLIMPSE_phase_v1.1.1) inside the docker run (-it) session 

### Citation

If you use GLIMPSE2, please cite:

> Rubinacci S, Hofmeister RJ, Sousa da Mota B, Delaneau O. **Imputation of low-coverage sequencing data from 150,119 UK Biobank genomes.** *Nature Genetics* 55, 1088-1090 (2023). [doi:10.1038/s41588-023-01438-3](https://doi.org/10.1038/s41588-023-01438-3)

### Acknowledgements

GLIMPSE is developed and maintained by Simone Rubinacci (University of Helsinki) and Olivier Delaneau (Regeneron).

We thank Tim Fennell and Fulcrum Genomics for the performance, memory, build and Docker work in 2.1, and Chris Kachulis and Michael Gatzen (Broad Institute) for cloud-computing support (checkpointing and retries) and bug fixes.

### License

GLIMPSE is available under the MIT license. For more information please see the [LICENSE](LICENSE).
