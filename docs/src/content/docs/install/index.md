---
title: "Download and install"
description: "Static binaries, Docker, bioconda, or build from source."
---

GLIMPSE2 runs on **Linux x86_64**, **Linux aarch64** and **macOS ARM64** (Apple Silicon). Pick whichever of the options below is easiest for you.

## Static binaries

Pre-built binaries are attached to every [release on GitHub](https://github.com/odelaneau/GLIMPSE/releases). They only support 64-bit Linux on x86_64 CPUs with AVX2. If your machine doesn't fit, use Docker or build from source.

## Bioconda

```bash
conda install -c bioconda glimpse-bio
```

## Docker

The repository has a `Dockerfile` that builds an image for `linux/amd64` and `linux/arm64` with all the GLIMPSE2 tools, `bcftools`, `bgzip` and `tabix`. The [Docker page](docker/) explains how to build and run it. Images of earlier releases are on [Docker Hub](https://hub.docker.com/r/simrub/glimpse/).

## Build from source

Install the libraries (HTSlib, Boost, libdeflate and a few others), clone the repository and run `make system`. The [Build from source](build/) page has the exact commands for Debian/Ubuntu, Fedora and macOS, and [Build the dependencies yourself](dependencies/) covers systems where you can't install packages or where the available versions are too old.
