---
layout: default
title: System requirements
nav_order: 1
parent: Build from source
grand_parent: Installation
permalink: /docs/installation/build_from_source/system_requirements
---
# System requirements
{: .no_toc .text-center }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## System requirements
GLIMPSE2 is a set of C++ tools covering the process of haplotype phasing in large datasets. GLIMPSE2 supports the following platforms:

- **Linux x86_64** (GCC or Clang)
- **Linux aarch64/ARM64** (GCC or Clang)
- **macOS ARM64 / Apple Silicon** (Clang)

### Linux

We require a version of GCC > 4.4 or a recent version of Clang. We recommend using the latest available version for your system.

{: .note }
> **On aarch64/ARM64, Clang is strongly recommended and is the build's default compiler.** Clang generates much faster NEON code for GLIMPSE2's SIMDe-translated kernels than GCC does (~2.5× faster phase imputation on AWS Graviton / Neoverse). Install it with `sudo apt install clang` (Debian/Ubuntu) or `sudo dnf install clang` (Fedora/RHEL). `make system` uses `clang++` automatically when present and falls back to `g++` otherwise. On x86_64, GCC remains the default.

For example running the following instruction on Ubuntu 20.04 focal:

<div class="code-example" markdown="1">
```bash
sudo apt install build-essential
```
</div>

will install the GNU g++ compiler version 9.2. To check the version of your g++ compiler, simply run:

<div class="code-example" markdown="1">
```bash
g++ --version
```
</div>

### macOS (Apple Silicon)

macOS builds use the Apple Clang compiler included with Xcode Command Line Tools. To install:

<div class="code-example" markdown="1">
```bash
xcode-select --install
```
</div>

To check your compiler version:

<div class="code-example" markdown="1">
```bash
clang++ --version
```
</div>


