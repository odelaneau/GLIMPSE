---
layout: default
title: Compile GLIMPSE2
parent: Build from source
grand_parent: Installation
permalink: /docs/installation/build_from_source/compile_glimpse2
---
# Compile GLIMPSE2
{: .no_toc .text-center }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Compile GLIMPSE2
Download the last version of the GLIMPSE2 code using:
<div class="code-example" markdown="1">
```bash
git clone --recursive https://github.com/odelaneau/glimpse.git
```
</div>

The `--recursive` flag is required to fetch the [SIMDe](https://github.com/simd-everywhere/simde) submodule used for portable SIMD support. If you have already cloned without `--recursive`, run:

<div class="code-example" markdown="1">
```bash
git submodule update --init
```
</div>

Navigate to the downloaded folder using `cd glimpse`.

You'll find there a folder containing all the software packages and other utility folders:

- **chunk**: program to define chunks where to run phasing and imputation.
- **common**: basic source files used by different tools.
- **concordance**: program to verify the accuracy of low-coverage imputation against high-coverage genomes
- **docker**: all scripts needed to build a docker file comprising all binaries
- **docs**: documentation in html
- **ligate**: ligate multiple imputed BCF/VCF files into a single chromosome-length file
- **maps**: genetics maps in b37 and b38
- **phase**: program to impute and phase low-coverage data.
- **split_reference**: program to create GLIMPSE2's reference file format, used by GLIMPSE2_phase
- **tutorial**: test datasets and scripts
- **versions**: versioning

Each program in the suite contains the same folder structure:

- `bin`: folder for the compiled binary.
- `obj`: folder with all binary objects.
- `src`: folder with source code.
- `makefile`: includes `../common.mk`, the shared build configuration.

Build configuration (compiler flags, library paths, build targets) is centralized in `common.mk` at the repository root. Each tool's `makefile` simply includes it.

### Quick build (system target)

If your libraries are installed in standard system locations, you can build all tools with:

<div class="code-example" markdown="1">
```bash
make system
```
</div>

The `system` target auto-detects your OS, architecture, and library locations. It searches standard installation paths across all major Linux distributions (Debian/Ubuntu, Fedora/RHEL, Arch, openSUSE, etc.) and macOS Homebrew. If `pkg-config` is installed, it is used for more precise HTSlib detection.

| Platform | Compiler | Detection method |
|----------|----------|-----------------|
| Linux (all distros) | g++ | Searches `/usr/lib`, `/usr/lib64`, `/usr/local/lib`, Debian multiarch paths |
| macOS ARM64 | clang++ | Searches `/opt/homebrew/lib`, plus pkg-config |

On x86_64, the phase module automatically enables AVX2/FMA SIMD instructions. On ARM64 platforms, SIMD is provided via NEON instructions through the SIMDe compatibility library.

### Custom build (other targets)

If the auto-detection doesn't find your libraries (e.g., they are installed in a non-standard location), you can override individual paths on the command line:

<div class="code-example" markdown="1">
```bash
make system SYS_BOOST_LIB=/path/to/boost/lib SYS_HTSLIB_LIBDIR=/path/to/htslib/lib
```
</div>

The overridable variables are:

- `SYS_HTSLIB_INC`: path to the HTSlib header files directory.
- `SYS_HTSLIB_LIBDIR`: path to the directory containing `libhts.a`.
- `SYS_BOOST_INC`: path to the Boost header files directory.
- `SYS_BOOST_LIB`: path to the directory containing the Boost static libraries.

Alternatively, you can use one of the other predefined targets (`desktop`, `docker`, etc.) or create your own by editing `common.mk`. Since all tools share the same `common.mk`, you only need to edit it once — the change applies to all tools automatically.

