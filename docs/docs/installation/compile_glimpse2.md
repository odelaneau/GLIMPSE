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

The `system` target auto-detects your OS and architecture:

| Platform | Compiler | Library paths |
|----------|----------|---------------|
| Linux x86_64 | g++ | `/usr/lib/x86_64-linux-gnu` and `/usr/local/lib` |
| Linux aarch64 | g++ | `/usr/lib/aarch64-linux-gnu` and `/usr/local/lib` |
| macOS ARM64 | clang++ | `/opt/homebrew/lib` (Homebrew) |

On x86_64, the phase module automatically enables AVX2/FMA SIMD instructions. On ARM64 platforms, SIMD is provided via NEON instructions through the SIMDe compatibility library.

### Custom build (other targets)

If your libraries are installed in non-standard locations, you can use one of the other predefined targets (`desktop`, `docker`, etc.) or create your own by editing the makefile. The following variables need to be set:

- `HTSSRC`: path to the root of the HTSlib library, the prefix for HTSLIB_INC and HTSLIB_LIB paths.
- `HTSLIB_INC`: path to the HTSlib header files.
- `HTSLIB_LIB`: path to the static HTSlib library (file `libhts.a`).
- `BOOST_INC`: path to the Boost header files (usually `/usr/include`).
- `BOOST_LIB_IO`: path to the static Boost iostreams library (file `libboost_iostreams.a`).
- `BOOST_LIB_PO`: path to the static Boost `program_options` library (file `libboost_program_options.a`).
- `BOOST_LIB_SE`: path to the static Boost serialization library (file `libboost_serialization.a`).

If installed at the system level, static libraries (*.a files) can be located with this command:

<div class="code-example" markdown="1">
```bash
locate libboost_program_options.a libboost_iostreams.a libhts.a
```
</div>

Once all paths are correctly set up, proceed with the compilation using `make <target>`. The binary can be found in the `bin/` folder of each tool and will have a name similar to `GLIMPSE2_phase`. Since all tools share the same `common.mk`, you only need to edit it once — the change applies to all tools automatically.

