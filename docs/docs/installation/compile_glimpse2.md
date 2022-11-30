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
git clone https://github.com/odelaneau/glimpse.git
```
</div>

Navigate to the downloaded folder using `cd glimpse`.

You'll find there a folder containing all the software packages are other utility folders:

- **chunk**: program to define chunks where to run phasing and imputation.
- **common**: basic source files used by different tools.
- **concordance**: program to verify the accuracy of low-coverage imputation against high-coverage genomes
- **docker**: all scripts needed to build a docker file comprising all binaries
- **docs**: documentation in html
- **ligate**: ligate multiple imputed BCF/VCF files into a single chromosome-length file
- **maps**: genetics maps in b37 and b38
- **phase**: program to impute and phase low-coverage data.
- **split_reference**: prorgram to create GLIMPSE2's reference file format, used by GLIMPSE2_phase
- **tutorial**: test datasets and scripts
- **versions**: versioning

Each program in the suite contains the same folder structure:

- `bin`: folder for the compiled binary.
- `obj`: folder with all binary objects.
- `src`: folder with source code.
- `makefile`: Makefile to compile the program.

In order to compile a specific tool, for example _GLIMPSE2_phase_, go in directory of the software (cd `phase`) and edit the specific makefile at lines so that the following variables are correctly set up (look at the paths already there for an example):

- `HTSSRC`: path to the root of the HTSlib library, the prefix for HTSLIB_INC and HTSLIB_LIB paths.
- `HTSLIB_INC`: path to the HTSlib header files.
- `HTSLIB_LIB`: path to the static HTSlib library (file `libhts.a`).
- `BOOST_INC`: path to the BOOST header files (usually `/usr/include`). 
- `BOOST_LIB_IO`: path to the static BOOST iostreams library (file `libboost_iostreams.a`). 
- `BOOST_LIB_PO`: path to the static BOOST `program_options` library (file `libboost_program_options.a`). 
- `BOOST_LIB_SE`: path to the static BOOST serialization library (file `libboost_serialization.a`).

If installed at the system level, static libraries (*.a files) can be located with this command:

<div class="code-example" markdown="1">
```bash
locate libboost_program_options.a libboost_iostreams.a libhts.a
```
</div>

Once all paths correctly set up, proceed with the compilation using `make`. The binary can be found in the `bin/` folder of each tool and will have a name similar to `GLIMPSE2_phase`. You will need to copy the modified makefile in each tool (folder) of GLIMPSE2.

