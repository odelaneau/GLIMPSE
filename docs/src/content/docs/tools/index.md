---
title: "Documentation"
description: "Command-line reference for the GLIMPSE2 tools."
---

GLIMPSE2 is a set of command-line tools. Each is installed as `GLIMPSE2_<name>`, for example `GLIMPSE2_phase`, and prints its options with `--help`.

| Tool | What it does |
|:-----|:-------------|
| [chunk](chunk/) | Defines the imputation chunks along a chromosome. |
| [split_reference](split_reference/) | Converts a reference panel to the binary format, one file per chunk. |
| [phase](phase/) | Imputes and phases samples from BAM/CRAM files or genotype likelihoods. |
| [ligate](ligate/) | Joins the imputed chunks into a chromosome-wide file. |
| [concordance](concordance/) | Computes genotyping and switch error rates against truth data. |
| [inspect](inspect/) | Summarises a binary reference panel. |

The first four are the imputation pipeline, in that order. `concordance` and `inspect` are not part of it. The [tutorial](../tutorial/getting-started/) runs the pipeline end to end.

In the option tables, settings marked **Expert**, **Recommended** or **Experimental** mean what they say. Everything else is a normal option, and the default is the value used when you leave it out.
