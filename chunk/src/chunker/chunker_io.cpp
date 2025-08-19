/*******************************************************************************
 * Copyright (C) 2022-2023 Simone Rubinacci
 * Copyright (C) 2022-2023 Olivier Delaneau
 *
 * MIT Licence
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include "chunker_header.h"

void chunker::readData(std::string fmain, std::string region, long int nthreads) {
    // Start timing and log start of reading process
    tac.clock();
    vrb.title("Reading input files");
    vrb.progress("  * Reading              : [" + fmain + "]");

    // Initialize synchronous BCF/VCF reader (htslib)
    bcf_srs_t * sr =  bcf_sr_init();
    sr->collapse = COLLAPSE_NONE; // Keep multi-allelic variants separate
    sr->require_index = 1;        // Require an index file (.tbi/.csi)

    // Enable multi-threaded reading if requested
    if (nthreads > 1) bcf_sr_set_threads(sr, nthreads);

    // Restrict reading to given region (chr or chr:start-end)
    if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1)
        vrb.error("Impossible to jump to region [" + region + "] in [" + fmain + "]");

    // Add the main input file to the reader
    if(!(bcf_sr_add_reader(sr, fmain.c_str()))) {
        // Report error depending on whether it's an open or index failure
        if (sr->errnum != idx_load_failed)
            vrb.error("Failed to open file: " + fmain);
        else
            vrb.error("Failed to load index of the file: " + fmain);
    }

    // If whole chromosome mode: get contig length from header/index
    if (whole_chr) {
        const char **seq;
        int i, nseq;
        tbx_t *tbx = NULL;
        hts_idx_t *idx = NULL;

        // Decide if input is VCF or BCF and load correct type of index
        if (hts_get_format(sr->readers[0].file)->format == vcf) {
            tbx = tbx_index_load(fmain.c_str());
            if (!tbx) vrb.error("Could not load index for VCF: " + fmain);
        }
        else if (hts_get_format(sr->readers[0].file)->format == bcf) {
            idx = bcf_index_load(fmain.c_str());
            if (!idx) vrb.error("Could not load index for BCF: " + fmain);
        }
        else vrb.error("Could not detect the file type as VCF or BCF: " + fmain);

        // Retrieve sequence (chromosome) names
        seq = tbx ? tbx_seqnames(tbx, &nseq)
                  : bcf_index_seqnames(idx, sr->readers[0].header, &nseq);

        // Find the matching sequence for the region and get its length
        for (i = 0; i < nseq; i++) {
            if (seq[i] != region) continue; // Skip non-matching contigs
            bcf_hrec_t *hrec = sr ? bcf_hdr_get_hrec(sr->readers[0].header, BCF_HL_CTG, "ID", seq[i], NULL) : NULL;
            int hkey = hrec ? bcf_hrec_find_key(hrec, "length") : -1;

            // Set contig length if available, otherwise use fallback value
            if (hkey >= 0) this->contig_len = atol(hrec->vals[hkey]);
            if (this->contig_len <= 0) this->contig_len = 1248956422; // fallback
            break;
        }

        // Clean up sequence/index memory
        free(seq);
        if (tbx) tbx_destroy(tbx);
        if (idx) hts_idx_destroy(idx);
    }

    // Counters and temporary storage for INFO fields
    int n_variants = 0;             // Total number of biallelic variants
    int n_comm_variants_cnt = 0;    // Number of common variants
    int *vAC = NULL, *vAN = NULL;   // Buffers for AC and AN
    int nAC = 0, nAN = 0;           // Sizes of AC/AN buffers
    int rAC = 0, rAN = 0;           // Return values from htslib calls

    bcf1_t * line_main; // Pointer to current variant record

    // Loop through all records in the specified region
    while (bcf_sr_next_line(sr)) {
        line_main = bcf_sr_get_line(sr, 0);

        // Keep only biallelic sites (n_allele == 2)
        if (line_main->n_allele == 2) {
            n_variants++;

            // Store chromosome name for the first common variant
            if (n_comm_variants_cnt < 1)
                chrID = bcf_hdr_id2name(sr->readers[0].header, line_main->rid);

            // Read INFO/AC and INFO/AN values
            rAC = bcf_get_info_int32(sr->readers[0].header, line_main, "AC", &vAC, &nAC);
            rAN = bcf_get_info_int32(sr->readers[0].header, line_main, "AN", &vAN, &nAN);

            // AC and AN must be scalar integers
            if ((nAC != 1) || (nAN != 1))
                vrb.error("VCF/BCF needs AC/AN INFO fields to be present");

            // Compute minor allele frequency
            double MAF = std::min(
                ((vAN[0] - vAC[0]) * 1.0f) / (vAN[0]),
                (vAC[0] * 1.0f) / (vAN[0])
            );

            // Classify as common if MAF >= sparse_maf
            const bool is_common = MAF >= this->sparse_maf;

            // Store position (convert 0-based VCF to 1-based)
            this->positions_all_mb.push_back(line_main->pos + 1);

            // Map position → index in positions_all_mb
            this->map_positions_all.insert(
                std::pair<int, int>(line_main->pos + 1, positions_all_mb.size() - 1)
            );

            // Map index in "all variants" → index in "common variants" (or current count)
            this->all2common.push_back(n_comm_variants_cnt);

            // If common, store in common variant structures
            if (is_common) {
                this->positions_common_mb.push_back(line_main->pos + 1);
                this->common2all.push_back(positions_all_mb.size() - 1);
                n_comm_variants_cnt++;
            }
        }
    }

    // Clean up reader and temporary buffers
    bcf_sr_destroy(sr);
    if (vAC) free(vAC);
    if (vAN) free(vAN);

    // Log final stats
    vrb.bullet("Input file read      : [" + fmain + "] (" + stb.str(tac.rel_time() * 1.0 / 1000, 2) + "s)");
    vrb.bullet("#variants            : [" + stb.str(n_variants) + " total | "
               + stb.str(n_variants - n_comm_variants_cnt) + " rare | "
               + stb.str(n_comm_variants_cnt) + " common]");
}
