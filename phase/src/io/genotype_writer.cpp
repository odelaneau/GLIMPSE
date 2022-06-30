/*******************************************************************************
 * Copyright (C) 2020 Olivier Delaneau, University of Lausanne
 * Copyright (C) 2020 Simone Rubinacci, University of Lausanne
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

#include <io/genotype_writer.h>
#include "../../versions/versions.h"

#define OFILE_VCFU	0
#define OFILE_VCFC	1
#define OFILE_BCFC	2

genotype_writer::genotype_writer(const haplotype_set & _H, const genotype_set & _G, const variant_map & _V): H(_H), G(_G), V(_V) {
}

genotype_writer::~genotype_writer() {
}

void genotype_writer::writeGenotypes(const string fname, const int output_start, const int output_stop, const int n_main, const int n_threads) const {
	// Init
	tac.clock();
	string file_format = "w";
	unsigned int file_type = OFILE_VCFU;
	if (fname.size() > 6 && fname.substr(fname.size()-6) == "vcf.gz") { file_format = "wz"; file_type = OFILE_VCFC; }
	if (fname.size() > 3 && fname.substr(fname.size()-3) == "bcf") { file_format = "wb"; file_type = OFILE_BCFC; }
	htsFile * fp = hts_open(fname.c_str(),file_format.c_str());
	if (n_threads > 1) hts_set_threads(fp, n_threads);
	bcf_hdr_t * hdr = bcf_hdr_init("w");
	bcf1_t *rec = bcf_init1();

	// Create VCF header
	bcf_hdr_append(hdr, string("##fileDate="+tac.date()).c_str());
	bcf_hdr_append(hdr, string("##source=GLIMPSE_phase v" + string(PHASE_VERSION)).c_str());
	bcf_hdr_append(hdr, string("##contig=<ID="+ V.vec_pos[0]->chr + ">").c_str());
	bcf_hdr_append(hdr, "##INFO=<ID=RAF,Number=A,Type=Float,Description=\"ALT allele frequency in the reference panel\">");
	bcf_hdr_append(hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"ALT allele frequency computed from DS/GP field across target samples\">");
	if (H.fploidy == 2) bcf_hdr_append(hdr, "##INFO=<ID=INFO,Number=A,Type=Float,Description=\"IMPUTE info quality score for diploid samples\">");
	bcf_hdr_append(hdr, "##INFO=<ID=BUF,Number=A,Type=Integer,Description=\"Is it a variant site falling within buffer regions? (0=no/1=yes)\">");
	bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Unphased genotypes\">");
	bcf_hdr_append(hdr, "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Genotype dosage\">");
	bcf_hdr_append(hdr, "##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Genotype posteriors\">");
	bcf_hdr_append(hdr, "##FORMAT=<ID=HS,Number=1,Type=Integer,Description=\"Sampled haplotype pairs packed into intergers (max: 15 pairs, see NMAIN header line)\">");
	bcf_hdr_append(hdr, string("##NMAIN="+stb.str(n_main)).c_str());
	bcf_hdr_append(hdr, string("##FPLOIDY="+stb.str(H.fploidy)).c_str());

	//Add samples
	vector < int > ptr_gps = vector < int > (G.n_ind, 0);		// Pointers to iterate over sparse GPs
	for (int i = 0 ; i < G.n_ind ; i ++) bcf_hdr_add_sample(hdr, G.vecG[i]->name.c_str());
	bcf_hdr_add_sample(hdr, NULL);      // to update internal structures
	if (bcf_hdr_write(fp, hdr) )  vrb.error("Failed to write header to file [" + string(fname) + "]");

	//Add records
	int * genotypes = (int*)malloc(bcf_hdr_nsamples(hdr)*H.max_ploidy*sizeof(int));
	float * dosages = (float*)malloc(bcf_hdr_nsamples(hdr)*1*sizeof(float));
	float * posteriors = (float*)malloc(bcf_hdr_nsamples(hdr)*(H.max_ploidy+1)*sizeof(float));
	int * haplotypes = (int*)malloc(bcf_hdr_nsamples(hdr)*sizeof(int));

	const std::size_t nsites = V.vec_pos.size();
	float prog_step = 1.0/nsites;
	float prog_bar = 0.0;

	for (int l = 0 ; l < nsites ; l ++) {
		// Clear current VCF record
		bcf_clear1(rec);

		// Update variant informations
		rec->rid = bcf_hdr_name2id(hdr, V.vec_pos[l]->chr.c_str());
		rec->pos = V.vec_pos[l]->bp - 1;
		bcf_update_id(hdr, rec, V.vec_pos[l]->id.c_str());
		string alleles = V.vec_pos[l]->ref + "," + V.vec_pos[l]->alt;
		bcf_update_alleles_str(hdr, rec, alleles.c_str());

		// Store individual data
		float ds_sum=0.0, ds2_sum=0.0, ds4_sum=0.0;
		std::fill(&genotypes[0], &genotypes[0] + H.max_ploidy*G.vecG.size(), bcf_gt_unphased(false));
		for (int i = 0 ; i < G.n_ind ; i++)
		{
			// Initialialize output data in case of Ref/Ref genotype
			int32_t hs = 0;
			float ds = 0.0f, gp0 = 1.0f, gp1 = 0.0f, gp2 = 0.0f;

			// Genotype is NOT 100% certain Ref/Ref
			if ((ptr_gps[i]<G.vecG[i]->stored_data.size()) && (G.vecG[i]->stored_data[ptr_gps[i]].idx==l))
			{
				genotypes[H.max_ploidy*i+0] = bcf_gt_unphased(G.vecG[i]->H0[l]);
				gp0 = G.vecG[i]->stored_data[ptr_gps[i]].gp0;
				gp1 = G.vecG[i]->stored_data[ptr_gps[i]].gp1;
				hs = G.vecG[i]->stored_data[ptr_gps[i]].hs;
				ds += gp1;
				if (G.vecG[i]->ploidy > 1)
				{
					genotypes[H.max_ploidy*i+1] = bcf_gt_unphased(G.vecG[i]->H1[l]);
					gp2 = G.vecG[i]->stored_data[ptr_gps[i]].getGp2();
					ds += 2.0f * gp2;
					alg.clamp(ds, 0.0f, 2.0f);
				}
				ptr_gps[i] ++;
			}

			// Store DS + GP rounded 4 decimals
			dosages[i] = roundf(ds * 1000.0) / 1000.0;
			posteriors[(H.max_ploidy+1)*i+0] = roundf(gp0 * 1000.0) / 1000.0;
			posteriors[(H.max_ploidy+1)*i+1] = roundf(gp1 * 1000.0) / 1000.0;
			if (H.max_ploidy>1)
			{
				if (G.vecG[i]->ploidy > 1) posteriors[(H.max_ploidy+1)*i+2] = roundf(gp2 * 1000.0) / 1000.0;
				else
				{
					genotypes[H.max_ploidy*i+1] = bcf_int32_vector_end;
					bcf_float_set(&posteriors[(H.max_ploidy+1)*i+2], bcf_float_vector_end);
				}
			}
			haplotypes[i] = hs;

			// Compute INFO/INFO statistics
			ds_sum += ds;
			if (H.fploidy == 2)
			{
				ds2_sum += ds * ds;
				ds4_sum += gp1 + 4.0*gp2;
			}
		}

		// Update INFO fields
		int buffer = ((V.vec_pos[l]->bp < output_start) || (V.vec_pos[l]->bp > output_stop));
		float freq_alt_refp = V.vec_pos[l]->calt * 1.0f / (V.vec_pos[l]->calt + V.vec_pos[l]->cref);
		float freq_alt_main = ds_sum / H.n_main_haps;
		float infoscore = (H.fploidy == 2 && freq_alt_main>0.0f && freq_alt_main<1.0f) ? (float)(1.0f - (ds4_sum - ds2_sum) / (H.n_main_haps * freq_alt_main * (1.0f - freq_alt_main))) : 1.0f;
		infoscore = roundf(std::max(infoscore,0.0f) * 1000.0) / 1000.0;
		bcf_update_info_float(hdr, rec, "RAF", &freq_alt_refp, 1);
		bcf_update_info_float(hdr, rec, "AF", &freq_alt_main, 1);
		if (H.fploidy == 2) bcf_update_info_float(hdr, rec, "INFO", &infoscore, 1);
		bcf_update_info_int32(hdr, rec, "BUF", &buffer, 1);

		// Update FORMAT fields
		bcf_update_genotypes(hdr, rec, genotypes, bcf_hdr_nsamples(hdr)*H.max_ploidy);
		bcf_update_format_float(hdr, rec, "DS", dosages, bcf_hdr_nsamples(hdr)*1);
		bcf_update_format_float(hdr, rec, "GP", posteriors, bcf_hdr_nsamples(hdr)*(H.max_ploidy+1));
		bcf_update_format_int32(hdr, rec, "HS", haplotypes, bcf_hdr_nsamples(hdr)*1);

		//Write record
		if (bcf_write1(fp, hdr, rec) ) vrb.error("Failed to write the record [" + V.vec_pos[l]->chr + ":" + std::to_string(V.vec_pos[l]->bp) + "_" + V.vec_pos[l]->ref + "_" + V.vec_pos[l]->alt +"] to file [" + string(fname) + "]");
		prog_bar+=prog_step;
		vrb.progress("  * VCF writing", prog_bar);
	}
	free(genotypes);
	free(dosages);
	free(posteriors);
	free(haplotypes);
	bcf_destroy1(rec);
	bcf_hdr_destroy(hdr);
	if (hts_close(fp)) vrb.error("Non zero status when closing VCF/BCF file descriptor");
	switch (file_type) {
	case OFILE_VCFU: vrb.bullet("VCF writing [Uncompressed / N=" + stb.str(G.n_ind) + " / L=" + stb.str(V.size()) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_VCFC: vrb.bullet("VCF writing [Compressed / N=" + stb.str(G.n_ind) + " / L=" + stb.str(V.size()) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_BCFC: vrb.bullet("BCF writing [Compressed / N=" + stb.str(G.n_ind) + " / L=" + stb.str(V.size()) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	}
}
