////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2018 Olivier Delaneau, University of Lausanne
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////
#include <io/genotype_writer.h>
#include <version/version.h>


#define OFILE_VCFU	0
#define OFILE_VCFC	1
#define OFILE_BCFC	2

genotype_writer::genotype_writer(haplotype_set & _H, genotype_set & _G, variant_map & _V): H(_H), G(_G), V(_V) {
}

genotype_writer::~genotype_writer() {
}

void genotype_writer::writeGenotypes(string fname, int output_start, int output_stop, int n_main) {
	// Init
	tac.clock();
	string file_format = "w";
	unsigned int file_type = OFILE_VCFU;
	if (fname.size() > 6 && fname.substr(fname.size()-6) == "vcf.gz") { file_format = "wz"; file_type = OFILE_VCFC; }
	if (fname.size() > 3 && fname.substr(fname.size()-3) == "bcf") { file_format = "wb"; file_type = OFILE_BCFC; }
	htsFile * fp = hts_open(fname.c_str(),file_format.c_str());
	bcf_hdr_t * hdr = bcf_hdr_init("w");
	bcf1_t *rec = bcf_init1();

	// Create VCF header
	bcf_hdr_append(hdr, string("##fileDate="+tac.date()).c_str());
	bcf_hdr_append(hdr, string("##source=LCC_phase v" + string(VERSION)).c_str());
	bcf_hdr_append(hdr, string("##contig=<ID="+ V.vec_pos[0]->chr + ">").c_str());
	bcf_hdr_append(hdr, "##INFO=<ID=AFref,Number=A,Type=Float,Description=\"Allele Frequency\">");
	bcf_hdr_append(hdr, "##INFO=<ID=AFmain,Number=A,Type=Float,Description=\"Allele Frequency\">");
	bcf_hdr_append(hdr, "##INFO=<ID=AFfull,Number=A,Type=Float,Description=\"Allele Frequency\">");
	bcf_hdr_append(hdr, "##INFO=<ID=CM,Number=A,Type=Float,Description=\"Interpolated cM position\">");
	bcf_hdr_append(hdr, "##INFO=<ID=BUF,Number=A,Type=Integer,Description=\"Is it a buffer specific variant site? (0=no/1=yes)\">");
	bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotypes\">");
	bcf_hdr_append(hdr, "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Genotype dosage\">");
	bcf_hdr_append(hdr, "##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Genotype posteriors\">");
	bcf_hdr_append(hdr, "##FORMAT=<ID=HS,Number=1,Type=Integer,Description=\"Sampled haplotype pairs packed into intergers (max: 16 pairs, see NMAIN header line)\">");
	bcf_hdr_append(hdr, string("##NMAIN="+stb.str(n_main)).c_str());

	//Add samples
	for (int i = 0 ; i < G.n_ind ; i ++) bcf_hdr_add_sample(hdr, G.vecG[i]->name.c_str());
	bcf_hdr_add_sample(hdr, NULL);      // to update internal structures
	bcf_hdr_write(fp, hdr);

	//Add records
	int * genotypes = (int*)malloc(bcf_hdr_nsamples(hdr)*2*sizeof(int));
	float * dosages = (float*)malloc(bcf_hdr_nsamples(hdr)*1*sizeof(float));
	float * posteriors = (float*)malloc(bcf_hdr_nsamples(hdr)*3*sizeof(float));
	int * haplotypes = (int*)malloc(bcf_hdr_nsamples(hdr)*sizeof(int));

	for (int l = 0 ; l < V.size() ; l ++) {
		int current_position = V.vec_pos[l]->bp;
		int buffer = ((current_position < output_start) || (current_position >= output_stop));
		bcf_clear1(rec);
		rec->rid = bcf_hdr_name2id(hdr, V.vec_pos[l]->chr.c_str());
		rec->pos = current_position - 1;
		bcf_update_id(hdr, rec, V.vec_pos[l]->id.c_str());
		string alleles = V.vec_pos[l]->ref + "," + V.vec_pos[l]->alt;
		bcf_update_alleles_str(hdr, rec, alleles.c_str());
		float count_alt = 0, count_alt_ref = V.vec_pos[l]->calt;
		for (int i = 0 ; i < G.n_ind ; i++) {
			bool a0 = G.vecG[i]->H0[l];
			bool a1 = G.vecG[i]->H1[l];
			genotypes[2*i+0] = bcf_gt_phased(a0);
			genotypes[2*i+1] = bcf_gt_phased(a1);
			float gp0 = G.vecG[i]->GP[2*l+0];
			float gp1 = G.vecG[i]->GP[2*l+1];
			float gp2 = abs(1.0 - gp0 - gp1);
			float ds = gp1 + 2 * gp2;
			count_alt += ds;
			gp0 = roundf(gp0 * 1000.0) / 1000.0;
			gp1 = roundf(gp1 * 1000.0) / 1000.0;
			gp2 = roundf(gp2 * 1000.0) / 1000.0;
			dosages[i] = roundf(ds * 1000.0) / 1000.0;
			posteriors[3*i+0] = gp0;
			posteriors[3*i+1] = gp1;
			posteriors[3*i+2] = gp2;
			haplotypes[i] = G.vecG[i]->HAP[l];
		}
		//INFOs
		float freq_alt_refp = count_alt_ref / (V.vec_pos[l]->calt + V.vec_pos[l]->cref);
		float freq_alt_main = count_alt / (2 * G.n_ind);
		float freq_alt_full = (count_alt + count_alt_ref) / (2 * G.n_ind + V.vec_pos[l]->calt + V.vec_pos[l]->cref);
		bcf_update_info_float(hdr, rec, "AFmain", &freq_alt_main, 1);
		bcf_update_info_float(hdr, rec, "AFref", &freq_alt_refp, 1);
		bcf_update_info_float(hdr, rec, "AFfull", &freq_alt_full, 1);
		bcf_update_info_int32(hdr, rec, "BUF", &buffer, 1);

		//FORMATs
		float val = (float)V.vec_pos[l]->cm;
		bcf_update_info_float(hdr, rec, "CM", &val, 1);
		bcf_update_genotypes(hdr, rec, genotypes, bcf_hdr_nsamples(hdr)*2);
		bcf_update_format_float(hdr, rec, "DS", dosages, bcf_hdr_nsamples(hdr)*1);
		bcf_update_format_float(hdr, rec, "GP", posteriors, bcf_hdr_nsamples(hdr)*3);
		bcf_update_format_int32(hdr, rec, "HS", haplotypes, bcf_hdr_nsamples(hdr)*1);

		//Write record
		bcf_write1(fp, hdr, rec);
		vrb.progress("  * VCF writing", (l+1)*1.0/V.size());
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
