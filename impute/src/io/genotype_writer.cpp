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

#define OFILE_VCFU	0
#define OFILE_VCFC	1
#define OFILE_BCFC	2

auto compLB =  [](auto const& m1, auto const& m2) { return m1->bp < m2; };

genotype_writer::genotype_writer(haplotype_set & _H, genotype_set & _G, variant_map & _V, std::vector < probability_set * >& _P, gmap_reader & _readerGM, std::string _region, float _maf_common): H(_H), G(_G), V(_V), P(_P), readerGM(_readerGM), region(_region), maf_common(_maf_common)
{
	//COUNTS
	n_variants = H.n_site;
	n_main_samples = H.n_targ/2;
	n_ref_samples = H.n_ref/2;
	e2sum = 0.0;
	fsum = 0.0;
	esum=0.0;

	genotypes = (int*)malloc(n_main_samples*2*sizeof(int));
	dosages = (float*)malloc(n_main_samples*1*sizeof(float));
	posteriors = (float*)malloc(n_main_samples*3*sizeof(float));
}

genotype_writer::~genotype_writer()
{
	free(genotypes);
	free(dosages);
	free(posteriors);
}

void genotype_writer::writeGenotypesAndImpute(std::string funphased, std::string freference, std::string fname, int start, int stop) {
	tac.clock();
	genotype_stream input_stream(region, maf_common, n_variants, n_main_samples, n_ref_samples);
	input_stream.openStream(funphased, freference);

	tac.clock();
	std::string file_format = "w";
	unsigned int file_type = OFILE_VCFU;
	if (fname.size() > 6 && fname.substr(fname.size()-6) == "vcf.gz") { file_format = "wz"; file_type = OFILE_VCFC; }
	if (fname.size() > 3 && fname.substr(fname.size()-3) == "bcf") { file_format = "wb"; file_type = OFILE_BCFC; }
	htsFile * fp = hts_open(fname.c_str(),file_format.c_str());
	bcf_hdr_t * hdr = bcf_hdr_init("w");
	bcf1_t *rec = bcf_init1();

	// Create VCF header
	bcf_hdr_append(hdr, std::string("##fileDate="+tac.date()).c_str());
	bcf_hdr_append(hdr, "##source=LCC impute");
	bcf_hdr_append(hdr, std::string("##contig=<ID="+ V.vec_pos[0]->chr + ">").c_str());
	bcf_hdr_append(hdr, "##INFO=<ID=INFO,Number=1,Type=Float,Description=\"IMPUTE2 info score, AF from expected genotype dosage\">");
	bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotypes\">");
	bcf_hdr_append(hdr, "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Genotype dosage\">");
	bcf_hdr_append(hdr, "##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Genotype posteriors\">");

	for (int i = 0 ; i < G.n_ind ; i ++) bcf_hdr_add_sample(hdr, G.vecG[i]->name.c_str());
	bcf_hdr_add_sample(hdr, NULL);      // to update internal structures
	static_cast<void>(bcf_hdr_write(fp, hdr));
	int32_t tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
	const double counter_limit = 1.0/V.size();

	variant& curr_variant = *V.vec_pos[0];
	auto it_minPosW = std::lower_bound(V.vec_pos.begin(), V.vec_pos.end(), start, compLB);
	if (it_minPosW == V.vec_pos.end())	vrb.error("Error finding first unbuffered position.");

	input_stream.i_next_common_variant = std::distance(V.vec_pos.begin(),it_minPosW);
	input_stream.i_common_variant =  input_stream.i_next_common_variant-1;

	int current_position;
	double est_af;
	float info, count_alt_ref;
	const double n_main_haps_d = 2*(double)n_main_samples;

	while (input_stream.readMarker())
	{
		if (input_stream.is_common_variant)
			curr_variant = *V.vec_pos[input_stream.i_common_variant];//COPY, maybe a reference is enough
		else curr_variant = input_stream.curr_variant;

		current_position = curr_variant.bp;
		if (current_position >= start && current_position < stop)
		{
			bcf_clear1(rec);

			esum=0.0;e2sum=0.0;fsum=0.0;
			count_alt_ref = curr_variant.calt;

			if (input_stream.is_common_variant) phase_marker(input_stream.i_common_variant, input_stream);
			else
					if (input_stream.i_common_variant<0) impute_marker_border(0, input_stream, true);
					else impute_marker(input_stream.i_common_variant, get_linear_interp_weight(input_stream.i_common_variant, current_position), input_stream);

			//est_af = (esum + (double)count_alt_ref) / (2 * (double)(n_main_samples+n_ref_samples)); //AF in main + ref
			est_af = esum / n_main_haps_d; //AF in main - Original INFO score
			info = (est_af>0.0 && est_af<1.0) ? (float)(1.0 - (fsum - e2sum) / (n_main_haps_d * est_af * (1.0 - est_af))) : 1; //info can be negative

			rec->rid = bcf_hdr_name2id(hdr, curr_variant.chr.c_str());
			rec->pos = current_position - 1;
			bcf_update_id(hdr, rec, curr_variant.id.c_str());
			bcf_update_alleles_str(hdr, rec, std::string(curr_variant.ref + "," + curr_variant.alt).c_str());

			bcf_update_filter(hdr,rec,&tmpi,1);
		    bcf_update_info_float(hdr, rec, "INFO", &info, 1);
			bcf_update_genotypes(hdr, rec, genotypes, bcf_hdr_nsamples(hdr)*2);
			bcf_update_format_float(hdr, rec, "DS", dosages, bcf_hdr_nsamples(hdr)*1);
			bcf_update_format_float(hdr, rec, "GP", posteriors, bcf_hdr_nsamples(hdr)*3);
			static_cast<void>(bcf_write1(fp, hdr, rec));
		}
		vrb.progress("  * VCF writing", (input_stream.i_common_variant+1)*counter_limit);
	}

	bcf_destroy1(rec);
	bcf_hdr_destroy(hdr);
	if (hts_close(fp)) vrb.error("Non zero status when closing VCF/BCF file descriptor");
	switch (file_type) {
	case OFILE_VCFU: vrb.bullet("VCF writing [Uncompressed / N=" + stb.str(G.n_ind) + " / L=" + stb.str(V.size()) + " / Ltot=" + stb.str(input_stream.i_variant) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_VCFC: vrb.bullet("VCF writing [Compressed / N=" + stb.str(G.n_ind) + " / L=" + stb.str(V.size()) + " / Ltot=" + stb.str(input_stream.i_variant) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_BCFC: vrb.bullet("BCF writing [Compressed / N=" + stb.str(G.n_ind) + " / L=" + stb.str(V.size()) + " / Ltot=" + stb.str(input_stream.i_variant) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	}
}

void genotype_writer::phase_marker(const int l, const genotype_stream& input_stream)
{
	bool a0, a1;
	float gp0, gp1, gp2, ds;

	for (int i = 0 ; i < G.n_ind ; i++) {
		a0 = G.vecG[i]->H0[l];
		a1 = G.vecG[i]->H1[l];
		genotypes[2*i+0] = bcf_gt_phased(a0);
		genotypes[2*i+1] = bcf_gt_phased(a1);
		gp0 = G.vecG[i]->GP[2*l+0];
		gp1 = G.vecG[i]->GP[2*l+1];
		gp2 = abs(1.0 - gp0 - gp1);
		ds = gp1 + 2 * gp2;
		dosages[i] = roundf(ds * 1000.0) / 1000.0;
		posteriors[3*i+0] = roundf(gp0 * 1000.0) / 1000.0;;
		posteriors[3*i+1] = roundf(gp1 * 1000.0) / 1000.0;;
		posteriors[3*i+2] = roundf(gp2 * 1000.0) / 1000.0;;

		esum+= posteriors[3*i + 1] + 2.0*posteriors[3*i + 2];
		e2sum+= (posteriors[3*i + 1] + 2.0*posteriors[3*i + 2]) * (posteriors[3*i + 1] + 2.0*posteriors[3*i + 2]);
		fsum+= posteriors[3*i + 1] + 4.0*posteriors[3*i + 2];
	}
}
void genotype_writer::impute_marker(const int l, const float w,  const genotype_stream& input_stream)
{
	const float w1M = std::max(1.0 - w,0.0);
	float sum = 0.0f;
	float post0, post1,post2;
	std::array<float,2> p1hL {0.0f,0.0f};
	std::array<float,2> p1hR {0.0f,0.0f};
	float p1h0Comb, p1h1Comb;
	int max_el;
	bool a0,a1;

	for (int i = 0 ; i < G.n_ind ; i++)
	{
		p1hL[0] = 0.0f;
		p1hL[1] = 0.0f;
		p1hR[0] = 0.0f;
		p1hR[1] = 0.0f;

		for (int h = 0; h < 2; ++ h)
		{
			//sum 1s for h
			const std::vector<int>& pS = P[i]->reference_haps[h][l];
			const std::vector<float>& pL = P[i]->prob_stateL[h][l];
			const std::vector<float>& pR = P[i]->prob_stateR[h][l];
			const size_t num_states = pS.size();

			for (size_t j =0; j<num_states; ++j)
				if (input_stream.ref_haps.get(pS[j]))
				{
					p1hL[h]+=pL[j];
					p1hR[h]+=pR[j];
				}
		}
		p1h0Comb = p1hL[0] * w + p1hR[0] * w1M;
		p1h1Comb = p1hL[1] * w + p1hR[1] * w1M;

		p1h0Comb = 0.9999 * p1h0Comb + 0.0001 * abs(1.0f - p1h0Comb);
		p1h1Comb = 0.9999 * p1h1Comb + 0.0001 * abs(1.0f - p1h1Comb);

		post0 = abs(1.0f - p1h0Comb) * abs(1.0f - p1h1Comb);
		post2 = p1h0Comb * p1h1Comb;
		post1 = abs(1.0f - post0 - post2);

		//curr_GL are assured to be .333 if flat likelihood or values with sum > 0.0
		///otherwise division problem!
		post0 *= input_stream.curr_GL[3*i + 0];
		post1 *= input_stream.curr_GL[3*i + 1];
		post2 *= input_stream.curr_GL[3*i + 2];

		sum = post0 + post1 + post2;

		post0/=sum;
		post1/=sum;
		post2/=sum;

		max_el = max3gt(post0,post1,post2);
		a0 = max_el > 0;
		a1 = max_el > 1;
		if (max_el ==1 && p1h0Comb < p1h1Comb)
		{
			a0 = false;
			a1 = true;
		}
		genotypes[2*i+0] = bcf_gt_phased(a0);
		genotypes[2*i+1] = bcf_gt_phased(a1);

		posteriors[3*i + 0] = roundf(post0 * 1000.0) / 1000.0;
		posteriors[3*i + 1] = roundf(post1 * 1000.0) / 1000.0;
		posteriors[3*i + 2] = roundf(post2 * 1000.0) / 1000.0;

		esum+= posteriors[3*i + 1] + 2.0*posteriors[3*i + 2];
		e2sum+= (posteriors[3*i + 1] + 2.0*posteriors[3*i + 2]) * (posteriors[3*i + 1] + 2.0*posteriors[3*i + 2]);
		fsum+= posteriors[3*i + 1] + 4.0*posteriors[3*i + 2];

		dosages[i] = roundf((posteriors[3*i + 1] + 2 * posteriors[3*i + 2]) * 1000.0) / 1000.0;
	}
}

void genotype_writer::impute_marker_border(const int l, const genotype_stream& input_stream, const bool left)
{
	float sum = 0.0f;
	float post0, post1,post2;
	std::array<float,2> p1h {0.0,0.0};
	float p1h0Comb, p1h1Comb;
	int max_el;
	bool a0,a1;

	for (int i = 0 ; i < G.n_ind ; i++)
	{
		p1h[0] = 0.0f;
		p1h[1] = 0.0f;

		for (size_t h = 0; h < 2; ++ h)
		{
			//sum 1s for h
			const std::vector<int>& pS = P[i]->reference_haps[h][l];
			const std::vector<float>& pLR = left ? P[i]->prob_stateL[h][l] : P[i]->prob_stateR[h][l];
			const size_t num_states = pS.size();

			for (size_t j =0; j<num_states; ++j) if (input_stream.ref_haps.get(pS[j])) p1h[h]+=pLR[j];
		}

		p1h0Comb = 0.9999 * p1h[0] + 0.0001 * abs(1.0f - p1h[0]);
		p1h1Comb = 0.9999 * p1h[1] + 0.0001 * abs(1.0f - p1h[1]);

		post0 = abs(1.0f - p1h0Comb) * (1.0f - p1h1Comb);
		post2 = p1h0Comb * p1h1Comb;
		post1 = abs(1.0f - post0 - post2);

		//curr_GL are assured to be .333 if flat likelihood or values with sum > 0.0
		///otherwise division problem!
		post0 *= input_stream.curr_GL[3*i + 0];
		post1 *= input_stream.curr_GL[3*i + 1];
		post2 *= input_stream.curr_GL[3*i + 2];

		sum = post0 + post1 + post2;

		post0/=sum;
		post1/=sum;
		post2/=sum;

		max_el = max3gt(post0,post1,post2);

		a0 = max_el > 0;
		a1 = max_el > 1;
		if (max_el ==1 && p1h[0] < p1h[1])
		{
			a0 = false;
			a1 = true;
		}
		genotypes[2*i+0] = bcf_gt_phased(a0);
		genotypes[2*i+1] = bcf_gt_phased(a1);

		posteriors[3*i + 0] = roundf(post0 * 1000.0) / 1000.0;
		posteriors[3*i + 1] = roundf(post1 * 1000.0) / 1000.0;
		posteriors[3*i + 2] = roundf(post2 * 1000.0) / 1000.0;

		esum+= posteriors[3*i + 1] + 2.0*posteriors[3*i + 2];
		e2sum+= (posteriors[3*i + 1] + 2.0*posteriors[3*i + 2]) * (posteriors[3*i + 1] + 2.0*posteriors[3*i + 2]);
		fsum+= posteriors[3*i + 1] + 4.0*posteriors[3*i + 2];

		dosages[i] = roundf((posteriors[3*i + 1] + 2 * posteriors[3*i + 2]) * 1000.0) / 1000.0;
	}
}

float genotype_writer::get_linear_interp_weight(const int l, const int current_position)
{
	float genpos_a = V.vec_pos[l]->cm;
	float genpos_b = V.vec_pos[l+1]->cm;
	float genpos_x = V.getCentiMorganPos(readerGM,current_position)-V.baseline;


	float genpos_ba = abs(genpos_b - genpos_a);
	float genpos_bx = abs(genpos_b - genpos_x);

	//need to be a bit careful about genpos to avoid nan values
	if (genpos_ba < 1e-7) genpos_ba = 1e-7;
	if (genpos_bx > genpos_ba) genpos_bx = genpos_ba;

	return (genpos_bx / genpos_ba);
}

int genotype_writer::max3gt(const float& a,const float& b, const float& c) const
{
	return a < b? (b < c? 2: 1): (a < c? 2: 0);
}
