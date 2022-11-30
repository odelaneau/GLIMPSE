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

#include "call_set_header.h"

call_set::call_set () : N(0), L(0), D(0), T(0), use_subset_samples(false), fploidy(2), n_haploid(0), n_diploid(0), ploidy(2), n_fields_in_group_files(0) {}

void call_set::initialize(std::vector < double > _maf_bins,double _T, int _D) {
	vrb.title("Initializing engine based on allele frequency bins:");
	use_subset_samples=false;
	T = _T;
	D = _D;
	bins = _maf_bins;
	L = bins.size() - 1;

	site2grp.clear();
	rsquared_str.clear();
	vrb.bullet("Validation GL Probs	>= " + stb.str(T));
	vrb.bullet("Validation Depth 	>= " + stb.str(D));
	vrb.bullet("#AF bins = " + stb.str(L));
	vrb.bullet("AF bins [" + stb.str(bins) + "]");
}

void call_set::initialize(double _T, int _D) {
	vrb.title("Initializing engine based on allele frequency bins:");
	use_subset_samples=false;
	T = _T;
	D = _D;

	site2grp.clear();
	rsquared_str.clear();
	vrb.bullet("Validation GL Probs	>= " + stb.str(T));
	vrb.bullet("Validation Depth 	>= " + stb.str(D));
	//vrb.bullet("#AF bins = " + stb.str(L));
	//vrb.bullet("AC bins [" + stb.str(ac_bins) + "]");
}




void call_set::initialize(std::string fgrps, double _T, int _D) {
	vrb.title("Initializing engine based on groups:");
	use_subset_samples=false;
	T = _T;
	D = _D;
	bins.clear();
	L = 0;

	//Read site2groups
	std::string buffer;
	input_file fd (fgrps);
	std::vector < std::string > tokens;
	std::map < std::string, std::pair < int, bool > > :: iterator itG;
	rsquared_str.clear();
	while (getline(fd, buffer, '\n'))
	{
		const int n_tokens = stb.split(buffer, tokens);
		if (n_tokens == 3 || n_tokens == 5)
		{
			if (n_fields_in_group_files == 0) n_fields_in_group_files = n_tokens;
			if (n_fields_in_group_files!=n_tokens) vrb.error("Group file is not well formatted. Requested 3 or 5 columns for all sites: CHR\tPOS\tGroupID or CHR\tPOS\tREF\tALT\tGroupID");

			const std::string uuid = (n_tokens == 5) ? tokens[0] + "_" + tokens[1] + "_" + tokens[2] + "_" + tokens[3] : tokens[0] + "_" + tokens[1];
			const std::string groupid = (n_tokens == 5) ? tokens[4] : tokens[2];

			itG = site2grp.find(uuid);
			if (itG == site2grp.end()) {
				//Search for group index
				int grp_idx = -1;
				for (int g  = 0 ; g < rsquared_str.size() && grp_idx < 0; g ++) if (rsquared_str[g] == groupid) grp_idx = g;
				if (grp_idx < 0) {
					grp_idx = rsquared_str.size();
					rsquared_str.push_back(groupid);
				}
				//Add new site of interest
				site2grp.insert(std::pair < std::string, std::pair < int, bool > > (uuid, std::pair < int, bool > (grp_idx, false)));
			}
		}
	}
	if (n_fields_in_group_files == 0) vrb.error("Group file is not well formatted. Requested 3 or 5 columns for all sites: CHR\tPOS\tGroupID or CHR\tPOS\tREF\tALT\tGroupID");
	vrb.bullet("Validation GL Probs	>= " + stb.str(T));
	vrb.bullet("Validation Depth 	>= " + stb.str(D));
	vrb.bullet("#sites  = " + stb.str(site2grp.size()));
	vrb.bullet("#groups = " + stb.str(rsquared_str.size()));
	fd.close();
}

void call_set::setTargets(std::string fsamples) {
	vrb.title("Reading the subset of samples to consider in the analysis");
	use_subset_samples=true;
	std::string buffer;
	input_file fd (fsamples);
	std::vector < std::string > tokens;
	std::pair<std::set<std::string>::iterator,bool> ret;

	while (getline(fd, buffer))
	{
		if (stb.split(buffer, tokens) < 1) vrb.error("Empty line found in sample file.");
		ret = subset_samples_set.insert(tokens[0]);
		if (ret.second) subset_samples.push_back(tokens[0]);

	}
	vrb.bullet("#samples  = " + stb.str(subset_samples.size()));
	fd.close();
	if (subset_samples.size()==0) vrb.error("No sample found in the samples file (--samples).");
}

void call_set::write_record(htsFile *out_fd, bcf_hdr_t * out_hdr, bcf_hdr_t * hdr_in, bcf1_t *line)
{
	//bcf_update_format(out_hdr, line, "GT", NULL, 0, BCF_HT_INT);
	if (bcf_write(out_fd, out_hdr, line) ) vrb.error("Failed to write the record output to file");
}

call_set::~call_set() {

}


