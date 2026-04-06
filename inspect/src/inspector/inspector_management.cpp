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

#include <inspector/inspector_header.h>
#include <boost/archive/binary_iarchive.hpp>
#include <htslib/vcf.h>
#include <sys/stat.h>

inspector::inspector() {
}

inspector::~inspector() {
}

void inspector::read_binary_panel() {
	std::string filename = options["input"].as < std::string > ();
	vrb.wait("  * Reading binary reference panel");
	tac.clock();

	std::ifstream ifs(filename, std::ios::binary | std::ios_base::in);
	if (!ifs.good()) vrb.error("Cannot open binary reference panel file: [" + filename + "]");

	try {
		boost::archive::binary_iarchive ia(ifs);
		ia >> H;
		ia >> V;
	} catch (std::exception& e) {
		std::stringstream err_str;
		err_str << "Problem reading binary reference panel (exception from boost archive). ";
		err_str << "Ensure you are using the same GLIMPSE and boost library version. ";
		err_str << e.what();
		vrb.error(err_str.str());
	}

	if (H.Ypacked.size() == 0) vrb.error("Problem reading binary file format. Empty PBWT detected.");
	vrb.bullet("Binary reference panel read (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

static std::string format_number(unsigned long n) {
	std::string s = std::to_string(n);
	int pos = s.length() - 3;
	while (pos > 0) {
		s.insert(pos, ",");
		pos -= 3;
	}
	return s;
}

static std::string format_bp(int bp) {
	double val = bp;
	if (val >= 1e9) return stb.str(val / 1e9, 2) + " Gbp";
	if (val >= 1e6) return stb.str(val / 1e6, 2) + " Mbp";
	if (val >= 1e3) return stb.str(val / 1e3, 2) + " kbp";
	return stb.str(bp) + " bp";
}

static std::string format_pct(unsigned long count, unsigned long total) {
	if (total == 0) return "0.0%";
	return stb.str(((double)count / total) * 100.0, 1) + "%";
}

void inspector::print_statistics() {
	std::string filename = options["input"].as < std::string > ();

	// File size
	struct stat st;
	long long file_size = 0;
	if (stat(filename.c_str(), &st) == 0) file_size = st.st_size;

	std::string size_str;
	if (file_size >= (long long)1024*1024*1024) size_str = stb.str((double)file_size / (1024.0*1024.0*1024.0), 2) + " GB";
	else if (file_size >= 1024*1024) size_str = stb.str((double)file_size / (1024.0*1024.0), 2) + " MB";
	else if (file_size >= 1024) size_str = stb.str((double)file_size / 1024.0, 2) + " KB";
	else size_str = stb.str(file_size) + " bytes";

	// Region stats (coords are 1-based inclusive, per htslib region-string convention)
	int input_span = V.input_stop - V.input_start + 1;
	int output_span = V.output_stop - V.output_start + 1;

	// Genetic map stats: scan all variants for cM range.
	// Note: the .bin format does not record whether split_reference was given a
	// real genetic map or fell back to the bp/1e6 synthetic map — both paths
	// baseline-subtract the first variant's cm, so we can't distinguish here.
	// We always print the observed cM range and let the caller interpret it.
	double input_min_cm = std::numeric_limits<double>::max();
	double input_max_cm = std::numeric_limits<double>::lowest();
	double output_min_cm = std::numeric_limits<double>::max();
	double output_max_cm = std::numeric_limits<double>::lowest();

	// Variant type counts (type field uses htslib VCF_* bitmask: VCF_SNP=1, VCF_MNP=2, VCF_INDEL=4, VCF_OTHER=8, etc.)
	unsigned long n_snps = 0, n_mnps = 0, n_indels = 0, n_other = 0;
	unsigned long n_lq = 0;
	unsigned long n_core = 0, n_buffer = 0;

	// Allele frequency bins (based on minor allele count / frequency)
	unsigned long af_singleton = 0;    // MAC == 1
	unsigned long af_mac_2_5 = 0;      // MAC 2-5
	unsigned long af_maf_lt_001 = 0;   // MAF < 0.01
	unsigned long af_maf_001_005 = 0;  // MAF 0.01-0.05
	unsigned long af_maf_005_050 = 0;  // MAF 0.05-0.50
	unsigned long af_monomorphic = 0;  // MAC == 0

	for (int i = 0; i < (int)V.vec_pos.size(); i++) {
		variant * v = V.vec_pos[i];

		// Genetic map
		if (v->cm < input_min_cm) input_min_cm = v->cm;
		if (v->cm > input_max_cm) input_max_cm = v->cm;
		if (v->bp >= V.output_start && v->bp <= V.output_stop) {
			if (v->cm < output_min_cm) output_min_cm = v->cm;
			if (v->cm > output_max_cm) output_max_cm = v->cm;
		}

		// Core vs buffer
		if (v->bp >= V.output_start && v->bp <= V.output_stop) n_core++;
		else n_buffer++;

		// Variant type (htslib bitmask)
		if (v->type & VCF_SNP) n_snps++;
		else if (v->type & VCF_MNP) n_mnps++;
		else if (v->type & VCF_INDEL) n_indels++;
		else n_other++;

		// Low quality
		if (v->LQ) n_lq++;

		// Allele frequency distribution. Use the per-variant allele-number
		// (cref + calt from the source VCF) as the denominator, so sites with
		// missing genotypes bucket correctly.
		unsigned int mac = v->getMAC();
		unsigned int an = v->cref + v->calt;
		double maf = (an > 0) ? (double)mac / an : 0.0;

		if (mac == 0) af_monomorphic++;
		else if (mac == 1) af_singleton++;
		else if (mac <= 5) af_mac_2_5++;
		else if (maf < 0.01) af_maf_lt_001++;
		else if (maf < 0.05) af_maf_001_005++;
		else af_maf_005_050++;
	}

	// Print everything
	vrb.title("Binary reference panel summary:");
	vrb.bullet("File                 : " + filename);
	vrb.bullet("File size            : " + size_str);
	vrb.bullet("Chromosome           : " + V.chrid);
	vrb.print("");

	vrb.bullet("Input region         : " + V.input_gregion + " (" + format_bp(input_span) + ")");
	vrb.bullet("Output region        : " + V.output_gregion + " (" + format_bp(output_span) + ")");

	if (input_min_cm <= input_max_cm)
		vrb.bullet("Genetic map (input)  : " + stb.str(input_min_cm, 4) + " - " + stb.str(input_max_cm, 4) + " cM (" + stb.str(input_max_cm - input_min_cm, 4) + " cM span)");
	if (output_min_cm <= output_max_cm)
		vrb.bullet("Genetic map (output) : " + stb.str(output_min_cm, 4) + " - " + stb.str(output_max_cm, 4) + " cM (" + stb.str(output_max_cm - output_min_cm, 4) + " cM span)");

	vrb.print("");
	vrb.bullet("Haplotypes           : " + format_number(H.n_ref_haps));
	vrb.bullet("Variants (total)     : " + format_number(H.n_tot_sites));
	vrb.bullet("  Common             : " + format_number(H.n_com_sites) + " (" + format_pct(H.n_com_sites, H.n_tot_sites) + ")");
	vrb.bullet("  Rare               : " + format_number(H.n_rar_sites) + " (" + format_pct(H.n_rar_sites, H.n_tot_sites) + ")");
	vrb.bullet("  Common HQ          : " + format_number(H.n_com_sites_hq) + " (" + format_pct(H.n_com_sites_hq, H.n_tot_sites) + ")");
	vrb.bullet("  Low quality        : " + format_number(n_lq) + " (" + format_pct(n_lq, H.n_tot_sites) + ")");

	vrb.print("");
	vrb.bullet("Variant types:");
	vrb.bullet("  SNPs               : " + format_number(n_snps) + " (" + format_pct(n_snps, H.n_tot_sites) + ")");
	if (n_mnps > 0)
		vrb.bullet("  MNPs               : " + format_number(n_mnps) + " (" + format_pct(n_mnps, H.n_tot_sites) + ")");
	vrb.bullet("  Indels             : " + format_number(n_indels) + " (" + format_pct(n_indels, H.n_tot_sites) + ")");
	if (n_other > 0)
		vrb.bullet("  Other              : " + format_number(n_other) + " (" + format_pct(n_other, H.n_tot_sites) + ")");

	vrb.print("");
	vrb.bullet("Allele frequency distribution:");
	if (af_monomorphic > 0)
		vrb.bullet("  Monomorphic        : " + format_number(af_monomorphic) + " (" + format_pct(af_monomorphic, H.n_tot_sites) + ")");
	vrb.bullet("  Singletons         : " + format_number(af_singleton) + " (" + format_pct(af_singleton, H.n_tot_sites) + ")");
	vrb.bullet("  MAC 2-5            : " + format_number(af_mac_2_5) + " (" + format_pct(af_mac_2_5, H.n_tot_sites) + ")");
	vrb.bullet("  MAF < 1%           : " + format_number(af_maf_lt_001) + " (" + format_pct(af_maf_lt_001, H.n_tot_sites) + ")");
	vrb.bullet("  MAF 1-5%           : " + format_number(af_maf_001_005) + " (" + format_pct(af_maf_001_005, H.n_tot_sites) + ")");
	vrb.bullet("  MAF 5-50%          : " + format_number(af_maf_005_050) + " (" + format_pct(af_maf_005_050, H.n_tot_sites) + ")");

	vrb.print("");
	vrb.bullet("Region breakdown:");
	vrb.bullet("  Core (output)      : " + format_number(n_core) + " variants");
	vrb.bullet("  Buffer only        : " + format_number(n_buffer) + " variants");
}

void inspector::inspect(std::vector < std::string > & args) {
	declare_options();
	parse_command_line(args);
	check_options();
	read_binary_panel();
	print_statistics();
}
