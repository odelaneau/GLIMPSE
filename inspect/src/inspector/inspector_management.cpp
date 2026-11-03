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

#include "../../versions/versions.h"
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
	// The .bin does not record whether split_reference used a genetic map, but without a usable
	// one every cM is exactly bp/1e6 minus the first variant's bp/1e6 (variant_map::setGeneticMap),
	// which interpolating a real map does not reproduce. Repeat the same arithmetic to compare.
	const double first_variant_cm = V.vec_pos.empty() ? 0.0 : V.vec_pos[0]->bp * 1.0 / 1e6;
	bool cm_is_bp_over_1e6 = true;
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
		if (v->cm != v->bp * 1.0 / 1e6 - first_variant_cm) cm_is_bp_over_1e6 = false;
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
	// All cM values are 0 with or without a map when every variant shares one position
	if (V.vec_pos.size() > 1 && V.vec_pos.front()->bp != V.vec_pos.back()->bp)
		vrb.bullet("Genetic map source   : " + std::string(cm_is_bp_over_1e6 ? "none, constant 1 cM/Mb" : "interpolated from a genetic map"));

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

void inspector::write_haplotypes() {
	const std::string fname = options["output"].as < std::string > ();
	const std::string fname_index = fname + ".csi";
	const bool include_buffers = options.count("include-buffers");

	// INFO/CM is anchored on the first core variant, so the overlapping buffers of adjacent
	// chunks can be used to put every chunk of a chromosome on one scale. Stored cM values
	// are already relative to the first variant of the bin, which is the fallback.
	double cm_offset = 0.0;
	bool core_has_variants = false;
	for (int i = 0 ; i < (int)V.vec_pos.size() && !core_has_variants ; i++) {
		if (V.vec_pos[i]->bp >= V.output_start && V.vec_pos[i]->bp <= V.output_stop) {
			cm_offset = V.vec_pos[i]->cm;
			core_has_variants = true;
		}
	}
	if (!core_has_variants) {
		if (include_buffers) vrb.warning("No variant in the core region: INFO/CM is relative to the first buffer variant instead");
		else vrb.warning("No variant in the core region: the output has no records");
	}

	vrb.wait("  * Writing reference haplotypes");
	tac.clock();

	// Rare alleles are stored per haplotype; invert to list the minor-allele carriers per site
	std::vector < std::vector < int > > rare_carriers(H.n_tot_sites);
	for (int h = 0 ; h < (int)H.n_ref_haps ; h++)
		for (const int i_site : H.ShapRef[h]) rare_carriers[i_site].push_back(h);

	htsFile * fp = hts_open(fname.c_str(), out_file_format.c_str());
	if (fp == NULL) vrb.error("Can't write to [" + fname + "]");
	if (out_file_indexed && options["threads"].as < int > () > 1) hts_set_threads(fp, options["threads"].as < int > ());

	// Every chunk of a panel gets an identical header, so the outputs concatenate with bcftools concat
	bcf_hdr_t * hdr = bcf_hdr_init("w");
	bcf_hdr_append(hdr, std::string("##source=GLIMPSE2_inspect v" + std::string(INSPECT_VERSION)).c_str());
	bcf_hdr_append(hdr, std::string("##contig=<ID=" + V.chrid + ">").c_str());
	bcf_hdr_append(hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"ALT allele count in the reference panel\">");
	bcf_hdr_append(hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in the reference panel\">");
	bcf_hdr_append(hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"ALT allele frequency in the reference panel\">");
	bcf_hdr_append(hdr, "##INFO=<ID=RARE,Number=0,Type=Flag,Description=\"Rare variant, stored sparsely by GLIMPSE2 (minor allele frequency below the --sparse-maf of GLIMPSE2_split_reference)\">");
	bcf_hdr_append(hdr, "##INFO=<ID=CM,Number=1,Type=Float,Description=\"Genetic position in cM relative to the first core-region variant of the GLIMPSE2 .bin file (1 cM/Mb if the panel was split without --map)\">");
	if (include_buffers) bcf_hdr_append(hdr, "##INFO=<ID=BUFFER,Number=0,Type=Flag,Description=\"Variant in a buffer region of the GLIMPSE2 .bin file, outside its core (output) region\">");
	bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Reference haplotype allele\">");

	// One haploid sample per haplotype, numbered from 1 and zero-padded so names sort in order
	const std::size_t name_width = std::to_string(H.n_ref_haps).size();
	for (unsigned int h = 0 ; h < H.n_ref_haps ; h++) {
		const std::string number = std::to_string(h + 1);
		const std::string name = "h" + std::string(name_width - number.size(), '0') + number;
		bcf_hdr_add_sample(hdr, name.c_str());
	}
	bcf_hdr_add_sample(hdr, NULL);
	if (bcf_hdr_write(fp, hdr)) vrb.error("Failed to write header to [" + fname + "]");
	if (out_file_indexed && bcf_idx_init(fp, hdr, 14, fname_index.c_str())) vrb.error("Failed to initialise index [" + fname_index + "]");

	bcf1_t * rec = bcf_init1();
	std::vector < int32_t > genotypes(H.n_ref_haps);
	const int32_t rid = bcf_hdr_name2id(hdr, V.chrid.c_str());
	unsigned long n_written = 0;
	int i_common = 0;
	for (int i = 0 ; i < (int)V.vec_pos.size() ; i++) {
		const variant * v = V.vec_pos[i];
		const bool is_common = H.flag_common[i];
		const bool in_core = v->bp >= V.output_start && v->bp <= V.output_stop;

		if (in_core || include_buffers) {
			bcf_clear1(rec);
			rec->rid = rid;
			rec->pos = v->bp - 1;
			bcf_update_id(hdr, rec, v->id.c_str());
			const char * alleles[2] = { v->ref.c_str(), v->alt.c_str() };
			bcf_update_alleles(hdr, rec, alleles, 2);

			const int32_t ac = v->calt;
			const int32_t an = v->cref + v->calt;
			const float af = (float)ac / an;
			const float cm = v->cm - cm_offset;
			bcf_update_info_int32(hdr, rec, "AC", &ac, 1);
			bcf_update_info_int32(hdr, rec, "AN", &an, 1);
			bcf_update_info_float(hdr, rec, "AF", &af, 1);
			if (!is_common) bcf_update_info_flag(hdr, rec, "RARE", NULL, 1);
			bcf_update_info_float(hdr, rec, "CM", &cm, 1);
			if (!in_core) bcf_update_info_flag(hdr, rec, "BUFFER", NULL, 1);

			if (is_common) {
				for (int h = 0 ; h < (int)H.n_ref_haps ; h++) genotypes[h] = bcf_gt_unphased(H.HvarRef.get(i_common, h));
			} else {
				std::fill(genotypes.begin(), genotypes.end(), bcf_gt_unphased(H.major_alleles[i]));
				for (const int h : rare_carriers[i]) genotypes[h] = bcf_gt_unphased(!H.major_alleles[i]);
			}
			bcf_update_genotypes(hdr, rec, genotypes.data(), H.n_ref_haps);
			if (bcf_write(fp, hdr, rec)) vrb.error("Failed to write record to [" + fname + "]");
			n_written++;
		}
		i_common += is_common;
	}
	bcf_destroy1(rec);

	if (out_file_indexed && bcf_idx_save(fp)) vrb.error("Failed to write index [" + fname_index + "]");
	if (hts_close(fp)) vrb.error("Failed to close [" + fname + "]");
	bcf_hdr_destroy(hdr);
	vrb.bullet("Wrote " + format_number(n_written) + " variants x " + format_number(H.n_ref_haps) + " haplotypes to [" + fname + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void inspector::inspect(std::vector < std::string > & args) {
	declare_options();
	parse_command_line(args);
	check_options();
	read_binary_panel();
	print_statistics();
	if (!out_file_format.empty()) write_haplotypes();
}
