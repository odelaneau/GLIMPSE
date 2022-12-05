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

#include <io/genotype_writer.h>
#include "../../versions/versions.h"

#define OFILE_VCFU	0
#define OFILE_VCFC	1
#define OFILE_BCFC	2

#define OFILE_BGENU	3
#define OFILE_BGENG	4
#define OFILE_BGENZ	5

genotype_writer::genotype_writer(const haplotype_set & _H, const genotype_set & _G, const variant_map & _V, const caller &  _C): H(_H), G(_G), V(_V), C(_C) {
}

genotype_writer::~genotype_writer() {
}

void genotype_writer::writeGenotypes(const std::string fname, OutputFormat output_fmt, OutputCompression output_compr, const int n_bits_bgen, const int n_main, const int n_threads, const std::string fai_fname) const {
	tac.clock();
	unsigned int file_type = OFILE_BCFC;
	std::string file_format = "wb";

	if (output_fmt == OutputFormat::VCF)
	{
		if (output_compr == OutputCompression::NONE) {file_format = "w"; file_type = OFILE_VCFU;}
		else {file_format = "wz"; file_type = OFILE_VCFC;}
	}

	htsFile * fp = hts_open(fname.c_str(),file_format.c_str());
	if (n_threads > 1) hts_set_threads(fp, n_threads);
	bcf_hdr_t * hdr = bcf_hdr_init("w");
	bcf1_t *rec = bcf_init1();

	std::string file_fmt = "";
	switch (file_type) {
		case OFILE_VCFU: file_fmt="VCF"; break;
		case OFILE_VCFC: file_fmt="VCF"; break;
		case OFILE_BCFC: file_fmt="BCF"; break;
	}

	// Create VCF header
	bcf_hdr_append(hdr, std::string("##fileDate="+tac.date()).c_str());
	bcf_hdr_append(hdr, std::string("##source=GLIMPSE_phase v" + std::string(PHASE_VERSION)).c_str());
	if (!fai_fname.empty()) update_header_from_fai(hdr,fai_fname);
	else bcf_hdr_append(hdr, std::string("##contig=<ID="+ V.chrid + ">").c_str());
	bcf_hdr_append(hdr, "##INFO=<ID=RAF,Number=A,Type=Float,Description=\"ALT allele frequency in the reference panel\">");
	bcf_hdr_append(hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"ALT allele frequency computed from DS/GP field across target samples\">");
	if (H.fploidy>0) bcf_hdr_append(hdr, "##INFO=<ID=INFO,Number=A,Type=Float,Description=\"IMPUTE info quality score for diploid samples\">");
	bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotypes\">");
	bcf_hdr_append(hdr, "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Genotype dosage\">");
	bcf_hdr_append(hdr, "##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Genotype posteriors\">");
	bcf_hdr_append(hdr, std::string("##NMAIN="+stb.str(n_main)).c_str());
	bcf_hdr_append(hdr, std::string("##FPLOIDY="+stb.str(H.fploidy)).c_str());

	//Add samples
	std::vector < int > ptr_gps = std::vector < int > (G.vecG.size(), 0);		// Pointers to iterate over sparse GPs
	for (int i = 0 ; i < G.vecG.size() ; i ++) bcf_hdr_add_sample(hdr, G.vecG[i]->name.c_str());
	bcf_hdr_add_sample(hdr, NULL);      // to update internal structures
	if (bcf_hdr_write(fp, hdr) )  vrb.error("Failed to write header to file [" + std::string(fname) + "]");

	//Add records
	int * genotypes = (int*)malloc(bcf_hdr_nsamples(hdr)*H.max_ploidy*sizeof(int));
	float * dosages = (float*)malloc(bcf_hdr_nsamples(hdr)*1*sizeof(float));
	float * posteriors = (float*)malloc(bcf_hdr_nsamples(hdr)*(H.max_ploidy+1)*sizeof(float));
	//int * haplotypes = (int*)malloc(bcf_hdr_nsamples(hdr)*sizeof(int));

	std::map<float,float*> map_ps;
	const std::size_t nsites = V.vec_pos.size();
	float prog_step = 1.0/nsites;
	float prog_bar = 0.0;
	int32_t chrid = bcf_hdr_name2id(hdr, V.chrid.c_str());

	for (int l = 0 ; l < nsites ; l ++) {
		// Clear current VCF record
		bcf_clear1(rec);

		// Update variant informations
		rec->rid = chrid; //bcf_hdr_name2id(hdr, chrid.c_str());
		rec->pos = V.vec_pos[l]->bp - 1;
		bcf_update_id(hdr, rec, V.vec_pos[l]->id.c_str());
		std::string alleles = V.vec_pos[l]->ref + "," + V.vec_pos[l]->alt;
		bcf_update_alleles_str(hdr, rec, alleles.c_str());

		// Store individual data
		float ds_sum=0.0, ds2_sum=0.0, ds4_sum=0.0;
		//std::fill(&genotypes[0], &genotypes[0] + H.max_ploidy*G.vecG.size(), bcf_gt_unphased(false));
		std::fill(&genotypes[0], &genotypes[0] + H.max_ploidy*G.vecG.size(), bcf_gt_phased(false));
		for (int i = 0 ; i < G.vecG.size() ; i++)
		{
			// Initialialize output data in case of Ref/Ref genotype
			float ds = 0.0f, gp0 = 1.0f, gp1 = 0.0f, gp2 = 0.0f;

			// Genotype is not mono
			if ((ptr_gps[i]<G.vecG[i]->stored_data.size()) && (G.vecG[i]->stored_data[ptr_gps[i]].idx==l))
			{
				if (G.vecG[i]->ploidy > 1)
				{
					gp0 = G.vecG[i]->stored_data[ptr_gps[i]].gp0;
					gp1 = G.vecG[i]->stored_data[ptr_gps[i]].gp1;
					gp2 = G.vecG[i]->stored_data[ptr_gps[i]].getGp2();
					ds = gp1 + 2.0f * gp2;
					if (gp1 > gp0 && gp1 > gp2)
					{
						const bool hds = G.vecG[i]->stored_data[ptr_gps[i]].hds;
						genotypes[H.max_ploidy*i+hds] = bcf_gt_phased(true);
						genotypes[H.max_ploidy*i+(1-hds)] = bcf_gt_phased(false);
					}
					else
					{
						genotypes[H.max_ploidy*i+0] = bcf_gt_phased(gp0 < gp2);
						genotypes[H.max_ploidy*i+1] = bcf_gt_phased(gp0 < gp2);
					}
				}
				else
				{
					gp0 = G.vecG[i]->stored_data[ptr_gps[i]].gp0;
					gp1 = G.vecG[i]->stored_data[ptr_gps[i]].gp1;
					ds = gp1;
					genotypes[H.max_ploidy*i+0] = bcf_gt_phased(gp0 < gp1);
					if (H.max_ploidy > 1) genotypes[H.max_ploidy*i+1] = bcf_int32_vector_end;
				}
				ptr_gps[i] ++;
			}

			// Store DS + GP rounded
			dosages[i] = std::roundf(ds * 1000.0) / 1000.0;
			posteriors[(H.max_ploidy+1)*i+0] = floorf(gp0 * 1000.0) / 1000.0;
			posteriors[(H.max_ploidy+1)*i+1] = floorf(gp1 * 1000.0) / 1000.0;
			map_ps.clear();
			map_ps.insert(std::make_pair(1.0f-(gp0-posteriors[(H.max_ploidy+1)*i+0]),&posteriors[(H.max_ploidy+1)*i+0]));
			map_ps.insert(std::make_pair(1.0f-(gp1-posteriors[(H.max_ploidy+1)*i+1]),&posteriors[(H.max_ploidy+1)*i+1]));

			if (H.max_ploidy>1)
			{
				if (G.vecG[i]->ploidy > 1)
				{
					posteriors[(H.max_ploidy+1)*i+2] = floorf(std::max(1.0f - (posteriors[(H.max_ploidy+1)*i+0]+posteriors[(H.max_ploidy+1)*i+1]), 0.0f)*1000.0)/1000.0;
					map_ps.insert(std::make_pair(1.0f-(gp2-posteriors[(H.max_ploidy+1)*i+2]),&posteriors[(H.max_ploidy+1)*i+2]));
				}
				else
				{
					genotypes[H.max_ploidy*i+1] = bcf_int32_vector_end;
					bcf_float_set(&posteriors[(H.max_ploidy+1)*i+2], bcf_float_vector_end);
				}
			}
			for (auto iter = map_ps.begin(); iter != map_ps.end() && std::accumulate(std::next(posteriors,(H.max_ploidy+1)*i+0), std::next(posteriors,(H.max_ploidy+1)*i+G.vecG[i]->ploidy+1), 0.0f)<0.9999f; ++iter) *iter->second += 0.001f;

			// Compute INFO/INFO statistics
			ds_sum += ds;
			if (H.fploidy == 2)
			{
				ds2_sum += ds * ds;
				ds4_sum += gp1 + 4.0*gp2;
			}
		}

		// Update INFO fields
		float freq_alt_refp = V.vec_pos[l]->calt * 1.0f / (V.vec_pos[l]->calt + V.vec_pos[l]->cref);
		float freq_alt_main = ds_sum / H.n_tar_haps;
		float infoscore = (H.fploidy == 2 && freq_alt_main>0.0 && freq_alt_main<1.0) ? (float)(1.0 - (ds4_sum - ds2_sum) / (H.n_tot_haps * freq_alt_main * (1.0 - freq_alt_main))) : 1.0f;
		infoscore = (infoscore<0.0f)?0.0f:infoscore;
		infoscore = roundf(infoscore * 1000.0) / 1000.0;
		bcf_update_info_float(hdr, rec, "RAF", &freq_alt_refp, 1);
		bcf_update_info_float(hdr, rec, "AF", &freq_alt_main, 1);
		if (H.fploidy>0) bcf_update_info_float(hdr, rec, "INFO", &infoscore, 1);

		// Update FORMAT fields
		bcf_update_genotypes(hdr, rec, genotypes, bcf_hdr_nsamples(hdr)*H.max_ploidy);
		bcf_update_format_float(hdr, rec, "DS", dosages, bcf_hdr_nsamples(hdr)*1);
		bcf_update_format_float(hdr, rec, "GP", posteriors, bcf_hdr_nsamples(hdr)*(H.max_ploidy+1));

		//Write record
		if (bcf_write1(fp, hdr, rec) ) vrb.error("Failed to write the record [" + V.chrid + ":" + std::to_string(V.vec_pos[l]->bp) + "_" + V.vec_pos[l]->ref + "_" + V.vec_pos[l]->alt +"] to file [" + std::string(fname) + "]");
		prog_bar+=prog_step;
		vrb.progress("  * " + file_fmt + " writing", prog_bar);
	}
	free(genotypes);
	free(dosages);
	free(posteriors);
	//free(haplotypes);
	bcf_destroy1(rec);
	bcf_hdr_destroy(hdr);
	if (hts_close(fp)) vrb.error("Non zero status when closing VCF/BCF file descriptor");

	switch (file_type) {
	case OFILE_VCFU: vrb.bullet("VCF writing [Uncompressed / N=" + stb.str(G.vecG.size()) + " / L=" + stb.str(V.size()) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_VCFC: vrb.bullet("VCF writing [Compressed / N=" + stb.str(G.vecG.size()) + " / L=" + stb.str(V.size()) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_BCFC: vrb.bullet("BCF writing [Compressed / N=" + stb.str(G.vecG.size()) + " / L=" + stb.str(V.size()) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	}

	vrb.bullet("Validating output and creating index");
	//create index using htslib (csi, using default bcftools option 14)
	if (!bcf_index_build3(fname.c_str(), NULL, 14, n_threads)) vrb.bullet("Index successfully created");
	else vrb.warning("Problem building the index for the output file. This can indicate a problem during the creation of the file. Try to build the index using tabix/bcftools.");
}

void genotype_writer::update_header_from_fai(bcf_hdr_t * hdr, const std::string fai_fname) const
{
    faidx_t *fai = fai_load3(fai_fname.c_str(),fai_fname.c_str(),NULL,FAI_FASTA);
    if ( !fai ) vrb.error("Could not parse " + fai_fname);

    // add any new contig lines
    int i, n = faidx_nseq(fai);
    for (i=0; i<n; i++)
    	bcf_hdr_append(hdr, std::string("##contig=<ID="+ std::string(faidx_iseq(fai,i)) + ",length=" + std::to_string(faidx_seq_len(fai,faidx_iseq(fai,i))) + ">").c_str());

    fai_destroy(fai);
}

#ifdef __BGEN__
void genotype_writer::writeGenotypesBgen(const std::string fname, OutputFormat oformat, OutputCompression ocompr, const int n_bits_bgen, const int n_main, const int n_threads) const
{
	// Init
	tac.clock();
	unsigned int file_type = 0;
	if (oformat == OutputFormat::BGEN && ocompr == OutputCompression::NONE) { file_type = OFILE_BGENU; }
	else if (oformat == OutputFormat::BGEN  && ocompr == OutputCompression::ZSTD) { file_type = OFILE_BGENZ; }
	else { file_type = OFILE_BGENG; }

	size_t n_sites_unbuf = 0;
	const std::size_t nsites = V.vec_pos.size();
	for (int l = 0 ; l < nsites ; l ++) if ((V.vec_pos[l]->bp >= V.output_start) && (V.vec_pos[l]->bp <= V.output_stop)) n_sites_unbuf++;

	std::vector < int > ptr_gps = std::vector < int > (G.vecG.size(), 0);		// Pointers to iterate over sparse GPs
	genfile::bgen::Context context;
	std::ios_base::openmode openmode = std::ostream::out;
	if (ocompr != OutputCompression::NONE) openmode |= std::ostream::binary;
	std::ofstream ofile(fname, openmode);
	initialiseBGEN(context, ocompr, n_sites_unbuf);
	initialiseStream(context,ofile);

	//Add records
	float prog_step = 1.0/n_sites_unbuf;
	float prog_bar = 0.0;
	for (int l = 0 ; l < nsites ; l ++) {
		if ((V.vec_pos[l]->bp < V.output_start) || (V.vec_pos[l]->bp > V.output_stop)) continue;
		std::vector< genfile::byte_t > buffer ;
		std::vector< genfile::byte_t > buffer2 ;
		const variant* curr_variant = V.vec_pos[l];
		genfile::byte_t* end =  genfile::bgen::write_snp_identifying_data(
				&buffer,
				context,
				V.chrid + ":" + stb.str(curr_variant->bp) + "_" + curr_variant->ref + "_" + curr_variant->alt,
				curr_variant->id,
				V.chrid,
				curr_variant->bp,
				2, [&curr_variant]( std::size_t i ) { if( i == 0 ) { return curr_variant->ref; } else { return curr_variant->alt ; }}
		);

		ofile.write( reinterpret_cast< char* >( &buffer[0] ), end - &buffer[0] ) ;

		// Store individual data
		genfile::bgen::GenotypeDataBlockWriter writer( &buffer, &buffer2, context, n_bits_bgen);
		writer.initialise(G.vecG.size(), 2 ) ; //nsamples, nalleles

		for (int i = 0 ; i < G.vecG.size() ; i++)
		{
			// Initialialize output data in case of Ref/Ref genotype
			float gp0 = 1.0f, gp1 = 0.0f, gp2 = 0.0f;

			// Genotype is NOT 100% certain Ref/Ref
			if ((ptr_gps[i]<G.vecG[i]->stored_data.size()) && (G.vecG[i]->stored_data[ptr_gps[i]].idx==l))
			{
				gp0 = G.vecG[i]->stored_data[ptr_gps[i]].gp0;
				gp1 = G.vecG[i]->stored_data[ptr_gps[i]].gp1;
				if (G.vecG[i]->ploidy > 1) gp2 = G.vecG[i]->stored_data[ptr_gps[i]].getGp2();
				ptr_gps[i] ++;
			}
			writer.set_sample(i) ;
			writer.set_number_of_entries( G.vecG[i]->ploidy, G.vecG[i]->ploidy+1, genfile::ePerUnorderedGenotype, genfile::eProbability ) ; //ploidy, ngenotypes
			writer.set_value( 0, gp0);
			writer.set_value( 1, gp1);
			if (G.vecG[i]->ploidy > 1) writer.set_value( 2, gp2);
		}
		writer.finalise() ;
		ofile.write( reinterpret_cast< char const* >( writer.repr().first ), writer.repr().second - writer.repr().first ) ;
		//Write record
		prog_bar+=prog_step;
		vrb.progress("  * BGEN writing", prog_bar);
	}
	ofile.close();
	switch (file_type) {
	case OFILE_VCFU: vrb.bullet("VCF writing [Uncompressed / N=" + stb.str(G.vecG.size()) + " / L=" + stb.str(V.size()) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;

		case OFILE_BGENU: vrb.bullet("BGEN writing completed. [Uncompressed / N=" + stb.str(G.vecG.size()) + " / L=" + stb.str(n_sites_unbuf) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
		case OFILE_BGENG: vrb.bullet("BGEN writing completed. [ZLIB compressed / N=" + stb.str(G.vecG.size()) + " / L=" + stb.str(n_sites_unbuf) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
		case OFILE_BGENZ: vrb.bullet("BGEN writing completed. [ZSTD compressed / N=" + stb.str(G.vecG.size()) + " / L=" + stb.str(n_sites_unbuf) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	}
}

void genotype_writer::initialiseBGEN(genfile::bgen::Context& context,OutputCompression ocompr, const size_t n_sites_unbuf) const
{
	context.number_of_variants = n_sites_unbuf;
	context.number_of_samples = G.vecG.size();

	//change the values of this line to change layout and compression method
	genfile::bgen::Compression compression;
	switch (ocompr)
	{
		case OutputCompression::NONE :
			compression = genfile::bgen::Compression::e_NoCompression;
			break;
		case OutputCompression::ZLIB :
			compression = genfile::bgen::Compression::e_ZlibCompression;
			break;
		case OutputCompression::ZSTD :
			compression = genfile::bgen::Compression::e_ZstdCompression;
			break;
	}
	context.flags = genfile::bgen::e_Layout2 | compression;
};

void genotype_writer::initialiseStream(genfile::bgen::Context& context,std::ostream& oStream) const
{
	context.flags |= genfile::bgen::e_SampleIdentifiers;
	set_sample_names_impl(context,oStream);
};

void genotype_writer::set_sample_names_impl(genfile::bgen::Context& context,std::ostream& oStream) const
{
	std::vector<std::string> sample_names(G.vecG.size());
	for (int i = 0 ; i < G.vecG.size() ; i ++) sample_names[i] = G.vecG[i]->name;

	uint32_t offset = context.header_size() ;
	update_offset_and_header_block(context,oStream, offset) ;
	if( context.flags && genfile::bgen::e_SampleIdentifiers ) {
		offset += genfile::bgen::write_sample_identifier_block(
			oStream,
			context,
			sample_names
		) ;
	}
	update_offset_and_header_block(context,oStream, offset) ;
	oStream.seekp( offset+4, std::ios_base::beg ) ;

};

void genotype_writer::update_offset_and_header_block(genfile::bgen::Context& context,std::ostream& oStream, uint32_t offset) const
{
	oStream.seekp( 0, std::ios_base::beg ) ;
	if( !oStream.bad() ) {
		genfile::bgen::write_offset( oStream, offset ) ;
		genfile::bgen::write_header_block(
				oStream,
				context
		) ;
	}
}
#endif
