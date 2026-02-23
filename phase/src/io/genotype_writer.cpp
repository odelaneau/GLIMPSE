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

void genotype_writer::writeGenotypes(const argument_set& A) const {
	tac.clock();
	unsigned int file_type = OFILE_BCFC;
	std::string file_format = "wb";

	if (A.mOutputFormat == OutputFormat::VCF)
	{
		if (A.mOutputCompression == OutputCompression::NONE) {file_format = "w"; file_type = OFILE_VCFU;}
		else {file_format = "wz"; file_type = OFILE_VCFC;}
	}

	htsFile * fp = hts_open(A.mOutputFilename.c_str(),file_format.c_str());
	if (A.mNumThreads > 1) hts_set_threads(fp, A.mNumThreads);
	bcf_hdr_t * hdr = bcf_hdr_init("w");
	bcf1_t *rec = bcf_init1();

	std::string file_fmt;
	switch (file_type) {
		case OFILE_VCFU: file_fmt="VCF"; break;
		case OFILE_VCFC: file_fmt="VCF"; break;
		case OFILE_BCFC: file_fmt="BCF"; break;
	}
	int32_t pass=bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");

	// Create VCF header
	bcf_hdr_append(hdr, std::string("##fileDate="+tac.date()).c_str());
	bcf_hdr_append(hdr, std::string("##source=GLIMPSE2_phase_v" + std::string(PHASE_VERSION)).c_str());
	if (!A.mContigFaiFilename.empty()) update_header_from_fai(hdr,A.mContigFaiFilename);
	else
	{
		if (H.contigs_header.size()>0)
		{
		    for (int i = 0; i < H.contigs_header.size(); ++i)
		    {
		    	bcf_hdr_append(hdr, H.contigs_header[i].c_str());
		    }
		}
		else bcf_hdr_append(hdr, std::string("##contig=<ID="+ V.chrid + ">").c_str());
	}
	bcf_hdr_append(hdr, "##INFO=<ID=RAF,Number=A,Type=Float,Description=\"ALT allele frequency in the reference panel\">");
	bcf_hdr_append(hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"ALT allele frequency computed from DS/GP field across target samples\">");
	bcf_hdr_append(hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"ALT allele count from GT across target samples\">");
	bcf_hdr_append(hdr, "##INFO=<ID=AN,Number=A,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
	bcf_hdr_append(hdr, "##INFO=<ID=INFO,Number=A,Type=Float,Description=\"IMPUTE info quality score for diploid samples\">");
	bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotypes\">");

	if (A.mPrintDS) bcf_hdr_append(hdr, "##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"Genotype dosage\">");
	if (A.mPrintGP)	bcf_hdr_append(hdr, "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype posterior probabilities\">");
	if (A.mPrintAP) bcf_hdr_append(hdr, "##FORMAT=<ID=AP,Number=2,Type=Float,Description=\"ALT allele probability of first haplotype and second haplotype\">");

	bcf_hdr_append(hdr, std::string("##NMAIN="+stb.str(A.mMain)).c_str());
	bcf_hdr_append(hdr, std::string("##FPLOIDY="+stb.str(H.fploidy)).c_str());

	//Add samples
	std::vector < int > ptr_gps = std::vector < int > (G.vecG.size(), 0);		// Pointers to iterate over sparse GPs
	for (int i = 0 ; i < G.vecG.size() ; i ++) bcf_hdr_add_sample(hdr, G.vecG[i]->name.c_str());
	bcf_hdr_add_sample(hdr, NULL);      // to update internal structures
	if (bcf_hdr_write(fp, hdr) )  vrb.error("Failed to write header to file [" + std::string(A.mContigFaiFilename) + "]");

	//Add records
	int * genotypes = (int*)malloc(bcf_hdr_nsamples(hdr)*H.max_ploidy*sizeof(int));
	float * dosages = (float*)malloc(bcf_hdr_nsamples(hdr)*1*sizeof(float));
	float * posteriors = (float*)malloc(bcf_hdr_nsamples(hdr)*(H.max_ploidy+1)*sizeof(float));
	float * alt_ap;
	if (A.mPrintAP) alt_ap = (float*)malloc(bcf_hdr_nsamples(hdr)*(H.max_ploidy)*sizeof(float));//max_ploidy
	//int * haplotypes = (int*)malloc(bcf_hdr_nsamples(hdr)*sizeof(int));

	unsigned int ac = 0;
	unsigned int an = H.n_tar_haps;
	float freq_alt_main = 0.0f;
	float infoscore = 0.0;
	const uint8_t ploidy_const2 = (H.max_ploidy > 1) ? 2:0;

	std::map<float,float*> map_ps;
	const std::size_t nsites = V.vec_pos.size();
	float prog_step = 1.0/nsites;
	float prog_bar = 0.0;
	int32_t chrid = bcf_hdr_name2id(hdr, V.chrid.c_str());
	std::array<float,3> gp;
	const double d_n_tar_haps= (double)H.n_tar_haps;

	for (int l = 0 ; l < nsites ; l ++)
	{
		double ds_sum=0.0, ds2_sum=0.0, ds4_sum=0.0;
		for (int i = 0 ; i < G.vecG.size() ; i++)
		{
			if ((ptr_gps[i]<G.vecG[i]->stored_data.size()) && (G.vecG[i]->stored_data[ptr_gps[i]].idx==l))
			{
				if (G.vecG[i]->ploidy > 1)
				{
					gp[0] = G.vecG[i]->stored_data[ptr_gps[i]].gp0;
					gp[1] = G.vecG[i]->stored_data[ptr_gps[i]].gp1;
					gp[2] = G.vecG[i]->stored_data[ptr_gps[i]].getGp2();
					if (gp[1] > gp[0] && gp[1] > gp[2])
					{
						const bool hds = G.vecG[i]->stored_data[ptr_gps[i]].hds;
						genotypes[H.max_ploidy*i+hds] = bcf_gt_phased(true);
						genotypes[H.max_ploidy*i+(1-hds)] = bcf_gt_phased(false);
						ac += 1;
					}
					else
					{
						genotypes[H.max_ploidy*i+0] = bcf_gt_phased(gp[0] < gp[2]);
						genotypes[H.max_ploidy*i+1] = bcf_gt_phased(gp[0] < gp[2]);
						ac += 2*(gp[0] < gp[2]);
					}
				}
				else
				{
					gp[0] = G.vecG[i]->stored_data[ptr_gps[i]].gp0;
					gp[1] = G.vecG[i]->stored_data[ptr_gps[i]].gp1;
					gp[2] = 0.0f;
					genotypes[H.max_ploidy*i+0] = bcf_gt_phased(gp[0] < gp[1]);
					if (H.max_ploidy > 1) genotypes[H.max_ploidy*i+1] = bcf_int32_vector_end;
					ac += (gp[0] < gp[2]);
				}
				if (A.mPrintAP)
				{
					for (size_t j=0; j<G.vecG[i]->ploidy; ++j) { 
                      double temp = static_cast<double>(G.vecG[i]->stored_alt_ds[ptr_gps[i]].ap[j]) * 10000.0;
					  double rounded = std::round(temp);
					  alt_ap[G.vecG[i]->ploidy*i + j] =  static_cast<float>(rounded / 10000.0); 
					  // std::roundf(10000.0f * G.vecG[i]->stored_alt_ds[ptr_gps[i]].ap[j])/10000.0f;
					}
				}
				ptr_gps[i] ++;
			}
			else
			{
				gp = {1.0f,0.0f,0.0f};
				genotypes[H.max_ploidy*i+0] = bcf_gt_phased(false);
				if (H.max_ploidy > 1) genotypes[H.max_ploidy*i+1] = bcf_gt_phased(false);
				if (A.mPrintAP)
				{
					alt_ap[H.max_ploidy*i + 0] = 0.0f;
					alt_ap[H.max_ploidy*i + 1] = 0.0f;
				}
			}
			dosages[i] = std::roundf((gp[1] + 2 * gp[2]) * 1000.0f) / 1000.0f;
			for (int j=0; j<H.max_ploidy;++j) posteriors[(H.max_ploidy+1)*i+j] = ::floorf(gp[j] * 1000.0f) / 1000.0f;
			map_ps.clear();
			map_ps.insert(std::make_pair(1.0f-(gp[0]-posteriors[(H.max_ploidy+1)*i+0]),&posteriors[(H.max_ploidy+1)*i+0]));//max_ploidy
			map_ps.insert(std::make_pair(1.0f-(gp[1]-posteriors[(H.max_ploidy+1)*i+1]),&posteriors[(H.max_ploidy+1)*i+1]));//max_ploidy
			if (H.max_ploidy>1)
			{
				if (G.vecG[i]->ploidy > 1)
				{
					posteriors[(H.max_ploidy+1)*i+2] = floorf(std::max(1.0f - (posteriors[(H.max_ploidy+1)*i+0]+posteriors[(H.max_ploidy+1)*i+1]), 0.0f)*1000.0)/1000.0;
					map_ps.insert(std::make_pair(1.0f-(gp[2]-posteriors[(H.max_ploidy+1)*i+2]),&posteriors[(H.max_ploidy+1)*i+2]));
				}
				else
				{
					genotypes[H.max_ploidy*i+1] = bcf_int32_vector_end;
					bcf_float_set(&posteriors[(H.max_ploidy+1)*i+2], bcf_float_vector_end);
				}
			}
			for (auto iter = map_ps.begin(); iter != map_ps.end() && std::accumulate(std::next(posteriors,(H.max_ploidy+1)*i+0), std::next(posteriors,(H.max_ploidy+1)*i+G.vecG[i]->ploidy+1), 0.0f)<0.999f; ++iter) *iter->second += 0.001f; //max_ploidy
			ds_sum += (double) posteriors[(H.max_ploidy+1)*i+1] + ploidy_const2*posteriors[(H.max_ploidy+1)*i+ploidy_const2];
			ds2_sum += (double) ( posteriors[(H.max_ploidy+1)*i+1] + ploidy_const2*posteriors[(H.max_ploidy+1)*i+ploidy_const2]) * ( posteriors[(H.max_ploidy+1)*i+1] + ploidy_const2*posteriors[(H.max_ploidy+1)*i+ploidy_const2]);
			ds4_sum += (double)  posteriors[(H.max_ploidy+1)*i+1] + 4*ploidy_const2*posteriors[(H.max_ploidy+1)*i+ploidy_const2];
		}
		float freq_alt_refp = V.vec_pos[l]->calt * 1.0f / (V.vec_pos[l]->calt + V.vec_pos[l]->cref);
		double theta = ds_sum / d_n_tar_haps;
		infoscore = (theta>0 && theta<1) ? (float)(1 - (ds4_sum - ds2_sum) / (d_n_tar_haps * theta * (1.0 - theta))) : 1;
		infoscore = roundf((infoscore * 1000.0f)) / 1000.0f;

		bcf_clear1(rec);
		rec->rid = chrid; //bcf_hdr_name2id(hdr, chrid.c_str());
		rec->pos = V.vec_pos[l]->bp - 1;
		std::string alleles = V.vec_pos[l]->ref + "," + V.vec_pos[l]->alt;
		bcf_update_id(hdr, rec, V.vec_pos[l]->id.c_str());
		bcf_update_alleles_str(hdr, rec, alleles.c_str());
		bcf_update_filter(hdr,rec,&pass,1);
		bcf_update_info_float(hdr, rec, "RAF", &freq_alt_refp, 1);
		bcf_update_info_int32(hdr, rec, "AC", &ac, 1);
		bcf_update_info_int32(hdr, rec, "AN", &an, 1);
		bcf_update_info_float(hdr, rec, "AF", &freq_alt_main, 1);
		bcf_update_info_float(hdr, rec, "INFO", &infoscore, 1);//fploidy
		bcf_update_genotypes(hdr, rec, genotypes, bcf_hdr_nsamples(hdr)*H.max_ploidy);
		if (A.mPrintDS) bcf_update_format_float(hdr, rec, "DS", dosages, bcf_hdr_nsamples(hdr)*1);
		if (A.mPrintAP) bcf_update_format_float(hdr, rec, "AP", alt_ap, bcf_hdr_nsamples(hdr)*H.max_ploidy);
		if (A.mPrintGP) bcf_update_format_float(hdr, rec, "GP", posteriors, bcf_hdr_nsamples(hdr)*(H.max_ploidy+1));
		if (bcf_write1(fp, hdr, rec) ) vrb.error("Failed to write the record [" + V.chrid + ":" + std::to_string(V.vec_pos[l]->bp) + "_" + V.vec_pos[l]->ref + "_" + V.vec_pos[l]->alt +"] to file [" + std::string(A.mOutputFilename) + "]");
		prog_bar+=prog_step;
		vrb.progress("  * " + file_fmt + " writing", prog_bar);
	}
	free(genotypes);
	free(dosages);
	if (A.mPrintGP) free(posteriors);
	if (A.mPrintAP)
	{
		free(alt_ap);
	}
	//free(haplotypes);
	bcf_destroy1(rec);
	bcf_hdr_destroy(hdr);
	if (hts_close(fp)) vrb.error("Non zero status when closing VCF/BCF file descriptor");

	switch (file_type) {
	case OFILE_VCFU: vrb.bullet("VCF writing \t[done]\t\t[RAW  / N=" + stb.str(G.vecG.size()) + " / L=" + stb.str(V.size()) + "]\t(" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_VCFC: vrb.bullet("VCF writing \t[done]\t\t[ZLIB / N=" + stb.str(G.vecG.size()) + " / L=" + stb.str(V.size()) + "]\t(" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_BCFC: vrb.bullet("BCF writing \t[done]\t\t[ZLIB / N=" + stb.str(G.vecG.size()) + " / L=" + stb.str(V.size()) + "]\t(" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	}

	if (A.mOutputIndex && file_type!=OFILE_VCFU)
	{
		tac.clock();
		vrb.progress("  * Indexing output");
		if (!bcf_index_build3(A.mOutputFilename.c_str(), NULL, 14, A.mNumThreads)) vrb.bullet("Indexing output\t[done]\t\t\t\t\t(" + stb.str(tac.rel_time()*0.001, 2) + "s)");
		else vrb.warning("Problem building the index for the output file. This can indicate a problem during the creation of the file. Try to build the index using tabix/bcftools.");
	}
	if (A.mOutputIndex && file_type==OFILE_VCFU)
		vrb.warning("Plain VCF file format with no compression as output: no CSI index can be created. Uncompressed VCF file format is strongly discouraged. Please use vcf.gz or bcf file format.");
}

void genotype_writer::update_header_from_fai(bcf_hdr_t * hdr, const std::string fai_fname) const
{
    faidx_t *fai = fai_load3(fai_fname.c_str(),fai_fname.c_str(),NULL,FAI_FASTA);
    if ( !fai ) vrb.error("Could not parse " + fai_fname);

    // add any new contig lines
    int i, n = faidx_nseq(fai);
    for (i=0; i<n; i++)
    {
    	std::string length = "";
    	if (faidx_seq_len(fai,faidx_iseq(fai,i)) > 0) length = ",length=" + std::to_string(faidx_seq_len(fai,faidx_iseq(fai,i)));
    	bcf_hdr_append(hdr, std::string("##contig=<ID="+ std::string(faidx_iseq(fai,i)) + length + ">").c_str());
    }
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
