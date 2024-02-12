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

#include <ligater/ligater_header.h>
#include <htslib/khash.h>
#include <sys/stat.h>
#include <utils/otools.h>
#include <utils/basic_stats.h>


#define OFILE_VCFU	0
#define OFILE_VCFC	1
#define OFILE_BCFC	2

#define GET(n,i)	(((n)>>(i))&1U)
#define TOG(n,i)	((n)^=(1UL<<(i)))

#define SWAP(type_t, a, b) { type_t t = a; a = b; b = t; }
KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

void ligater::remove_info(bcf_hdr_t *hdr, bcf1_t *line)
{
    // remove all INFO fields
    if ( !(line->unpacked & BCF_UN_INFO) ) bcf_unpack(line, BCF_UN_INFO);

    for (int i=0; i<line->n_info; i++)
    {
        bcf_info_t *inf = &line->d.info[i];
        const char *key = bcf_hdr_int2id(hdr,BCF_DT_ID,inf->key);
        if ( key[0]=='A' && key[1]=='N' && !key[2] ) continue;
        if ( key[0]=='A' && key[1]=='C' && !key[2] ) continue;

        if ( inf->vptr_free )
        {
            free(inf->vptr - inf->vptr_off);
            inf->vptr_free = 0;
        }
        line->d.shared_dirty |= BCF1_DIRTY_INF;
        inf->vptr = NULL;
        inf->vptr_off = inf->vptr_len = 0;
    }
}

void ligater::remove_format(bcf_hdr_t *hdr, bcf1_t *line)
{
    // remove all FORMAT fields except GT
    if ( !(line->unpacked & BCF_UN_FMT) ) bcf_unpack(line, BCF_UN_FMT);

    for (int i=0; i<line->n_fmt; i++)
    {
        bcf_fmt_t *fmt = &line->d.fmt[i];
        const char *key = bcf_hdr_int2id(hdr,BCF_DT_ID,fmt->id);
        if ( key[0]=='G' && key[1]=='T' && !key[2] ) continue;

        if ( fmt->p_free )
        {
            free(fmt->p - fmt->p_off);
            fmt->p_free = 0;
        }
        line->d.indiv_dirty = 1;
        fmt->p = NULL;
    }
}

void ligater::phase_update(bcf_hdr_t *hdr, bcf1_t *line, const bool uphalf)
{
	int nGTs = bcf_get_genotypes(hdr, line, &GTa, &mGTa);
	if ( nGTs <= 0 ) return;    // GT field is not present
	for (auto i=0; i<bcf_hdr_nsamples(hdr); i++)
    {
		if ( !swap_phase[uphalf][i] ) continue;
		int *gt = &GTa[i*2];
		if ( bcf_gt_is_missing(gt[0]) || gt[1]==bcf_int32_vector_end ) continue;
		if (!bcf_gt_is_phased(gt[1])) continue;
        const int gt0 = bcf_gt_phased(bcf_gt_allele(gt[1])==1);
        const int gt1 = bcf_gt_phased(bcf_gt_allele(gt[0])==1);
        gt[0] = gt0;
        gt[1] = gt1;
    }
    bcf_update_genotypes(hdr,line,GTa,nGTs);
}

void ligater::update_distances()
{
	for (int i = 0 ; i < nsamples; i++)
	{
	    int *gta = &GTa[i*2];
	    int *gtb = &GTb[i*2];
	    if ( gta[1]==bcf_int32_vector_end || gtb[1]==bcf_int32_vector_end ) continue;
	    if ( bcf_gt_is_missing(gta[0]) || bcf_gt_is_missing(gta[1]) || bcf_gt_is_missing(gtb[0]) || bcf_gt_is_missing(gtb[1]) ) continue;
	    if ( !bcf_gt_is_phased(gta[1]) || !bcf_gt_is_phased(gtb[1]) ) continue;
	    if ( bcf_gt_allele(gta[0])==bcf_gt_allele(gta[1]) || bcf_gt_allele(gtb[0])==bcf_gt_allele(gtb[1]) ) continue;
	    if ( bcf_gt_allele(gta[0])==bcf_gt_allele(gtb[0]) && bcf_gt_allele(gta[1])==bcf_gt_allele(gtb[1]) )
	    {
	        if ( swap_phase[0][i] ) nmism[i]++; else nmatch[i]++;
	    }
	    if ( bcf_gt_allele(gta[0])==bcf_gt_allele(gtb[1]) && bcf_gt_allele(gta[1])==bcf_gt_allele(gtb[0]) )
	    {
	        if ( swap_phase[0][i] ) nmatch[i]++; else nmism[i]++;
	    }
	}
}

void ligater::write_record(htsFile *fd, bcf_hdr_t * out_hdr, bcf_hdr_t * hdr_in, bcf1_t *line, const bool uphalf)
{
	if ( nswap[uphalf] ) phase_update(hdr_in, line, uphalf);
	//remove_info(out_hdr,line);
	//remove_format(out_hdr,line);
	bcf_translate(out_hdr, hdr_in, line);
	if (bcf_write(fd, out_hdr, line) ) vrb.error("Failed to write the record output to file");
}

void ligater::scan_overlap(const int ifname, const char * seek_chr, int seek_pos)
{
	bcf_srs_t * sr =  bcf_sr_init();
	sr->require_index = 1;
	sr->collapse = COLLAPSE_NONE;
	sr->max_unpack = BCF_UN_FMT;

	int n_threads = options["threads"].as < int > ();
	if (n_threads > 1) bcf_sr_set_threads(sr, n_threads);

	if (!bcf_sr_add_reader (sr, filenames[ifname-2].c_str())) vrb.error("Problem opening/creating index file for [" + filenames[ifname-2] + "]");
	if (!bcf_sr_add_reader (sr, filenames[ifname-1].c_str())) vrb.error("Problem opening/creating index file for [" + filenames[ifname-1] + "]");

	int nset = 0;
	int n_sites_buff = 0;
	int n_sites_tot = 0;
	int last_pos=seek_pos;

	bcf1_t * line0 = NULL, * line1 = NULL;
	bcf_sr_seek(sr, seek_chr, seek_pos);
	while ((nset = bcf_sr_next_line (sr)))
	{
		if (nset==1)
		{
			if ( !bcf_sr_has_line(sr,0) && bcf_sr_region_done(sr,0)) break;  // no input from the first reader
			++n_sites_tot;
			continue;
		}
		line0 =  bcf_sr_get_line(sr, 0);
		if (line0->n_allele != 2) continue;

		line1 =  bcf_sr_get_line(sr, 1);
		//bcf_unpack(line0,BCF_UN_FMT);
		//bcf_unpack(line1,BCF_UN_FMT);
		int nGTsa = bcf_get_genotypes(sr->readers[0].header, line0, &GTa, &mGTa);
		int nGTsb = bcf_get_genotypes(sr->readers[1].header, line1, &GTb, &mGTb);
		if ( nGTsa <= 0 || nGTsb <= 0 )
			vrb.error("GT field is not present in overlap at position: " + std::to_string(line0->pos + 1));

		update_distances();
		last_pos = line0->pos;
		++n_sites_buff;
		++n_sites_tot;
	}
	bcf_sr_destroy(sr);

	stats1D stats_all;
	stats1D phaseq;

	nswap[1]=0;
	for (int i = 0 ; i < nsamples; i++)
	{
		swap_phase[1][i] = nmatch[i] < nmism[i];
		nswap[1] += swap_phase[1][i];

		stats_all.push(nmatch[i] + nmism[i]);

		float f = 99;
        if ( nmatch[i] && nmism[i] )
        {
            // Entropy-inspired quality. The factor 0.7 shifts and scales to (0,1)
           float f0 = (float)nmatch[i]/(nmatch[i]+nmism[i]);
           f = (99*(0.7 + f0*logf(f0) + (1-f0)*logf(1-f0))/0.7);
        }
        phaseq.push(f);

		nmatch[i] = 0;
		nmism[i]  = 0;
	}

	if (n_sites_buff <=0) vrb.error("Overlap is empty");
	nsites_buff_d2.push_back(n_sites_buff/2);
	vrb.print("Buf " + stb.str(nsites_buff_d2.size() -1) + " ["+std::string(seek_chr)+":"+stb.str(seek_pos+1)+"-"+stb.str(last_pos+1)+"] [L_isec=" + stb.str(n_sites_buff) + " / L_tot=" + stb.str(n_sites_tot) + "] [Avg #hets=" + stb.str(stats_all.mean()) + "] [Switch rate=" + stb.str(nswap[1]*1.0 / nsamples) + "] [Avg phaseQ=" + stb.str(phaseq.mean()) + "]");
}

void ligater::ligate() {
	tac.clock();
	vrb.title("Ligating chunks");

	//Create all input file descriptors
	vrb.bullet("Creating file descriptor");

	std::string file_format = "w";
	std::string fname = options["output"].as < std::string > ();
	unsigned int file_type = OFILE_VCFU;
	if (fname.size() > 6 && fname.substr(fname.size()-6) == "vcf.gz") { file_format = "wz"; file_type = OFILE_VCFC; }
	if (fname.size() > 3 && fname.substr(fname.size()-3) == "bcf") { file_format = "wb"; file_type = OFILE_BCFC; }
	bcf_srs_t * sr =  bcf_sr_init();
	sr->require_index = 1;
	int n_threads = options["threads"].as < int > ();
	if (n_threads > 1) if (bcf_sr_set_threads(sr, n_threads) < 0) vrb.error("Failed to create threads");

	bcf_hdr_t * out_hdr = NULL;
	bcf1_t *line_t = bcf_init();
	std::vector<int> start_pos(nfiles);

	for (int f = 0, prev_chrid = -1 ; f < nfiles ; f ++)
	{
		htsFile *fp = hts_open(filenames[f].c_str(), "r"); if ( !fp ) vrb.error("Failed to open: " + filenames[f] + ".");
		bcf_hdr_t *hdr = bcf_hdr_read(fp); if ( !hdr ) vrb.error("Failed to parse header: " + filenames[f] +".");
		out_hdr = bcf_hdr_merge(out_hdr,hdr);
        if ( bcf_hdr_nsamples(hdr) != bcf_hdr_nsamples(out_hdr) ) vrb.error("Different number of samples in " + filenames[f] + ".");
        for (int j=0; j<bcf_hdr_nsamples(hdr); j++)
        	if ( std::string(out_hdr->samples[j]) != std::string(hdr->samples[j]) )  vrb.error("Different sample names in " + filenames[f] + ".");

        int ret = bcf_read(fp, hdr, line_t);
		if ( ret!=0 ) vrb.error("Empty file detected: " + filenames[f] +".");
        else
        {
            int chrid = bcf_hdr_id2int(out_hdr,BCF_DT_CTG,bcf_seqname(hdr,line_t));
            start_pos[f] = chrid==prev_chrid ? line_t->pos : -1;
            prev_chrid = chrid;
        }
        bcf_hdr_destroy(hdr);
        if ( hts_close(fp)!=0 ) vrb.error("Close failed: " + filenames[f] + ".");
	}
	bcf_destroy(line_t);

    for (int i=1; i<nfiles; i++) if ( start_pos[i-1]!=-1 && start_pos[i]!=-1 && start_pos[i]<start_pos[i-1] ) vrb.error("The files not in ascending order");
    int i = 0, nrm = 0;
    /*
    while ( i<out_hdr->nhrec )
    {
    	bcf_hrec_t *hrec = out_hdr->hrec[i];
		// remove everything in INFO and FORMAT except INFO/AC, INFO/AN and FORMAT/GT
    	const std::string strk(hrec->key);
    	if (strk != "INFO" && strk!="FORMAT") { i++; continue; }

		int id = bcf_hrec_find_key(hrec, "ID");
		if ( id>=0 )
		{
			if (strk=="FORMAT" && std::string(hrec->vals[id])=="GT" ) { i++; continue; }
			if (strk=="INFO" && (std::string(hrec->vals[id])=="AC" || std::string(hrec->vals[id])=="AN") ) { i++; continue; }

			vdict_t *d = (vdict_t*)out_hdr->dict[BCF_DT_ID];
			khint_t k = kh_get(vdict, d, out_hdr->hrec[i]->vals[id]);
			kh_val(d, k).hrec[2] = NULL;
			kh_val(d, k).info[2] |= 0xf;

	        nrm++;
	        out_hdr->nhrec--;
	        if ( i < out_hdr->nhrec ) memmove(&out_hdr->hrec[i],&out_hdr->hrec[i+1],(out_hdr->nhrec-i)*sizeof(bcf_hrec_t*));
	        bcf_hrec_destroy(hrec);
		}
    }
    if ( nrm ) if (bcf_hdr_sync(out_hdr) < 0) vrb.error("Failed to update header");
    */
	nsamples = bcf_hdr_nsamples(out_hdr);

	nswap = {0,0};
	swap_phase = {std::vector<bool>(nsamples, false), std::vector<bool>(nsamples, false)};
	nmatch = std::vector < int > (nsamples, 0);
	nmism = std::vector < int > (nsamples, 0);

	GTa = GTb = NULL;
	mGTa = 0, mGTb=0;

	htsFile * out_fp = hts_open(fname.c_str(),file_format.c_str());
	if ( out_fp == NULL ) vrb.error("Can't write to " + fname + ".");
	if (n_threads > 1) hts_set_opt(out_fp, HTS_OPT_THREAD_POOL, sr->p);
	bcf_hdr_add_sample(out_hdr, NULL);
	if (bcf_hdr_write(out_fp, out_hdr)) vrb.error("Failed to write header to output file");

	std::string fnidx= (file_type == OFILE_BCFC)? fname +".csi" : fname+".tbi";
	if (file_type!=OFILE_VCFU)
	{
		std::ifstream file_idx(fnidx);
		if (file_idx.good())
		{
			file_idx.close();
			if (std::remove(fnidx.c_str()))
			{
				vrb.warning("Detected index file. Were not able to delete it. Trying to skip index creation.");
				file_type=OFILE_VCFU;
			}
		}
		if (file_type!=OFILE_VCFU)
		{
			if (bcf_idx_init(out_fp, out_hdr, 14, fnidx.c_str())) vrb.error("Initializing index");
		}
	}
	bcf1_t* line = NULL;

	int n_variants = 0;
	int n_variants_at_start_cnk = 0;
	int chunk_counter=0;
	int n_sites_buff = 0;

	int prev_readers_size = 0;
    std::string prev_chr = "";
    std::array<int,2> prev_pos = {0,0};
    int first_pos = 0;
    int ifname = 0;

	vrb.bullet("#samples = " + stb.str(nsamples));
	vrb.print("");
	tac.clock();

    // keep only two open files at a time
    while ( ifname < nfiles )
    {
        int new_file = 0;
        while ( sr->nreaders < 2 && ifname < nfiles )
        {
            if ( !bcf_sr_add_reader (sr, filenames[ifname].c_str())) vrb.error("Failed to open " + filenames[ifname] + ".");
            new_file = 1;
            ifname++;
            if ( start_pos[ifname-1]==-1 ) break;   // new chromosome, start with only one file open
            if ( ifname < nfiles && start_pos[ifname]==-1 ) break; // next file starts on a different chromosome
        }

        // is there a line from the previous run? Seek the newly opened reader to that position
        int seek_pos = -1;
        int seek_chr = -1;
        if ( bcf_sr_has_line(sr,0) )
        {
            bcf1_t *line0 = bcf_sr_get_line(sr,0);
            bcf_sr_seek(sr, bcf_seqname(sr->readers[0].header,line0), line0->pos);
            seek_pos = line0->pos;
            seek_chr = bcf_hdr_name2id(out_hdr, bcf_seqname(sr->readers[0].header,line0));
        }
        else if ( new_file ) bcf_sr_seek(sr,NULL,0);  // set to start

        int nret;
        while ( (nret = bcf_sr_next_line(sr)) )
        {
            if ( !bcf_sr_has_line(sr,0) ) if ( bcf_sr_region_done(sr,0) )  bcf_sr_remove_reader(sr, 0);

            // Get a line to learn about current position
            for (i=0; i<sr->nreaders; i++) if ( bcf_sr_has_line(sr,i) ) break;
            line = bcf_sr_get_line(sr,i);

            // This can happen after bcf_sr_seek: indel may start before the coordinate which we seek to.
            if ( seek_chr>=0 && seek_pos>line->pos && seek_chr==bcf_hdr_name2id(out_hdr, bcf_seqname(sr->readers[i].header,line)) ) continue;
            seek_pos = seek_chr = -1;

            //  Check if the position overlaps with the next, yet unopened, reader
            int must_seek = 0;
            while ( ifname < nfiles && start_pos[ifname]!=-1 && line->pos >= start_pos[ifname] )
            {
                must_seek = 1;
                if ( !bcf_sr_add_reader(sr, filenames[ifname].c_str())) vrb.error("Failed to open " + filenames[ifname] + ".");
                if  (sr->nreaders>2) vrb.error("Three files overlapping at position: " + std::to_string(line->pos+1));
                ifname++;
            }
            if ( must_seek )
            {
                bcf_sr_seek(sr, bcf_seqname(sr->readers[i].header,line), line->pos);
                seek_pos = line->pos;
                seek_chr = bcf_hdr_name2id(out_hdr, bcf_seqname(sr->readers[i].header,line));
                continue;
            }

            if ( sr->nreaders>1 && ((!bcf_sr_has_line(sr,0) && !bcf_sr_region_done(sr,0)) || (!bcf_sr_has_line(sr,1) && !bcf_sr_region_done(sr,1))) )
            {
            	const bool uphalf = !bcf_sr_has_line(sr,0);
            	line = bcf_sr_get_line(sr,uphalf);
				write_record(out_fp, out_hdr, sr->readers[0].header, line,uphalf);
                prev_pos[uphalf]=line->pos;
                prev_readers_size = sr->nreaders;
				n_variants++;
            	continue;
            }

            if (sr->nreaders<2)
            {
            	if (prev_readers_size == 0)
            	{
            		n_variants_at_start_cnk = n_variants;
            		prev_chr = std::string(bcf_seqname(sr->readers[i].header,line));
            		first_pos = (int)bcf_sr_get_line(sr,i)->pos + 1;
            		vrb.wait("Cnk " + stb.str(ifname-1) + " [" + prev_chr + ":" + stb.str(first_pos) + "-]");
            	}
            	else if (prev_readers_size == 2)
    			{
            		n_variants_at_start_cnk = n_variants;
            		prev_chr = std::string(bcf_seqname(sr->readers[i].header,line));
            		first_pos = (int)bcf_sr_get_line(sr,i)->pos + 1;
    				vrb.wait("Cnk " + stb.str(ifname-1) + " [" + prev_chr + ":" + stb.str(first_pos) + "-]");
            		//after a buffer we go back to one reader. Chunk 1 is now chunk 0.
    				n_sites_buff = 0;
            		nswap[0]=nswap[1];
            		swap_phase[0] = swap_phase[1];
    			}
				line = bcf_sr_get_line(sr,i);
    			write_record(out_fp, out_hdr, sr->readers[0].header, line,i);
            	prev_pos[i]=line->pos;
            }
            else
            {
            	if (n_sites_buff==0)
            	{
            		prev_chr = std::string(bcf_seqname(sr->readers[i].header,line));
    				vrb.print("Cnk " + stb.str(ifname-2) + " [" + prev_chr + ":" + stb.str(first_pos) + "-" + stb.str(prev_pos[0] + 1) + "] [L=" + stb.str(n_variants-n_variants_at_start_cnk) + "]" );
            		scan_overlap(ifname, bcf_seqname(sr->readers[i].header,line), line->pos);
            	}

				const bool uphalf = n_sites_buff >= nsites_buff_d2.back();
				line = bcf_sr_get_line(sr,uphalf);
				write_record(out_fp, out_hdr, sr->readers[uphalf].header, line, uphalf);
				++n_sites_buff;
	            prev_pos[0]=prev_pos[1]=line->pos;
            }
            prev_readers_size = sr->nreaders;
    		n_variants++;
        }
        if ( sr->nreaders ) while ( sr->nreaders ) bcf_sr_remove_reader(sr, 0);
    }
	vrb.print("Cnk " + stb.str(ifname-1) + " [" + prev_chr + ":" + stb.str(first_pos) + "-" + stb.str(prev_pos[0] + 1) + "] [L=" + stb.str(n_variants-n_variants_at_start_cnk) + "]" );
	if (file_type!=OFILE_VCFU)
	{
		if (bcf_idx_save(out_fp)) vrb.warning("Error writing index");
	}

	if (hts_close(out_fp)) vrb.error("Non zero status when closing VCF/BCF file descriptor");
	out_fp=NULL;
	bcf_sr_destroy(sr);
	sr=NULL;
	bcf_hdr_destroy(out_hdr);
	out_hdr=NULL;
	if (GTa) free(GTa);
	GTa=NULL;
	if (GTb) free(GTb);
	GTb=NULL;
	if (n_variants == 0) vrb.error("No variants to be phased in files");
	vrb.title("Writing completed [L=" + stb.str(n_variants) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

