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

#include <pir/pir_header.h>

typedef struct {     										// auxiliary data structure
	samFile * fp;											// the file handle
	bam_hdr_t * hdr;										// the file header
	hts_itr_t * iter;										// NULL if a region not specified
	int mq;													// mapping quality threshold
	bool keep_orphan,check_orientation,check_proper_pair;	// parameters for filtering
	int fflag;												// Filter flag
	unsigned int mtlen;										// Maximum observed fragment length
} aux_t;


void pir::extract() {
	readPositions (options["input"].as < string > (), options["region"].as < string > ());
	readSequences(options["bam"].as < string > (), options["region"].as < string > ());
}

void pir::readPositions (string fref, string region) {
	tac.clock();
	vrb.title("Reading variable positions in VCF/BCF file [" + fref + "]");

	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + region + "] in [" + fref + "]");

	//Opening files
	if (!(bcf_sr_add_reader (sr, fref.c_str()))) {
		if (sr->errnum != idx_load_failed) vrb.error("Failed to open file: " + fref + "");
		else vrb.error("Failed to load index of the file: " + fref + "");
	}

	//Parsing file
	int nset;
	bcf1_t * line;
	while ((nset = bcf_sr_next_line (sr))) {
		line =  bcf_sr_get_line(sr, 0);
		if (line->n_allele == 2) {
			bcf_unpack(line, BCF_UN_STR);
			positions.push_back(line->pos);
			refA.push_back(string(line->d.allele[0]));
			altA.push_back(string(line->d.allele[1]));
		}
	}
	bcf_sr_destroy(sr);

	vrb.bullet("Region = " + region + " / #variants = " + stb.str(positions.size()));
	vrb.bullet("Elapsed time = " + stb.str(tac.rel_time()*1.0/1000, 2) + "seconds");
}


void pir::readSequences(string fbam, string region) {
	tac.clock();

	vrb.title("Piling up sequencing reads from BAM file [" + fbam + "]");
	aux_t * data = (aux_t *) malloc (sizeof(aux_t));
	//vrb.bullet("Allocation successful");

	//Opening BAM file
	data->fp = sam_open(fbam.c_str(), "r");
	if (data->fp == 0) vrb.error("Cannot open BAM file!");
	//vrb.bullet("File opened successfully");

	//Reading header
    data->hdr = sam_hdr_read(data->fp);
    if (data->hdr == 0) vrb.error("Cannot parse BAM header!");
    //vrb.bullet("Header read successfully");

    //Loading index file
    hts_idx_t *idx = sam_index_load(data->fp, fbam.c_str());
    if (idx == NULL) vrb.error("Cannot load BAM index!");
    //vrb.bullet("Index file loaded successfully");

    //Setting up filters
    data->fflag = (BAM_FUNMAP | BAM_FSECONDARY);
    data->fflag |= BAM_FQCFAIL;
    data->fflag |= BAM_FDUP;
	data->fflag  |= BAM_FSUPPLEMENTARY;
	data->mq = min_mapq;
	data->keep_orphan = 0;
	data->check_orientation = 1;
	data->check_proper_pair = 1;
	data->mtlen = 0;
	//vrb.bullet("Read filters set up successfully");

	//Seek to region of interest
	data->iter = sam_itr_querys(idx, data->hdr, region.c_str());
	if (data->iter == NULL) vrb.error("Problem jumping to region [" + region + "] in BAM file");
	//else vrb.bullet("Jumping to region [" + region + "] done");

	//Parse BAM file
	parseBam((void *) data);
	vrb.bullet("Largest fragment observed = " + stb.str(data->mtlen) + " bp");
	vrb.bullet("Total number of reads retained = " + stb.str(dataPIR.size()));


	//Filtering
	int longest_pir = 0;
	vector < vector < allele > > TMP;
	for (vector < vector < allele > > :: iterator it = dataPIR.begin(); it != dataPIR.end() ; ++it) {
		longest_pir = max((int)it->size(), longest_pir);
		if (it->size() > 1) TMP.push_back(*it);
	}
	dataPIR = TMP;
	vrb.bullet("Total number of PIRs retained = " + stb.str(dataPIR.size()));
	vrb.bullet("Longest PIRs = " + stb.str(longest_pir));

	//Finalize
    bam_hdr_destroy(data->hdr);
    hts_idx_destroy(idx);
    if (data->fp) sam_close(data->fp);
    hts_itr_destroy(data->iter);
    free(data);
	vrb.bullet("Elapsed time = " + stb.str(tac.rel_time()*1.0/1000, 2) + " seconds");
}

static int read_bam(void *data, bam1_t *b) {
	aux_t * aux = (aux_t*) data;
	int ret, skip = 0;
	do{
		ret = (aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b));
		if (ret < 0) break;
		else if (b->core.tid < 0 || (b->core.flag & aux->fflag) || b->core.qual < aux->mq) {skip = 1; continue;} // unmapped
		else if(b->core.flag & BAM_FPAIRED){ //paired read
			if (b->core.flag & BAM_FMUNMAP){ //mate unmapped
				if (!aux->keep_orphan) { skip = 1;continue;}
			}else{
				if (aux->check_orientation && ((b->core.tid != b->core.mtid) ||  //must be on the same chr
					((b->core.flag & BAM_FREVERSE) && (b->core.flag & BAM_FMREVERSE)) || // read 1 and 2 cannot be on the same strand
					((b->core.pos < b->core.mpos) && (b->core.flag & BAM_FREVERSE)) || //read 1 must be on the +ve strand
					((b->core.pos > b->core.mpos) && !(b->core.flag & BAM_FREVERSE)) //read 2 must be on the -ve strand
				   )) { skip = 1; continue;}
				if (aux->check_proper_pair && !(b->core.flag & BAM_FPROPER_PAIR)) { skip = 1; continue;} //not properly paired
			}
		}
		unsigned int isize = abs(b->core.isize); // fragment length
		if(isize > aux->mtlen) aux->mtlen = isize;
		skip = 0;
	} while (skip);
	return ret;
}

void pir::parseBam(void * d){
	tac.clock();
	aux_t * data = (aux_t *) d;

	//Pile up reads
	const bam_pileup1_t * v_plp;
	int n_plp = 0, tid, pos;
	bam_plp_t s_plp = bam_plp_init(read_bam, (void*)data);
	bam_plp_set_maxcnt(s_plp, 1000);		//Maximum of 1000 reads per position


	//Getting start and end of region
	int beg = data->iter->beg;
	int end = data->iter->end;
	vrb.bullet("Window = [" + stb.str(beg) + " bp -> " + stb.str(end) + " bp]");


	//Looping through variable positions
	snp_index = 0;
	unsigned long basecount = 0, varcount = 0;
	while (((v_plp = bam_plp_auto(s_plp, &tid, &pos, &n_plp)) != 0) && (snp_index < positions.size())) {

		//Seek to correct variable position
		basecount ++;
		//int curr_pos = positions[snp_index]-1;	//0-based
		int curr_pos = positions[snp_index];	//0-based
		while ((pos > curr_pos) && (snp_index < (positions.size() - 1))) {
			snp_index++;
			//curr_pos = positions[snp_index]-1;	//0-based
			curr_pos = positions[snp_index];	//0-based
		}
		if (pos != curr_pos) continue;

		//Looping through reads at that variable position
		varcount ++;
		unsigned int b_ref = 0, b_alt = 0, b_dis = 0;
		for (int iread = 0 ; iread < n_plp ; iread ++) {
			const bam_pileup1_t * p = v_plp + iread;

			//Check if indel
			if (p->is_del) { continue; };
			if (p->is_refskip) { continue; };

			//Check base quality
			int baseq = p->qpos < p->b->core.l_qseq ? bam_get_qual(p->b)[p->qpos] : 0;
			//if (illumina13) baseq = baseq > 31 ? baseq - 31 : 0;
			//cout << baseq << " " << min_baseq << endl;
			if (baseq < min_baseq) { continue; };

			//Check base matching to site
			char base = getBase(bam_seqi(bam_get_seq(p->b), p->qpos));
			//cout << base << " vs [" << refA[snp_index][0] << " " << altA[snp_index][0] << "]" << endl;
			if (base != refA[snp_index][0] && base != altA[snp_index][0]) {
				continue;
			};

			//Everything is okay if we reach this point, so let's add this into PIR
			//Find/Add read in PIR collection
			string read_name = string(bam_get_qname(p->b));
			int read_index = 0;
			unordered_map < string, int > :: iterator it = mapPIR.find(read_name);
			if (it == mapPIR.end()) {
				read_index = dataPIR.size();
				mapPIR.insert(pair < string, int > (read_name, read_index));
				dataPIR.push_back(vector < allele >());
			} else read_index = it->second;

			//cout << "Push ! " << read_index << " " << snp_index << " " << baseq << " " << base << endl;
			dataPIR[read_index].push_back(allele(snp_index, baseq, base==altA[snp_index][0]));
		}
	}

	bam_plp_reset(s_plp);
	bam_plp_destroy(s_plp);
	vrb.bullet(stb.str(basecount) + " bases read");
	vrb.bullet(stb.str(varcount) + " variable positions read");
}
