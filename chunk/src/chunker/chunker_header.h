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

#ifndef _CHUNKER_H
#define _CHUNKER_H

#include "otools.h"
#include "gmap_reader.h"

struct chunk_info
{
	std::string chr;
	std::vector<long int> buf_start;
	std::vector<long int> buf_stop;
	std::vector<long int> cnk_start;
	std::vector<long int> cnk_stop;
	std::vector<float> curr_window_cm_size;
	std::vector<long int> curr_window_mb_size;
	std::vector<long int> curr_window_count;
	std::vector<long int> curr_window_common_count;

	chunk_info(){};
	chunk_info(std::string _chr) : chr(_chr){};

	void reset()
	{
		buf_start.clear();
		buf_stop.clear();
		cnk_start.clear();
		cnk_stop.clear();
		curr_window_cm_size.clear();
		curr_window_mb_size.clear();
		curr_window_count.clear();
		curr_window_common_count.clear();
	}

	void add_chunk(long int _buf_start, long int _buf_stop, long int _cnk_start, long int _cnk_stop, float _curr_window_cm_size, long int _curr_window_mb_size, long int _curr_window_count, long int _curr_window_common_count)
	{
		assert(_buf_start < _buf_stop);
		assert(_cnk_start < _cnk_stop);
		assert(_buf_start <= _cnk_start);
		assert(_buf_stop >= _buf_stop);

		buf_start.push_back(_buf_start);
		buf_stop.push_back(_buf_stop);
		cnk_start.push_back(_cnk_start);
		cnk_stop.push_back(_cnk_stop);
		curr_window_cm_size.push_back(_curr_window_cm_size);
		curr_window_mb_size.push_back(_curr_window_mb_size);
		curr_window_count.push_back(_curr_window_count);
		curr_window_common_count.push_back(_curr_window_common_count);
	}

	void output_to_file(output_file& fd)
	{
		for (long int i=0; i<buf_start.size(); ++i)
		{
			if (i<buf_start.size()-1)
			{
				long int mean_curr_stop_next_start = (cnk_stop[i] + cnk_start[i+1]) /2;
				cnk_stop[i] = mean_curr_stop_next_start;
				cnk_start[i + 1] = mean_curr_stop_next_start+1;
			}
			fd << i << "\t" << chr << "\t" << chr << ":"<< buf_start[i] << "-" << buf_stop[i] << "\t" << chr << ":" << cnk_start[i] << "-" << cnk_stop[i] << "\t" << curr_window_cm_size[i] << "\t" << curr_window_mb_size[i] << "\t" << curr_window_count[i] <<"\t" << curr_window_common_count[i] << std::endl;
		}
	}
};

class chunker {
public:
	//COMMAND LINE OPTIONS
	bpo::options_description descriptions;
	bpo::variables_map options;

	//long intERNAL DATA
	//variant_map V;

	float sparse_maf;
	std::string chrID;
	std::vector < float > positions_all_cm;
	std::vector < long int > positions_all_mb;

	std::vector < float > positions_common_cm;
	std::vector < long int > positions_common_mb;

	std::vector < long int > all2common;
	std::vector < long int > common2all;

	std::multimap < long int, long int > map_positions_all;	//associative container of variant with position in bp

	std::vector<float> chunk_cm_length;
	std::vector<long int> chunk_mb_length;
	std::vector<long int> chunk_common_count;

	chunk_info cnk_info;


	//PARAMETERS
	std::string gregion;
	long int start;
	long int stop;
	bool whole_chr;
	long int contig_len;

	float window_cm;
	float window_mb;
	long int window_count;
	float buffer_cm;
	float buffer_mb;
	long int buffer_count;

	//CONSTRUCTOR
	chunker();
	~chunker();

		/**
	 * @brief Reads variant data from a VCF/BCF file into chunker data structures.
	 *
	 * This function opens the specified VCF/BCF file, restricts reading to the given genomic
	 * region, and iterates through all biallelic variants. For each variant, it extracts allele
	 * counts (AC) and allele numbers (AN) from the INFO fields, computes the minor allele frequency (MAF),
	 * and classifies the variant as common or rare based on the @ref sparse_maf threshold.
	 *
	 * Variant positions, along with mapping tables between all variants and common variants,
	 * are stored in the chunker's internal member variables for later use in imputation/phasing.
	 *
	 * If @ref whole_chr is enabled, the function will also retrieve and store the contig length
	 * for the given chromosome from the VCF/BCF header or index. Requires the input file to have
	 * an accompanying index (.tbi/.csi) file.
	 *
	 * @param fmain    Path to the input VCF or BCF file (must be indexed).
	 * @param region   Genomic region string in format "CHR" or "CHR:START-END".
	 * @param nthreads Number of threads to use for reading (if > 1, enables multi-threaded reading).
	 *
	 * @throws std::runtime_error If:
	 *         - The file cannot be opened.
	 *         - The index cannot be loaded.
	 *         - The region cannot be set.
	 *         - The AC or AN INFO fields are missing or invalid.
	 *         - The file type is not recognized as VCF or BCF.
	 *
	 * @note Only biallelic variants (n_allele == 2) are processed. Multiallelic variants are skipped.
	 * @note AC and AN must be present in the INFO field and be single integer values.
	 * @note Positions are stored in 1-based coordinates (VCF format is 0-based).
	 * @note Updates the following member variables:
	 *       - positions_all_mb
	 *       - positions_common_mb
	 *       - map_positions_all
	 *       - all2common
	 *       - common2all
	 *       - contig_len (if whole_chr is true)
	 *       - chrID (chromosome ID of the first common variant)
	 *
	 * @details
	 * **Algorithm workflow:**
	 * 1. **Initialization**
	 *    - Start a timer for performance logging.
	 *    - Initialize a synchronous VCF/BCF reader (`bcf_srs_t`) from htslib.
	 *    - Configure reader to not collapse multi-allelic sites and to require an index.
	 *    - Enable multi-threading if `nthreads > 1`.
	 *
	 * 2. **Set region filter**
	 *    - Restrict reading to the given `region` (chromosome or range).
	 *    - Abort with an error if the region cannot be set.
	 *
	 * 3. **Open input file**
	 *    - Add the file as a reader.
	 *    - If opening or index loading fails, terminate with an error message.
	 *
	 * 4. **Optional: Get contig length** (only if `whole_chr == true`)
	 *    - Detect if the file is VCF or BCF and load the appropriate index (`.tbi` or `.csi`).
	 *    - Retrieve sequence names from the index.
	 *    - Match the `region` chromosome and attempt to read its length from the VCF header.
	 *    - If unavailable, assign a fallback contig length.
	 *
	 * 5. **Prepare buffers**
	 *    - Initialize counters for total and common variants.
	 *    - Allocate arrays for reading AC and AN INFO field values.
	 *
	 * 6. **Read and process variants**
	 *    - Loop over all lines in the file/region:
	 *      - Keep only **biallelic** variants (`n_allele == 2`).
	 *      - For the first common variant, store the chromosome ID (`chrID`).
	 *      - Retrieve AC and AN values and validate them.
	 *      - Compute **minor allele frequency (MAF)**:
	 *        \f[
	 *        MAF = \min\left(\frac{AN - AC}{AN}, \frac{AC}{AN}\right)
	 *        \f]
	 *      - Determine if the variant is **common** (`MAF >= sparse_maf`).
	 *      - Store the variant position (converted from 0-based to 1-based).
	 *      - Update:
	 *        - `positions_all_mb` with all variant positions.
	 *        - `map_positions_all` mapping position → index.
	 *        - `all2common` mapping all variants to common variant indices.
	 *      - If common, also update:
	 *        - `positions_common_mb` with positions of common variants.
	 *        - `common2all` mapping common variant indices back to all variant indices.
	 *
	 * 7. **Cleanup**
	 *    - Destroy the BCF/VCF reader.
	 *    - Free AC/AN buffers.
	 *
	 * 8. **Final reporting**
	 *    - Log total, rare, and common variant counts.
	 *    - Log elapsed processing time.
	 *
	 * @see sparse_maf
	 * @see whole_chr
	 */
	void readData(std::string fmain, std::string region, long int nthreads);

	//PARAMETERS
	void declare_options();
	void parse_command_line(std::vector < std::string > &);
	void check_options();
	void verbose_options();
	void verbose_files();

	//FILE I/O
	void split_recursive(output_file &, long int &, std::string &, long int, long int);
	void split_recursive_no_reset(output_file &, long int &, std::string &, long int, long int);
	void split_sequential(output_file &, long int &, std::string &, long int, long int, const bool);
	void add_buffer(const long int, const long int, long int&, long int&);
	void chunk();
	void chunk(std::vector < std::string > &);
	void buildCoordinates();
	void parseRegion();


	std::vector < long int > getByPos(const long int, const std::multimap < long int, long int >& map_positions);
	void setGeneticMap(const gmap_reader&, const std::multimap < long int, long int >&, const std::vector<long int>& positions_mb, std::vector<float>& positions_cm);
	void setGeneticMap(const std::multimap < long int, long int >&, const std::vector<long int>& positions_mb, std::vector<float>& positions_cm);
	long int setCentiMorgan(const std::vector < long int > & pos_bp, const std::vector < double > & pos_cM, const std::multimap < long int, long int >& map_positions, std::vector<float> & positions_cm);
	long int interpolateCentiMorgan(const std::vector < long int > & pos_bp, const std::vector < double > & pos_cM, const std::multimap < long int, long int >&, const std::vector<long int> & positions_mb, std::vector<float>& positions_cm);
	unsigned long int length() const;
	double lengthcM() const;

};


#endif


