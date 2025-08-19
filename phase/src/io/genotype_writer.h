/**
 * @file genotype_writer.h
 * @brief Declaration of the genotype_writer class for outputting genotype data to various formats.
 * 
 * This header defines the genotype_writer class which is responsible for writing
 * phased genotype data, variant annotations, and associated metadata to
 * genetics file formats such as VCF, BCF, and optionally BGEN.
 * 
 * The class interfaces with core data structures including haplotype sets,
 * genotype likelihoods, variant maps, and variant callers to efficiently
 * output imputed or phased genotype results.
 * 
 * Supported features include:
 * - Writing VCF/BCF files with optional compression and multi-threading.
 * - Handling of genotype posteriors and dosages.
 * - Updating headers from external FASTA index (.fai) files.
 * - Optional BGEN format support with configurable compression.
 * 
 * Dependencies:
 * - htslib for VCF/BCF reading and writing.
 * - genfile::bgen for BGEN output (conditional).
 */

#ifndef _GENOTYPE_WRITER_H
#define _GENOTYPE_WRITER_H

#include "otools.h"

#include "variant_map.h"
#include "containers/haplotype_set.h"
#include "containers/genotype_set.h"
#include "caller/caller_header.h"

#ifdef __BGEN__
#include "genfile/bgen/bgen.hpp"
#endif

/**
 * @class genotype_writer
 * @brief Handles output of genotype data to various file formats (VCF, BCF, BGEN).
 * 
 * This class is responsible for writing phased genotype data, variant information,
 * and associated metadata into common genetics file formats. It interfaces with
 * internal data containers such as haplotype sets, genotype sets, variant maps,
 * and caller information.
 * 
 * The writer supports compression options and multi-threaded output to optimize performance.
 * 
 * @note This class depends on external libraries such as htslib for VCF/BCF handling and
 * optionally genfile::bgen for BGEN file format support.
 * 
 * @section Data Members
 * - `H`: Reference to the haplotype_set containing phased haplotype data.
 * - `G`: Reference to the genotype_set containing genotype likelihood and posterior data.
 * - `V`: Reference to the variant_map containing variant position and allele information.
 * - `C`: Reference to a caller object responsible for variant calling data and parameters.
 * 
 * @section Constructors and Destructor
 * - `genotype_writer(const haplotype_set &, const genotype_set &, const variant_map &, const caller & C)`: Initializes the writer with references to data containers.
 * - `~genotype_writer()`: Destructor to clean up resources.
 * 
 * @section Member Functions
 * - `writeGenotypes(...)`: Writes genotype data to file with specified format, compression, and threading.
 * - `update_header_from_fai(...)`: Updates VCF header with contig info from a FASTA index (.fai) file.
 */
class genotype_writer {
public:
	//DATA
	const haplotype_set & H;   /**< Reference to haplotype data */
	const genotype_set & G;    /**< Reference to genotype data */
	const variant_map & V;     /**< Reference to variant metadata */
	const caller & C;          /**< Reference to variant caller */

	//CONSTRUCTORS/DESCTRUCTORS
	/**
	 * @brief Constructor initializes writer with references to core data objects.
	 * @param haplotype_set Reference to haplotype data container.
	 * @param genotype_set Reference to genotype data container.
	 * @param variant_map Reference to variant metadata.
	 * @param caller Reference to variant caller instance.
	 */
	genotype_writer(const haplotype_set &, const genotype_set &, const variant_map &, const caller & C);
	
	/**
	 * @brief Destructor handles cleanup.
	 */
	~genotype_writer();

	/**
	 * @brief Write genotypes and associated info to file.
	 * @param fname Output filename.
	 * @param oformat Output format (VCF, BCF, etc).
	 * @param ocompr Compression option (none, zlib, zstd).
	 * @param n_bits_bgen Number of bits used for BGEN format (if applicable).
	 * @param n_main Number of main samples or relevant parameter.
	 * @param n_threads Number of threads to use for writing.
	 * @param fai_fname Filename of fasta index (.fai) to update contig info.
	 */
	void writeGenotypes(const std::string fname, OutputFormat oformat, OutputCompression ocompr, const int n_bits_bgen, const int n_main, const int n_threads, const std::string fai_fname) const;
	
	/**
	 * @brief Update VCF header contig information from a fasta index (.fai).
	 * @param hdr Pointer to VCF header structure.
	 * @param fai_fname Path to fasta index file.
	 */
	void update_header_from_fai(bcf_hdr_t * hdr, const std::string fai_fname) const;

#ifdef __BGEN__
	/**
	 * @brief Write genotype data to a VCF/BCF output file.
	 * 
	 * This function writes genotype data stored in memory to a specified file
	 * in VCF or BCF format, optionally compressed and threaded.
	 * 
	 * The output includes variant records with genotype calls, dosages, and genotype posteriors.
	 * 
	 * @param[in] fname        Output file name.
	 * @param[in] output_fmt   Format of output (VCF, BCF).
	 * @param[in] output_compr Compression type for output.
	 * @param[in] n_bits_bgen  Number of bits for BGEN output (not used in this code).
	 * @param[in] n_main       Number of main samples.
	 * @param[in] n_threads    Number of threads to use.
	 * @param[in] fai_fname    Optional reference fasta index filename to update header.
	 * 
	 * @details
	 * The function performs the following major steps:
	 * 
	 * @mermaid
	 * flowchart TD
	 *     A[Start: Call writeGenotypes()] --> B[Open output file (VCF/BCF) with htslib]
	 *     B --> C[Create and initialize VCF header]
	 *     C --> D[Add sample names to header]
	 *     D --> E[Write header to file]
	 *     E --> F[Allocate memory buffers for genotypes, dosages, posteriors]
	 *     F --> G[For each variant site in V.vec_pos:]
	 *     G --> H[Clear previous VCF record]
	 *     H --> I[Update variant info: position, IDs, alleles]
	 *     I --> J[For each sample:]
	 *     J --> K[Check if genotype data exists for sample at site]
	 *     K --> L{If genotype exists}
	 *     L --> M[Extract genotype probabilities (GP0, GP1, GP2) and dosage (DS)]
	 *     M --> N[Determine phased genotype calls]
	 *     N --> O[Fill genotype buffer with phased GT values]
	 *     L -->|No| P[Fill genotype buffer with default/reference genotypes]
	 *     O --> Q[Round dosage and posteriors, store in buffers]
	 *     P --> Q
	 *     Q --> R[Update INFO fields (RAF, AF, INFO) in VCF record]
	 *     R --> S[Update FORMAT fields (GT, DS, GP) in VCF record]
	 *     S --> T[Write VCF record to output file]
	 *     T --> U[Increment progress bar]
	 *     U --> G
	 *     G --> V[Free allocated buffers]
	 *     V --> W[Destroy VCF record and header]
	 *     W --> X[Close output file]
	 *     X --> Y[Build index for output file]
	 *     Y --> Z[Finish and report success]
	 * @endmermaid
	 * 
	 * @note The function uses htslib for handling VCF/BCF input/output and indexing.
	 */
	void writeGenotypesBgen(const std::string fname, OutputFormat oformat, OutputCompression ocompr, const int n_bits_bgen, const int n_main, const int n_threads) const;
	
	/**
	 * @brief Initialize a BGEN context structure for genotype data output.
	 * 
	 * This function prepares the BGEN context before writing genotype data by:
	 * - Setting the total number of variants (sites) to be written.
	 * - Setting the total number of samples (individuals) in the dataset.
	 * - Configuring the compression algorithm and data layout flags.
	 * 
	 * @param[out] context        Reference to the BGEN context object to initialize.
	 * @param[in]  ocompr         Enum indicating the output compression method to use.
	 * @param[in]  n_sites_unbuf  Number of variant sites that will be output (unbuffered).
	 * 
	 * @details
	 * BGEN is a binary format for storing large-scale genotype data efficiently. To write
	 * data correctly, the context must specify how many variants and samples exist, as well
	 * as compression and layout options. This function:
	 * 
	 * 1. Sets `context.number_of_variants` to the provided count of variants to write.
	 * 2. Sets `context.number_of_samples` to the number of samples in the genotype set (`G.vecG.size()`).
	 * 3. Chooses the compression method based on the `OutputCompression` enum:
	 *    - `NONE`: no compression.
	 *    - `ZLIB`: zlib compression.
	 *    - `ZSTD`: Zstandard compression.
	 * 4. Sets `context.flags` by combining a fixed layout flag (`e_Layout2`) with the chosen compression flag.
	 * 
	 * @note
	 * The flag `e_Layout2` specifies the BGEN v1.2 layout, which supports various compression schemes
	 * and is widely compatible with BGEN tools.
	 * 
	 * @warning
	 * Ensure that the number of variants and samples match exactly the data to be written, 
	 * or downstream readers may fail to parse the BGEN file correctly.
	 * 
	 * @mermaid
	 * flowchart TD
	 *     A[Start initialiseBGEN] --> B[Set number_of_variants to n_sites_unbuf]
	 *     B --> C[Set number_of_samples to G.vecG.size()]
	 *     C --> D{Compression type from ocompr}
	 *     D -->|NONE| E[e_NoCompression]
	 *     D -->|ZLIB| F[e_ZlibCompression]
	 *     D -->|ZSTD| G[e_ZstdCompression]
	 *     E --> H[Set flags = e_Layout2 | compression]
	 *     F --> H
	 *     G --> H
	 *     H --> I[End initialiseBGEN]
	 */
	void initialiseBGEN(genfile::bgen::Context& context,OutputCompression ocompr, const size_t n_sites_unbuf) const;
	
	/**
	 * @brief Initialize the BGEN context stream for writing sample identifiers.
	 * 
	 * This function configures the BGEN context to include sample identifiers in the output
	 * stream. It modifies the context flags to enable writing sample names and then
	 * initializes the sample names in the provided output stream.
	 * 
	 * @param[in,out] context  Reference to the BGEN context object to configure.
	 * @param[in,out] oStream  Output stream where the BGEN data, including sample IDs, will be written.
	 * 
	 * @details
	 * BGEN files optionally include sample identifiers to map genotype data to individual samples.
	 * To include these identifiers:
	 * 
	 * 1. The function sets the `e_SampleIdentifiers` flag in the context's flags field.
	 * 2. It calls a helper function `set_sample_names_impl` that writes the sample names
	 *    into the output stream according to the BGEN specification.
	 * 
	 * This setup is necessary before writing genotype data if sample identifiers are
	 * required by downstream tools or workflows.
	 * 
	 * @note
	 * The method assumes the BGEN context and output stream are properly initialized before this call.
	 * 
	 * @warning
	 * Ensure that `set_sample_names_impl` is implemented to correctly write sample IDs,
	 * or the resulting BGEN file may be invalid or incompatible.
	 * 
	 * @mermaid
	 * flowchart TD
	 *     A[Start initialiseStream] --> B[Add e_SampleIdentifiers flag to context.flags]
	 *     B --> C[Call set_sample_names_impl(context, oStream)]
	 *     C --> D[Sample names written to stream]
	 *     D --> E[End initialiseStream]
	 */
	void initialiseStream(genfile::bgen::Context& context,std::ostream& oStream) const ;

	/**
	 * @brief Write sample names into the BGEN output stream and update context offsets.
	 * 
	 * This method writes the sample identifiers (sample names) into the BGEN file stream
	 * according to the BGEN specification, and updates the context's header information
	 * to reflect the new offsets after writing the sample block.
	 * 
	 * @param[in,out] context  Reference to the BGEN context, which maintains metadata including header size and offsets.
	 * @param[in,out] oStream  Output stream where the BGEN file content is being written.
	 * 
	 * @details
	 * The function performs the following operations:
	 * 1. Extracts sample names from the genotype set (`G.vecG`) and stores them in a vector.
	 * 2. Retrieves the current header size from the BGEN context to determine where to write.
	 * 3. Updates internal offsets in the BGEN header to prepare for writing the sample block.
	 * 4. If the `e_SampleIdentifiers` flag is set in the context, writes the sample identifier block
	 *    to the output stream using the genfile BGEN library function `write_sample_identifier_block`.
	 * 5. Updates the header offsets again after writing the sample block.
	 * 6. Moves the output stream write position to just after the header plus 4 bytes.
	 * 
	 * This process ensures that the sample names are correctly embedded in the BGEN file,
	 * which is essential for downstream tools to associate genotype data with samples.
	 * 
	 * @note
	 * The function relies on the genfile BGEN library's ability to write sample blocks and manage offsets.
	 * 
	 * @warning
	 * Proper synchronization of offsets is critical; failing to update offsets correctly
	 * will result in corrupted or unreadable BGEN files.
	 * 
	 * @mermaid
	 * flowchart TD
	 *     A[Start set_sample_names_impl] --> B[Extract sample names from G.vecG]
	 *     B --> C[Get current header size from context]
	 *     C --> D[Update offsets before writing sample block]
	 *     D --> E{Check if sample identifiers flag is set}
	 *     E -->|Yes| F[Write sample identifier block to oStream]
	 *     E -->|No| G[Skip writing sample identifiers]
	 *     F --> H[Update offsets after writing sample block]
	 *     G --> H
	 *     H --> I[Seek stream position to offset + 4]
	 *     I --> J[End set_sample_names_impl]
	 */
	void set_sample_names_impl(genfile::bgen::Context& context,std::ostream& oStream) const;
	
	/**
	 * @brief Update the offset and header block in the BGEN output stream.
	 * 
	 * This function updates the header portion of the BGEN file with the current offset value,
	 * ensuring that the header metadata remains consistent as the file is being written.
	 * 
	 * @param[in,out] context  Reference to the BGEN context containing header information.
	 * @param[in,out] oStream  Output stream where the BGEN file content is being written.
	 * @param[in] offset       The current offset in the file, indicating where data blocks start.
	 * 
	 * @details
	 * The function performs the following steps:
	 * 1. Moves the output stream's write position to the beginning of the file (offset 0).
	 * 2. Checks if the stream is in a good state (no errors).
	 * 3. Writes the updated offset value to the stream using `genfile::bgen::write_offset`.
	 * 4. Writes the updated header block, reflecting the current state of the context, using `genfile::bgen::write_header_block`.
	 * 
	 * Keeping the header and offset updated is essential for file integrity, so that downstream
	 * readers can correctly locate and interpret the genotype data blocks.
	 * 
	 * @note
	 * The BGEN file format requires accurate header offsets for parsing and random access.
	 * 
	 * @warning
	 * If the output stream is in a bad state (e.g., due to previous write failures), the function
	 * will silently skip writing to avoid further corruption.
	 * 
	 * @mermaid
	 * flowchart TD
	 *     A[Start update_offset_and_header_block] --> B[Seek stream to start (pos 0)]
	 *     B --> C{Check if stream is good}
	 *     C -->|Yes| D[Write current offset to stream]
	 *     D --> E[Write updated header block to stream]
	 *     C -->|No| F[Skip writing due to stream error]
	 *     E --> G[End update_offset_and_header_block]
	 *     F --> G
	 */
	void update_offset_and_header_block(genfile::bgen::Context& context,std::ostream& oStream, uint32_t offset) const;
#endif
};

#endif
