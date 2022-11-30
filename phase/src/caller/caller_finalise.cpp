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

#include <caller/caller_header.h>

#include <io/genotype_writer.h>

void caller::write_files_and_finalise() {
	vrb.title("Finalization:");

	//step0: multi-threading
	if (options["threads"].as < int > () > 1) pthread_mutex_destroy(&mutex_workers);

	std::string out_file = options["output"].as < std::string > ();
	//step1: writing best guess haplotypes in VCF/BCF file
	if (output_fmt != OutputFormat::BGEN)
		genotype_writer(H, G, V, *this).writeGenotypes(options["output"].as < std::string > (), output_fmt, output_compr, 0, options["main"].as < int > (), options["threads"].as < int > (), options.count("contigs-fai")?options["contigs-fai"].as <std::string> () : "");
	else
	{
		#ifdef __BGEN__
			genotype_writer(H, G, V, *this).writeGenotypesBgen(options["output"].as < std::string > (), output_fmt, output_compr, bgen_bits, options["main"].as < int > (), options["threads"].as < int > ());
		#else
			//should have already thrown an error at the beginning. But ok let's not waste computation in the we arrive here
			vrb.warning("Output in bgen format but the program seems to not support BGEN output. Switching to BCF (extension .bcf will added to output file)");
			genotype_writer(H, G, V, *this).writeGenotypes(options["output"].as < std::string > () + ".bcf", OutputFormat::BCF, OutputCompression::ZLIB, 0, options["main"].as < int > (), options["threads"].as < int > (),options.count("contigs-fai")?options["contigs-fai"].as <std::string> () : "");
			out_file += ".bcf";
		#endif
	}
	vrb.print("");
	vrb.bullet("Output file [" + out_file + "]");

	//step2: Measure overall running time
	vrb.bullet("Total running time [" + stb.str(tac.abs_time_display()) + "] (" + stb.str(tac.abs_time()) + "s)");
}
