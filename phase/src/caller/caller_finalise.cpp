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
	if (A.mNumThreads > 1) pthread_mutex_destroy(&mutex_workers);

	//step1: writing best guess haplotypes in VCF/BCF file
	if (A.mOutputFormat != OutputFormat::BGEN)
		genotype_writer(H, G, V, *this).writeGenotypes(A);
	else
	{
		#ifdef __BGEN__
			genotype_writer(H, G, V, *this).writeGenotypesBgen(A.mOutputFilename, A.mOutputFormat, A.mOutputCompression, A.mBgenNBits,  A.mMain, A.mNumThreads);
		#else
			//should have already thrown an error at the beginning. But ok let's not waste computation in the we arrive here
			vrb.error("Output in BGEN format but the program seems to not support BGEN output.");
			//genotype_writer(H, G, V, *this).writeGenotypes(A.mOutputFilename + ".bcf", OutputFormat::BCF, OutputCompression::ZLIB, 0, A.mMain, A.mNumThreads, A.mContigFaiFilename);
		#endif
	}
	vrb.print("");
	verbose_output_file();

	//step2: Measure overall running time
	vrb.title("Total running time [" + stb.str(tac.abs_time_display()) + "] (" + stb.str(tac.abs_time()) + "s)");
}
