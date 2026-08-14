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

#include <containers/hvar_cache.h>
#include <io/retry_io.h>
#include <chrono>

//Builds a phase-only, haplotype-major cache of an existing binary reference panel's
//common-variant matrix (see containers/hvar_cache.h for why: the on-disk .bin's
//row-major layout can't be mmap'd for a memory win, because state selection touches
//scattered haplotype columns across nearly every page of every row regardless of
//panel size). This is a one-time, explicit, opt-in maintenance step -- the source
//.bin is only ever read, never modified or rebuilt from source data, so it can be run
//against reference panels that are too expensive to regenerate.
void caller::build_hvar_cache() {
	vrb.title("[GLIMPSE2] Building hvar cache for GLIMPSE2_phase");

	std::string reference_filename = options["reference"].as < std::string > ();
	std::string cache_path = options.count("hvar-cache-file") ? options["hvar-cache-file"].as < std::string > () : default_hvar_cache_path(reference_filename);

	vrb.bullet("Reference panel      : [" + reference_filename + "]");
	vrb.bullet("Cache file (output)  : [" + cache_path + "]");

	//Force a full, row-major load of H, ignoring any existing cache -- we are
	//(re)building it from the source .bin, not consuming a previous cache.
	vrb.wait("  * Binary reference panel parsing");
	tac.clock();
	retry_with_backoff("reading binary reference panel [" + reference_filename + "]", 3, std::chrono::seconds(1), [&]() -> attempt_result {
		std::string err_msg;
		bool non_retryable = false;
		const bool ok = read_binary_reference_panel(reference_filename, err_msg, non_retryable, /*skip_cache=*/true);
		return { ok, non_retryable, err_msg };
	});
	vrb.bullet("Binary reference panel parsing [done] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	print_ref_panel_info("Binary");

	vrb.bullet("Transposing common-variant matrix to haplotype-major layout for caching...");
	tac.clock();
	bitmatrix HvarRefT;
	HvarRefT.allocate(H.n_ref_haps, H.n_com_sites);
	H.HvarRef.transpose(HvarRefT);
	vrb.bullet("Transpose done (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");

	std::string err_msg;
	if (!write_hvar_cache(cache_path, reference_filename, H.n_ref_haps, H.n_com_sites, HvarRefT, err_msg))
		vrb.error("Failed to write hvar cache: " + err_msg);

	vrb.bullet("Hvar cache written to [" + cache_path + "]");
	vrb.bullet("GLIMPSE2_phase runs against this --reference panel will now automatically use this cache to reduce memory usage.");
}
