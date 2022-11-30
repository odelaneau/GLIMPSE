/*******************************************************************************
 * Copyright (C) 2021 Rick Wertenbroek, University of Lausanne (UNIL),
 * University of Applied Sciences and Arts Western Switzerland (HES-SO),
 * School of Management and Engineering Vaud (HEIG-VD).
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

/**
 * @file c_api.h
 * @author Rick Wertenbroek
 * @date 23 Sep 2021
 * @brief C API for xSqueezeIt
 *
 */
#ifndef __C_API_H__
#define __C_API_H__

#include "htslib/vcf.h"
#include "htslib/synced_bcf_reader.h"

typedef void* c_xcf;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Allocates a new C compatible Xcf class
 *
 */
c_xcf *c_xcf_new();

// /**
//  * @brief Adds a reader to the given Xcf class
//  *
//  */
// void c_xcf_add_reader(c_xcf *x, bcf_sr_t* reader, int reader_id);

/**
 * @brief Adds all readers to given Xcf class
 * 
 */
void c_xcf_add_readers(c_xcf *x, bcf_srs_t* readers);

/**
 * @brief Update Xcf class with readers
 * 
 */
void c_xcf_update_readers(c_xcf *x, bcf_srs_t* readers);

/**
 * @brief Get the sample name
 * 
 */
const char* c_xcf_sample_name(c_xcf *x, int reader_id, const bcf_hdr_t *hdr, int sample_id);

/**
 * @brief Get the number of samples from file
 * 
 */
int c_xcf_nsamples(const char* fname);

/**
 * @brief equivalent with bcf_get_genotypes but also compatible with xSqueezeIt format
 *        This function will check if the given reader (id) is VCF/BCF or xSqueezeIt in
 *        the Xcf class and call the appropriate methods accordingly
 * 
 */
#define c_xcf_get_genotypes(x,reader_id,hdr,line,dst,ndst) __c__xcf__get__genotypes__void(x,reader_id,hdr,line,(void**)(dst),ndst)
int __c__xcf__get__genotypes__void(c_xcf *x, int reader_id, const bcf_hdr_t *hdr, bcf1_t *line, void **dst, int *ndst);

/**
 * @brief Deallocates the given Xcf class
 *
 */
void c_xcf_delete(c_xcf *x);

#ifdef __cplusplus
}
#endif

#endif /* __C_API_H__ */
