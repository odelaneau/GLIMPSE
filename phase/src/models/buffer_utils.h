/*******************************************************************************
 * Copyright (C) 2026 Fulcrum Genomics LLC
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

#ifndef _BUFFER_UTILS_H
#define _BUFFER_UTILS_H

#include <algorithm>
#include <cstddef>
#include <vector>

/**
 * Resizes a scratch buffer whose elements the caller always writes before reading them.
 *
 * Unlike std::vector::resize, growing past the capacity releases the old storage before allocating
 * the new and copies nothing, so the two allocations never coexist. The HMM buffers this is used for
 * reach several GB, and that transient copy would otherwise set the peak memory of the whole run.
 * Capacity still grows geometrically, as std::vector's does, so sizes that creep up between
 * iterations rarely reallocate. Contents after the call are unspecified.
 */
template <typename T, typename Allocator>
void resize_discarding_contents(std::vector<T, Allocator> & buffer, const std::size_t size) {
	if (size > buffer.capacity()) {
		const std::size_t capacity = std::max(2 * buffer.capacity(), size);
		std::vector<T, Allocator>().swap(buffer);
		buffer.reserve(capacity);
	}
	buffer.resize(size);
}

#endif
