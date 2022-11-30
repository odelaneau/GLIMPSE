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

#define _GLOBAL
#include <objects/variant.h>

variant::variant() : bp(0), type(0), idx(0), cref(0), calt(0), cm(-1), LQ(false)
{
}
variant::variant(const int32_t _bp, const std::string & _id, const std::string & _ref, const std::string & _alt, const uint8_t _type, const int32_t _idx,const uint32_t _cref, const uint32_t _calt, const bool _hq) :
		bp(_bp), id(_id), ref(_ref), alt(_alt), type(_type), idx(_idx), cref(_cref), calt(_calt), cm(-1), LQ(!_hq)
{
}
variant::variant(const int32_t _bp, const char* _id, const char* _ref, const char* _alt, const uint8_t _type, const int32_t _idx, const uint32_t _cref, const uint32_t _calt, const bool _hq) :
		bp(_bp), id(std::string(_id)), ref(std::string(_ref)), alt(std::string(_alt)), type(_type), idx(_idx), cref(_cref), calt(_calt), cm(-1), LQ(!_hq)
{
}


variant::~variant() {
}

unsigned int variant::getMAC() const {
	return std::min(cref, calt);
}

bool variant::isSingleton() const {
	return (calt == 1 || cref == 1);
}

bool variant::isMonomorphic() const {
	return (calt == 0 || cref == 0);
}

std::string variant::toString() const {
	return std::to_string(bp) + "_" + ref + "_" + alt;
}
