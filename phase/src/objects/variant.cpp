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

#define _GLOBAL
#include <objects/variant.h>

variant::variant(const string & _chr,const int _bp, const string & _id, const string & _ref, const string & _alt,const int _idx,const int _cref, const int _calt) :
	chr(_chr), bp(_bp), id(_id), ref(_ref), alt(_alt), idx(_idx), cref(_cref), calt(_calt), cm(-1)
{
}

variant::~variant() {
}

unsigned int variant::getMAC() const {
	return min(cref, calt);
}

bool variant::isSingleton() const {
	return (calt == 1 || cref == 1);
}

bool variant::isMonomorphic() const {
	return (calt == 0 || cref == 0);
}

string variant::toString() const {
	return chr + ":" + std::to_string(bp) + "_" + ref + "_" + alt;
}
