////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2018 Olivier Delaneau, University of Lausanne
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////
#define _GLOBAL
#include <objects/variant.h>

variant::variant()
{
}

variant::variant(std::string & chr, int bp, std::string & id, std::string & ref, std::string & alt, int idx) {
	this->chr = chr;
	this->bp = bp;
	this->id = id;
	this->ref = ref;
	this->alt = alt;
	this->idx = idx;
	this->cref = 0;
	this->calt = 0;
	this->cmis = 0;
	this->fbin = 0;
	cm = -1;
}

variant::~variant() {
}

unsigned int variant::getMAC() {
	return std::min(cref, calt);
}

bool variant::isSingleton() {
	return (calt == 1 || cref == 1);
}

bool variant::isMonomorphic() {
	return (calt == 0 || cref == 0);
}
