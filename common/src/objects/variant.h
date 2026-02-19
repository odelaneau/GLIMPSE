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

#ifndef _VARIANT_H
#define _VARIANT_H

#include <utils/otools.h>
#include "boost/serialization/serialization.hpp"

class variant {
public :

	int32_t bp;
	std::string id;
	std::string ref;
	std::string alt;
	uint8_t type;
	int32_t idx;
	uint32_t cref;
	uint32_t calt;
	double cm;
	bool LQ;		//Low quality site

	//CONSTRUCTOR/DESTRUCTOR
	variant();
	variant(const int32_t bp, const std::string & id, const std::string & ref, const std::string & alt, const uint8_t type, const int32_t idx, const uint32_t cref, const uint32_t calt, const bool hq);
	variant(const int32_t bp, const char* id, const char* ref, const char* alt, const uint8_t type, const int32_t idx, const uint32_t cref, const uint32_t calt, const bool hq);

	~variant();

	bool isSingleton() const;
	bool isMonomorphic() const;
	unsigned int getMAC() const;
	std::string toString() const;

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & bp;
		ar & id;
		ar & ref;
		ar & alt;
		ar & type;
		ar & idx;
		ar & cref;
		ar & calt;
		ar & cm;
		ar & LQ;
	}

};

#endif
