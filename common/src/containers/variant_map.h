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

#ifndef _SNP_SET_H
#define _SNP_SET_H

#include <utils/otools.h>
#include <objects/variant.h>
#include <io/gmap_reader.h>
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/serialization/map.hpp"



class variant_map {
public :
	//DATA
	//GENOMIC REGION
	std::string chrid;
	int input_start;
	int input_stop;
	int output_start;
	int output_stop;
	std::string input_gregion;
	std::string output_gregion;

	std::vector < variant * > vec_pos;			//vector of variants ordered by position in bp
	std::multimap < int, variant * > map_pos;	//associative container of variant with position in bp

	//CONSTRUCTOR/DESTRUCTOR
	variant_map();
	~variant_map();

	//METHODS
	std::size_t size() const;
	std::vector < variant * > getByPos(const int);
	std::vector < variant * > getByRef(const int, const std::string &, const std::string &);
	variant * getByIndex(const int);
	void push(variant *);
	void setGeneticMap(const gmap_reader&);
	void setGeneticMap();
	int setCentiMorgan(const std::vector < long int > & pos_bp, const std::vector < double > & pos_cM);
	int interpolateCentiMorgan(const std::vector < long int > & pos_bp, const std::vector < double > & pos_cM);
	unsigned int length() const;
	double lengthcM() const;

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & chrid;
		ar & input_start;
		ar & input_stop;
		ar & output_start;
		ar & output_stop;
		ar & input_gregion;
		ar & output_gregion;
		ar & vec_pos;
		ar & map_pos;
	}
};

#endif
