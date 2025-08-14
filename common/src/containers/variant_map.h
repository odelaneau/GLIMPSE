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

#include "otools.h"
#include "checksum_utils.h"
#include "variant.h"
#include "gmap_reader.h"
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/serialization/map.hpp"

/**
 * @class variant_map
 * @brief Container for managing a collection of genomic variants with positional indexing.
 *
 * This class stores variants ordered by base pair position in a vector and a multimap,
 * allowing efficient access by position and reference allele.
 * It also tracks genomic region boundaries and supports genetic map interpolation.
 */
class variant_map {
public :
	/// Chromosome identifier (e.g., "chr1")
	std::string chrid;

	/// Start position of the input region (base pairs)
	int input_start;

	/// End position of the input region (base pairs)
	int input_stop;

	/// Start position of the output region (base pairs)
	int output_start;

	/// End position of the output region (base pairs)
	int output_stop;

	/// Input genomic region name or descriptor
	std::string input_gregion;

	/// Output genomic region name or descriptor
	std::string output_gregion;

	/// Vector of pointers to variants ordered by base pair position
	std::vector < variant * > vec_pos;

	/// Multimap associating base pair positions to variant pointers (allows multiple variants per position)
	std::multimap < int, variant * > map_pos;

	/// Constructor
	variant_map();

	/// Destructor (frees variant pointers if applicable)
	~variant_map();

	/**
	 * @brief Get the total number of variants in the map.
	 * @return Number of variants.
	 */
	std::size_t size() const;

	/**
	 * @brief Retrieve all variants at a specific base pair position.
	 * @param pos Base pair position.
	 * @return Vector of pointers to variants at position pos.
	 */
	std::vector < variant * > getByPos(const int pos);

	/**
	 * @brief Retrieve variants by position and reference/alternate alleles.
	 * @param pos Base pair position.
	 * @param ref Reference allele string.
	 * @param alt Alternate allele string.
	 * @return Vector of pointers to variants matching the criteria.
	 */
	std::vector < variant * > getByRef(const int pos, const std::string & ref, const std::string & alt);

	/**
	 * @brief Retrieve variant by index in the vector.
	 * @param index Index into vec_pos.
	 * @return Pointer to variant at index.
	 */
	variant * getByIndex(const int index);

	/**
	 * @brief Add a variant pointer to the container.
	 * 
	 * The variant is added to both vec_pos and map_pos.
	 * 
	 * @param var Pointer to variant to add.
	 */
	void push(variant * var);

	/**
	 * @brief Set genetic map positions using a genetic map reader.
	 * @param gmap Genetic map reader object.
	 */
	void setGeneticMap(const gmap_reader& gmap);

	/**
	 * @brief Set genetic map positions without parameters.
	 * May initialize or reset genetic map data.
	 */
	void setGeneticMap();

	/**
	 * @brief Compute centiMorgan positions from base pair positions and known genetic map.
	 * 
	 * @param pos_bp Vector of base pair positions.
	 * @param pos_cM Vector of genetic map positions in centiMorgans.
	 * @return Number of variants successfully assigned centiMorgan positions.
	 */
	int setCentiMorgan(const std::vector < long int > & pos_bp, const std::vector < double > & pos_cM);

	/**
	 * @brief Interpolate centiMorgan positions for variants given base pair and cM vectors.
	 * 
	 * @param pos_bp Vector of base pair positions.
	 * @param pos_cM Vector of centiMorgan positions.
	 * @return Number of variants with interpolated centiMorgan positions.
	 */
	int interpolateCentiMorgan(const std::vector < long int > & pos_bp, const std::vector < double > & pos_cM);

	/**
	 * @brief Get the length of the variant map in base pairs (output_stop - output_start).
	 * @return Length in base pairs.
	 */
	unsigned int length() const;

	/**
	 * @brief Get the length of the variant map in centiMorgans.
	 * @return Length in cM.
	 */
	double lengthcM() const;

	/// Boost Serialization access
	friend class boost::serialization::access;

	/**
	 * @brief Serialize or deserialize the variant_map object.
	 * @tparam Archive Archive type.
	 * @param ar Archive instance.
	 * @param version Serialization version.
	 */
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

	/**
	 * @brief Update a checksum object with the data in the variant_map.
	 * 
	 * Processes string fields and each variant's checksum.
	 * 
	 * @param crc Checksum object to update.
	 */
	void update_checksum(checksum &crc) const
	{
		crc.process_data(chrid);
		crc.process_data(input_start);
		crc.process_data(input_stop);
		crc.process_data(output_start);
		crc.process_data(output_stop);
		crc.process_data(input_gregion);
		crc.process_data(output_gregion);
		for (variant * var : vec_pos) {
			var->update_checksum(crc);
		}
		for (const auto pair : map_pos) {
			crc.process_data(pair.first);
			pair.second->update_checksum(crc);
		}
	}
};

#endif
