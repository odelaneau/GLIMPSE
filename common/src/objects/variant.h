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

#include "otools.h"
#include "checksum_utils.h"
#include "boost/serialization/serialization.hpp"

/**
 * @class variant
 * @brief Represents a genetic variant with position, alleles, and metadata.
 *
 * Contains information about the variant's base pair position, ID, reference and alternate alleles,
 * variant type, index, encoded reference and alternate alleles, genetic map position (cM),
 * and a low-quality flag.
 */
class variant {
public:

	/** Base pair position of the variant */
	int32_t bp;

	/** Variant identifier */
	std::string id;

	/** Reference allele sequence */
	std::string ref;

	/** Alternate allele sequence */
	std::string alt;

	/** Variant type (encoded as uint8_t) */
	uint8_t type;

	/** Index of the variant */
	int32_t idx;

	/** Encoded reference allele */
	uint32_t cref;

	/** Encoded alternate allele */
	uint32_t calt;

	/** Genetic map position in centiMorgans */
	double cm;

	/** Flag indicating if the variant is low quality */
	bool LQ;

	/**
	 * @brief Default constructor initializes variant with default values.
	 */
	variant();

	/**
	 * @brief Constructor with parameters.
	 * @param bp Base pair position
	 * @param id Variant identifier string
	 * @param ref Reference allele string
	 * @param alt Alternate allele string
	 * @param type Variant type
	 * @param idx Variant index
	 * @param cref Encoded reference allele
	 * @param calt Encoded alternate allele
	 * @param hq High quality flag (true if high quality)
	 */
	variant(const int32_t bp, const std::string & id, const std::string & ref, const std::string & alt, const uint8_t type, const int32_t idx, const uint32_t cref, const uint32_t calt, const bool hq);

	/**
	 * @brief Constructor with C-string parameters.
	 * @param bp Base pair position
	 * @param id Variant identifier C-string
	 * @param ref Reference allele C-string
	 * @param alt Alternate allele C-string
	 * @param type Variant type
	 * @param idx Variant index
	 * @param cref Encoded reference allele
	 * @param calt Encoded alternate allele
	 * @param hq High quality flag (true if high quality)
	 */
	variant(const int32_t bp, const char* id, const char* ref, const char* alt, const uint8_t type, const int32_t idx, const uint32_t cref, const uint32_t calt, const bool hq);

	/** Destructor */
	~variant();

	/**
	 * @brief Check if variant is a singleton (occurs once in dataset).
	 * @return True if singleton, else false
	 */
	bool isSingleton() const;

	/**
	 * @brief Check if variant is monomorphic (no alternate alleles).
	 * @return True if monomorphic, else false
	 */
	bool isMonomorphic() const;

	/**
	 * @brief Get minor allele count (MAC) of variant.
	 * @return Minor allele count as unsigned int
	 */
	unsigned int getMAC() const;

	/**
	 * @brief Get string representation of variant.
	 * @return String describing variant information
	 */
	std::string toString() const;

	/** Serialization method for variant */
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

	/**
	 * @brief Update checksum with variant data for data integrity.
	 * @param crc Checksum object
	 */
	void update_checksum(checksum &crc) const
	{
		crc.process_data(bp);
		crc.process_data(id);
		crc.process_data(ref);
		crc.process_data(alt);
		crc.process_data(type);
		crc.process_data(idx);
		crc.process_data(cref);
		crc.process_data(calt);
		crc.process_data(cm);
		crc.process_data(LQ);
	}
};

#endif
