/// \file
/// \brief The amino_acid class header

/// \copyright
/// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 2011, Orengo Group, University College London
///
/// This program is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef AMINO_ACID_H_INCLUDED
#define AMINO_ACID_H_INCLUDED

#include <boost/operators.hpp>
#include <boost/optional.hpp>
#include <boost/throw_exception.hpp>

#include "common/type_aliases.h"
#include "exception/invalid_argument_exception.h"
#include "file/pdb/pdb_record.h"
#include "structure/structure_type_aliases.h"

#include <iosfwd>
// #include <map>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cath {

	/// \brief TODOCUMENT
	class amino_acid final : private boost::equivalent      < amino_acid,
	                                 boost::totally_ordered < amino_acid > > {
	private:
		static const char_str_str_tpl_vec & LETTER_CODE_AND_NAME_LIST();

		using char_size_unordered_map   = std::unordered_map<char,        size_t>;
		using string_size_unordered_map = std::unordered_map<std::string, size_t>;

		static const char_size_unordered_map   & INDEX_OF_LETTER();
		static const string_size_unordered_map & INDEX_OF_CODE();
		static const string_size_unordered_map & INDEX_OF_NAME();

		/// \brief TODOCUMENT
		str_opt raw_string;

		/// \brief TODOCUMENT
		size_opt index;

		template <typename T, size_t I> static std::unordered_map<T, size_t> build_index_unordered_map();
		template <typename T, size_t I> static T get_label(const size_t &);

		size_t get_letter_index(const char &);
		void set_letter_code_or_name(const std::string &);

		void check_is_proper_amino_acid() const;

	public:
		amino_acid(const std::string &,
		           const file::pdb_record &);
		explicit amino_acid(const std::string &);
		explicit amino_acid(const char &);

		bool is_proper_amino_acid() const;

		char        get_letter() const;
		std::string get_code() const;
		std::string get_name() const;

//		static const std::string UNKNOWN_AMINO_ACID_NAME;
	};

	amino_acid_vec make_amino_acids_of_chars(const char_vec &);

	/// \brief TODOCUMENT
	///
	/// \todo Since, some part of this arose as taking excessive & non-trivial time in some profile
	///       it is worth seeing if it can be made faster, probably by using a std::array and constexpr_find()
	///       and constexpr_for_n() to generate lookups with minimal runtime overhead.
	inline const char_str_str_tpl_vec & amino_acid::LETTER_CODE_AND_NAME_LIST() {
		static const char_str_str_tpl_vec letter_code_and_name_list = {
			std::make_tuple( 'A', "ALA", "Alanine"                            ),
			std::make_tuple( 'B', "ASX", "Ambiguous Asparagine/Aspartic Acid" ), // eg PDBs 156b, 1kp0, 1pgk, 2atc, 2fmd, 2rxn, 3atc, 3bcl and 3e2o
			std::make_tuple( 'C', "CYS", "Cysteine"                           ),
			std::make_tuple( 'D', "ASP", "Aspartic Acid"                      ),
			std::make_tuple( 'E', "GLU", "Glutamic Acid"                      ),
			std::make_tuple( 'F', "PHE", "Phenylalanine"                      ),
			std::make_tuple( 'G', "GLY", "Glycine"                            ),
			std::make_tuple( 'H', "HIS", "Histidine"                          ),
			std::make_tuple( 'I', "ILE", "Isoleucine"                         ),
			std::make_tuple( 'J', "XLE", "Leucine/Isoleucine"                 ), // eg PDBs?
			std::make_tuple( 'K', "LYS", "Lysine"                             ),
			std::make_tuple( 'L', "LEU", "Leucine"                            ),
			std::make_tuple( 'M', "MET", "Methionine"                         ),
			std::make_tuple( 'N', "ASN", "Asparagine"                         ),
			std::make_tuple( 'O', "PYL", "Pyrrolysine"                        ),
			std::make_tuple( 'P', "PRO", "Proline"                            ),
			std::make_tuple( 'Q', "GLN", "Glutamine"                          ),
			std::make_tuple( 'R', "ARG", "Arginine"                           ),
			std::make_tuple( 'S', "SER", "Serine"                             ),
			std::make_tuple( 'T', "THR", "Threonine"                          ),
			std::make_tuple( 'U', "SEC", "Selenocysteine"                     ), // eg PDBs 1aa6, 1cc1, 1fdi, 1fdo, 1h0h, 1kqf, 1kqg, 1pae, 1pfp, 2bc7, 2bc8, 2iv2, 2wpn, 2xsk, 3ean, 3eao, 3fwf, 3fwi, 3fwj, 3u5s, 3ze7, 3ze8, 3ze9, 3zea, 4kl8, 4kn9, 4ko1, 4ko2, 4ko3, 4ko4
			std::make_tuple( 'V', "VAL", "Valine"                             ),
			std::make_tuple( 'W', "TRP", "Tryptophan"                         ),
			std::make_tuple( 'X', "UNK", "Unknown"                            ),
			std::make_tuple( 'Y', "TYR", "Tyrosine"                           ),
			std::make_tuple( 'Z', "GLX", "Ambiguous Glutamine/Glutamic Acid"  )  // eg PDBs 156b, 1kp0, 1pgk, 2rxn, 2tnc, 3bcl and 4cpa
		};
		return letter_code_and_name_list;
	}

	/// \brief TODOCUMENT
	template <typename T, size_t I>
	inline std::unordered_map<T, size_t> amino_acid::build_index_unordered_map() {
		std::unordered_map<T, size_t> index_map;
		for (size_t amino_acid_ctr = 0; amino_acid_ctr < LETTER_CODE_AND_NAME_LIST().size(); ++amino_acid_ctr) {
			const T &amino_acid_label = std::get<I>( LETTER_CODE_AND_NAME_LIST()[ amino_acid_ctr ] );
			index_map.insert( std::make_pair( amino_acid_label, amino_acid_ctr ) );
		}
		return index_map;
	}

	/// \brief TODOCUMENT
	template <typename T, size_t I>
	inline T amino_acid::get_label(const size_t &arg_index ///< TODOCUMENT
	                               ) {
		if (arg_index >= amino_acid::LETTER_CODE_AND_NAME_LIST().size()) {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Amino acid index is out of range"));
		}
		return std::get<I>( amino_acid::LETTER_CODE_AND_NAME_LIST()[ arg_index ] );
	}

	/// \brief TODOCUMENT
	inline size_t amino_acid::get_letter_index(const char &arg_letter ///< The 1 letter code
	                                           ) {
		const auto index_itr = INDEX_OF_LETTER().find( arg_letter );
		if ( index_itr == common::cend( INDEX_OF_LETTER() ) ) {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception(
				"Amino acid letter \"" + std::string{ 1, arg_letter } + "\" is not a recognised letter (currently case-sensitive)"
			));
		}
		return index_itr->second;
	}

	/// \brief Ctor for amino_acid
	inline amino_acid::amino_acid(const std::string      &arg_string,    ///< TODOCUMENT
	                              const file::pdb_record &arg_pdb_record ///< TODOCUMENT
	                              ) {
		switch ( arg_pdb_record ) {
			case ( file::pdb_record::ATOM   ) : {
				set_letter_code_or_name( arg_string );
				break;
			}
			case ( file::pdb_record::HETATM ) : {
				if ( arg_string.length() != 3 ) {
					BOOST_THROW_EXCEPTION(common::invalid_argument_exception(
						"Cannot create a HETATM amino acid from a string that is not 3 characters long"
					));
				}
				raw_string = arg_string;
				break;
			}
			default : {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception(
					"Value of arg_pdb_record not recognised whilst constructing an amino_acid"
				));
			}
		}
		assert(   raw_string ||   index );
		assert( ! raw_string || ! index );
	}

	/// \brief Ctor for amino_acid
	inline amino_acid::amino_acid(const char &arg_letter ///< The 1 letter code
	                              ) : index{ get_letter_index( arg_letter ) } {
	}

	/// \brief TODOCUMENT
	inline bool amino_acid::is_proper_amino_acid() const {
		return static_cast<bool>( index );
	}

	/// \brief TODOCUMENT
	inline char amino_acid::get_letter() const {
		check_is_proper_amino_acid();
		return get_label<char, 0>( *index );
	}

	/// \brief TODOCUMENT
	inline std::string amino_acid::get_code() const {
		return is_proper_amino_acid() ? get_label<std::string, 1>( *index )
		                              : *raw_string;
	}

	/// \brief TODOCUMENT
	inline std::string amino_acid::get_name() const {
		check_is_proper_amino_acid();
		return get_label<std::string, 2>( *index );
	}

	/// \brief Simple less-than operator for amino-acid
	///
	/// This can be helpful for ordering amino acids (eg as a key to a map in sequence_similarity_score).
	///
	/// This is extended to the full range of <=, ==, !=, >, >= with Boost Operators (equivalent and totally_ordered).
	///
	/// \relates amino_acid
	inline bool operator<(const amino_acid &arg_amino_acid_1, ///< The first amino acid to compare
	                      const amino_acid &arg_amino_acid_2  ///< The second amino acid to compare
	                      ) {
		return ( arg_amino_acid_1.get_letter() < arg_amino_acid_2.get_letter() );
	}

	std::string get_code_of_amino_acid_letter(const char &);

	char get_letter_of_amino_acid_code(const std::string &);

	std::istream & operator>>(std::istream &,
	                          amino_acid &);
}

#endif
