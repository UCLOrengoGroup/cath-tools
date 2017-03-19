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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_PROTEIN_AMINO_ACID_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_PROTEIN_AMINO_ACID_H

#include <boost/operators.hpp>
#include <boost/optional.hpp>
#include <boost/range/irange.hpp>
#include <boost/throw_exception.hpp>
#include <boost/utility/string_ref.hpp>
#include <boost/variant.hpp>

#include "common/algorithm/contains.hpp"
#include "common/char_arr_type_aliases.hpp"
#include "common/type_aliases.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "exception/out_of_range_exception.hpp"
#include "file/pdb/pdb_record.hpp"
#include "structure/protein/dna_atom.hpp"
#include "structure/structure_type_aliases.hpp"

#include <iosfwd>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cath {
	namespace detail { struct aa_code_getter; }
	namespace detail { struct aa_type_getter; }

	/// \brief Represent the type of amino_acid record
	enum class amino_acid_type {
		AA,      ///< A standard amino_acid amino acid record (!)
		HETATOM, ///< A HETATM amino acid record
		DNA      ///< A DNA/RNA amino acid record
	};

	/// \brief TODOCUMENT
	///
	/// \todo Make get_code() return char_3_arr and make LETTER_CODE_AND_NAME_LIST store char_3_arr
	class amino_acid final : private boost::equivalent      < amino_acid,
	                                 boost::totally_ordered < amino_acid > > {
	private:
		friend detail::aa_code_getter;
		friend detail::aa_type_getter;

		/// \brief The number of chars used to store a HETATM
		static constexpr size_t NUM_HETATM_CHARS = 3;

		static const char_str_str_tpl_vec & LETTER_CODE_AND_NAME_LIST();

		using char_size_unordered_map   = std::unordered_map<char,        uint>;
		using string_size_unordered_map = std::unordered_map<std::string, uint>;

		static const char_size_unordered_map   & INDEX_OF_LETTER();
		static const string_size_unordered_map & INDEX_OF_CODE();
		static const string_size_unordered_map & INDEX_OF_NAME();

		using aa_variant_t     = uint;
		using hetatm_variant_t = std::array<char, NUM_HETATM_CHARS>;
		using dna_variant_t    = dna_atom;

		/// \brief The data for the amino acid, one of three different types in the variant
		boost::variant<
			aa_variant_t,
			hetatm_variant_t,
			dna_variant_t
		> data;

		template <typename T, size_t I> static std::unordered_map<T, uint> build_index_unordered_map();
		template <typename T, size_t I> static T get_label(const size_t &);

		static uint get_letter_index(const char &);
		void set_letter_code_or_name(const std::string &);

		const aa_variant_t & check_is_proper_amino_acid() const;

	public:
		amino_acid(const std::string &,
		           const file::pdb_record &);
		explicit amino_acid(const std::string &);
		explicit amino_acid(const char &);

		amino_acid_type get_type() const;

		const char_3_arr & get_hetatm_chars() const;

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
	inline std::unordered_map<T, uint> amino_acid::build_index_unordered_map() {
		std::unordered_map<T, uint> index_map;
		for (const uint &amino_acid_ctr : boost::irange( 0u, static_cast<uint>( LETTER_CODE_AND_NAME_LIST().size() ) ) ) {
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
	inline uint amino_acid::get_letter_index(const char &arg_letter ///< The 1 letter code
	                                         ) {
		const auto index_itr = INDEX_OF_LETTER().find( arg_letter );
		if ( index_itr == common::cend( INDEX_OF_LETTER() ) ) {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception(
				"Amino acid letter \"" + std::string{ arg_letter } + "\" is not a recognised letter (currently case-sensitive)"
			));
		}
		return index_itr->second;
	}

	/// \brief Check that this is a proper amino acid (rather than a HETATM record or DNA/RNA pseudo-amino-acid)
	///        and throw an exception if not, otherwise return the index value
	inline auto amino_acid::check_is_proper_amino_acid() const -> const aa_variant_t & {
		const aa_variant_t * const aa_ptr = boost::get<aa_variant_t>( &data );
		if ( aa_ptr == nullptr ) {
			BOOST_THROW_EXCEPTION(common::out_of_range_exception("Cannot use a generic HETATM amino_acid or DNA value as a proper, ATOM-record amino acid"));
		}
		return *aa_ptr;
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
				if ( common::contains( INDEX_OF_CODE(), arg_string ) ) {
					set_letter_code_or_name( arg_string );
				}
				else {
					data = char_3_arr{ {
						arg_string[ 0 ],
						arg_string[ 1 ],
						arg_string[ 2 ]
					} };
				}
				break;
			}
			default : {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception(
					"Value of arg_pdb_record not recognised whilst constructing an amino_acid"
				));
			}
		}
	}

	/// \brief Ctor for amino_acid
	inline amino_acid::amino_acid(const char &arg_letter ///< The 1 letter code
	                              ) : data{ get_letter_index( arg_letter ) } {
	}

	namespace detail {

		/// \brief A visitor to get the amino_acid_type from the variant in amino_acid
		struct aa_type_getter : public boost::static_visitor<amino_acid_type> {
			amino_acid_type operator()(const amino_acid::aa_variant_t     &) const {
				return amino_acid_type::AA;
			}
			amino_acid_type operator()(const amino_acid::hetatm_variant_t &) const {
				return amino_acid_type::HETATOM;
			}
			amino_acid_type operator()(const amino_acid::dna_variant_t    &) const {
				return amino_acid_type::DNA;
			}
		};

		/// \brief A visitor to get the three-letter code from the variant in amino_acid
		struct aa_code_getter : public boost::static_visitor<std::string> {
			std::string operator()(const amino_acid::aa_variant_t     &x) const {
				return amino_acid::get_label<std::string, 1>( x );
			}
			std::string operator()(const amino_acid::hetatm_variant_t &x) const {
				return { common::cbegin( x ), common::cend( x ) };
			}
			std::string operator()(const amino_acid::dna_variant_t    &x) const {
				return to_three_char_str( x );
			}
		};
	} // namespace detail

	/// \brief Get the amino_acid_type of this amino_acid
	inline amino_acid_type amino_acid::get_type() const {
		return boost::apply_visitor( detail::aa_type_getter{}, data );
	}

	/// \brief Getter for the raw string for HETATM amino_acids
	inline const char_3_arr & amino_acid::get_hetatm_chars() const {
		const hetatm_variant_t * const raw_string_ptr = boost::get<hetatm_variant_t>( &data );
		if ( raw_string_ptr == nullptr ) {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Unable to get raw_string from non-HETATM amino_acid"));
		}
		return *raw_string_ptr;
	}

	/// \brief TODOCUMENT
	inline char amino_acid::get_letter() const {
		return get_label<char, 0>( check_is_proper_amino_acid() );
	}

	/// \brief TODOCUMENT
	inline std::string amino_acid::get_code() const {
		return boost::apply_visitor( detail::aa_code_getter{}, data );
	}

	/// \brief TODOCUMENT
	inline std::string amino_acid::get_name() const {
		return get_label<std::string, 2>( check_is_proper_amino_acid() );
	}

	namespace detail {

		/// \brief Make a less-than comparator for the specified amino_acid
		inline std::tuple<uint8_t, char, char_3_arr> make_amino_acid_lt_comparator(const amino_acid &arg_aa ///< The amino_acid to query
		                                                                           ) {
			/// \todo Come GCC 6 as oldest compiler, remove this type alias and just return using braces
			using uint8_char_char_3_arr_tpl = std::tuple<uint8_t, char, char_3_arr>;
			switch ( arg_aa.get_type() ) {
				case ( amino_acid_type::AA ) : {
					return uint8_char_char_3_arr_tpl{ 0, arg_aa.get_letter(), char_3_arr{ { 0, 0, 0 } } };
				}
				case ( amino_acid_type::HETATOM ) : {
					return uint8_char_char_3_arr_tpl{ 1, 0,                   arg_aa.get_hetatm_chars() };
				}
				case ( amino_acid_type::DNA ) : {
					return uint8_char_char_3_arr_tpl{ 2, 0,                   char_3_arr{ { 0, 0, 0 } } };
				}
				default : {
					BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Value of arg_aa.get_type() not recognised whilst in make_amino_acid_lt_comparator()"));
					return {}; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
				}
			}
		}
	} // namespace detail

	/// \brief Simple less-than operator for amino-acid
	///
	/// This can be helpful for ordering amino acids (eg as a key to a map in sequence_similarity_score).
	///
	/// This evaluates:
	///  * proper amino acids as less than HETATM amino acids
	///  * proper amino acids as less than each other based on their single letters
	///  * non-proper amino acids as less than each other based on their raw strings
	///
	/// This is extended to the full range of <=, ==, !=, >, >= with Boost Operators (equivalent and totally_ordered).
	///
	/// \relates amino_acid
	inline bool operator<(const amino_acid &arg_amino_acid_1, ///< The first  amino acid to compare
	                      const amino_acid &arg_amino_acid_2  ///< The second amino acid to compare
	                      ) {
		return (
			detail::make_amino_acid_lt_comparator( arg_amino_acid_1 )
			<
			detail::make_amino_acid_lt_comparator( arg_amino_acid_2 )
		);
	}

	/// \brief Return whether the specified amino_acid is a proper amino_acid (rather than a HETATM record or DNA/RNA pseudo-amino-acid)
	inline bool is_proper_amino_acid(const amino_acid &arg_amino_acid ///< The amino_acid to query
	                                 ) {
		return ( arg_amino_acid.get_type() == amino_acid_type::AA );
	}

	std::string get_code_of_amino_acid_letter(const char &);

	char get_letter_of_amino_acid_code(const std::string &);

	std::ostream & operator<<(std::ostream &,
	                          const amino_acid &);

	std::istream & operator>>(std::istream &,
	                          amino_acid &);
} // namespace cath

#endif
