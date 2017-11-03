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
#include <boost/throw_exception.hpp>
#include <boost/utility/string_ref.hpp>
#include <boost/variant.hpp>

#include "common/algorithm/contains.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/char_arr_type_aliases.hpp"
#include "common/cpp17/invoke.hpp"
#include "common/function/ident.hpp"
#include "common/optional/make_optional_if.hpp"
#include "common/string/char_arr_to_string.hpp"
#include "common/type_aliases.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/exception/out_of_range_exception.hpp"
#include "file/pdb/pdb_record.hpp"
#include "structure/protein/dna_atom.hpp"
#include "structure/structure_type_aliases.hpp"

#include <iosfwd>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cath {
	namespace detail { struct aa_code_getter;            }
	namespace detail { struct tolerant_aa_letter_getter; }
	namespace detail { struct aa_type_getter;            }

	using char_char_3_arr_str_tpl     = std::tuple<char, char_3_arr, std::string>;
	using char_char_3_arr_str_tpl_vec = std::vector<char_char_3_arr_str_tpl>;

	/// \brief Represent the type of amino_acid record
	enum class amino_acid_type : char {
		AA,      ///< A standard amino_acid amino acid record (!)
		HETATOM, ///< A HETATM amino acid record
		DNA      ///< A DNA/RNA amino acid record
	};

	std::string to_string(const amino_acid_type &);

	/// \brief TODOCUMENT
	///
	/// \todo Make get_code() return char_3_arr and make LETTER_CODE_AND_NAME_LIST store char_3_arr
	class amino_acid final : private boost::equivalent      < amino_acid,
	                                 boost::totally_ordered < amino_acid > > {
	private:
		friend detail::aa_code_getter;
		friend detail::aa_type_getter;
		friend detail::tolerant_aa_letter_getter;

		/// \brief The number of chars used to store a HETATM
		static constexpr size_t NUM_HETATM_CHARS = 3;

		static const char_char_3_arr_str_tpl_vec & LETTER_CODE_AND_NAME_LIST();

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

		template <typename T, size_t I, typename Proj = common::ident>
		static std::unordered_map<T, uint> build_index_unordered_map(Proj && = Proj{});
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

		char_opt    get_letter_if_amino_acid() const;
		char        get_letter_tolerantly() const;
		char_3_arr  get_code() const;
		std::string get_name() const;

//		static const std::string UNKNOWN_AMINO_ACID_NAME;
	};

	amino_acid_vec make_amino_acids_of_chars(const char_vec &);

	/// \brief TODOCUMENT
	///
	/// \todo Since, some part of this arose as taking excessive & non-trivial time in some profile
	///       it is worth seeing if it can be made faster, probably by using a std::array and constexpr_find()
	///       and constexpr_for_n() to generate lookups with minimal runtime overhead.
	inline const char_char_3_arr_str_tpl_vec & amino_acid::LETTER_CODE_AND_NAME_LIST() {
		static const char_char_3_arr_str_tpl_vec letter_code_and_name_list = {
			std::make_tuple( 'A', char_3_arr{{ 'A','L','A' }}, "Alanine"                            ),
			std::make_tuple( 'B', char_3_arr{{ 'A','S','X' }}, "Ambiguous Asparagine/Aspartic Acid" ), // eg PDBs 156b, 1kp0, 1pgk, 2atc, 2fmd, 2rxn, 3atc, 3bcl and 3e2o
			std::make_tuple( 'C', char_3_arr{{ 'C','Y','S' }}, "Cysteine"                           ),
			std::make_tuple( 'D', char_3_arr{{ 'A','S','P' }}, "Aspartic Acid"                      ),
			std::make_tuple( 'E', char_3_arr{{ 'G','L','U' }}, "Glutamic Acid"                      ),
			std::make_tuple( 'F', char_3_arr{{ 'P','H','E' }}, "Phenylalanine"                      ),
			std::make_tuple( 'G', char_3_arr{{ 'G','L','Y' }}, "Glycine"                            ),
			std::make_tuple( 'H', char_3_arr{{ 'H','I','S' }}, "Histidine"                          ),
			std::make_tuple( 'I', char_3_arr{{ 'I','L','E' }}, "Isoleucine"                         ),
			std::make_tuple( 'J', char_3_arr{{ 'X','L','E' }}, "Leucine/Isoleucine"                 ), // eg PDBs?
			std::make_tuple( 'K', char_3_arr{{ 'L','Y','S' }}, "Lysine"                             ),
			std::make_tuple( 'L', char_3_arr{{ 'L','E','U' }}, "Leucine"                            ),
			std::make_tuple( 'M', char_3_arr{{ 'M','E','T' }}, "Methionine"                         ),
			std::make_tuple( 'N', char_3_arr{{ 'A','S','N' }}, "Asparagine"                         ),
			std::make_tuple( 'O', char_3_arr{{ 'P','Y','L' }}, "Pyrrolysine"                        ),
			std::make_tuple( 'P', char_3_arr{{ 'P','R','O' }}, "Proline"                            ),
			std::make_tuple( 'Q', char_3_arr{{ 'G','L','N' }}, "Glutamine"                          ),
			std::make_tuple( 'R', char_3_arr{{ 'A','R','G' }}, "Arginine"                           ),
			std::make_tuple( 'S', char_3_arr{{ 'S','E','R' }}, "Serine"                             ),
			std::make_tuple( 'T', char_3_arr{{ 'T','H','R' }}, "Threonine"                          ),
			std::make_tuple( 'U', char_3_arr{{ 'S','E','C' }}, "Selenocysteine"                     ), // eg PDBs 1aa6, 1cc1, 1fdi, 1fdo, 1h0h, 1kqf, 1kqg, 1pae, 1pfp, 2bc7, 2bc8, 2iv2, 2wpn, 2xsk, 3ean, 3eao, 3fwf, 3fwi, 3fwj, 3u5s, 3ze7, 3ze8, 3ze9, 3zea, 4kl8, 4kn9, 4ko1, 4ko2, 4ko3, 4ko4
			std::make_tuple( 'V', char_3_arr{{ 'V','A','L' }}, "Valine"                             ),
			std::make_tuple( 'W', char_3_arr{{ 'T','R','P' }}, "Tryptophan"                         ),
			std::make_tuple( 'X', char_3_arr{{ 'U','N','K' }}, "Unknown"                            ),
			std::make_tuple( 'Y', char_3_arr{{ 'T','Y','R' }}, "Tyrosine"                           ),
			std::make_tuple( 'Z', char_3_arr{{ 'G','L','X' }}, "Ambiguous Glutamine/Glutamic Acid"  )  // eg PDBs 156b, 1kp0, 1pgk, 2rxn, 2tnc, 3bcl and 4cpa
		};
		return letter_code_and_name_list;
	}

	/// \brief TODOCUMENT
	template <typename T, size_t I, typename Proj>
	inline std::unordered_map<T, uint> amino_acid::build_index_unordered_map(Proj &&arg_proj ///< An optional projection function to apply to each of the elements
	                                                                         ) {
		std::unordered_map<T, uint> index_map;
		for (const uint &amino_acid_ctr : common::indices( static_cast<uint>( LETTER_CODE_AND_NAME_LIST().size() ) ) ) {
			index_map.emplace(
				common::invoke(
					std::forward<Proj>( arg_proj ),
					std::get<I>( LETTER_CODE_AND_NAME_LIST()[ amino_acid_ctr ] )
				),
				amino_acid_ctr
			);
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
			const auto code = get_code();
			BOOST_THROW_EXCEPTION(common::out_of_range_exception(
				R"(Cannot use a generic HETATM amino_acid or DNA value as a proper, ATOM-record amino acid. Problem was a )"
				+ to_string( get_type() )
				+ R"( with code ")"
				+ std::string( code.begin(), code.end() )
				+ R"(".)"
			));
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
				return;
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
				return;
			}
		}
		BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Value of arg_pdb_record not recognised whilst constructing an amino_acid"));
	}

	/// \brief Ctor for amino_acid
	inline amino_acid::amino_acid(const char &arg_letter ///< The 1 letter code
	                              ) : data{ get_letter_index( arg_letter ) } {
	}

	namespace detail {

		/// \brief A visitor to get the three-letter code from the variant in amino_acid
		struct aa_code_getter : public boost::static_visitor<char_3_arr> {
			char_3_arr operator()(const amino_acid::aa_variant_t     &x) const {
				return amino_acid::get_label<char_3_arr, 1>( x );
			}
			char_3_arr operator()(const amino_acid::hetatm_variant_t &x) const {
				return x;
			}
			char_3_arr operator()(const amino_acid::dna_variant_t    &x) const {
				return to_three_char_arr( x );
			}
		};

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

		/// \brief A visitor to get the one-letter code from the variant in amino_acid.
		///        This just returns 'X' for DNA/HETATM amino_acids.
		struct tolerant_aa_letter_getter : public boost::static_visitor<char> {
			char operator()(const amino_acid::aa_variant_t     &x) const {
				return amino_acid::get_label<char, 0>( x );
			}
			char operator()(const amino_acid::hetatm_variant_t &x) const {
				// Some of the most common codes, as extracted from PDB dir with command like:
				//
				//     ls -1 | xargs grep -hPB99999 '^TER   ' | grep -P '^HETATM' | awk '{print substr( $0, 18, 3 )}' | sort | uniq -c | sort -g
				//
				// Codes from pages like : http://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/MSE
				//
				// Default to 'X'
				return
					( x == char_3_arr{ { ' ', 'I', 'C' } } ) ? 'C' :
					( x == char_3_arr{ { ' ', 'I', 'G' } } ) ? 'G' :
					( x == char_3_arr{ { 'H', 'Y', 'P' } } ) ? 'P' :
					( x == char_3_arr{ { 'L', 'C', 'G' } } ) ? 'G' :
					( x == char_3_arr{ { 'M', 'L', 'Y' } } ) ? 'K' :
					( x == char_3_arr{ { 'M', 'S', 'E' } } ) ? 'M' :
					( x == char_3_arr{ { 'N', 'C', 'X' } } ) ? 'N' :
					( x == char_3_arr{ { 'O', 'M', 'G' } } ) ? 'G' :
					( x == char_3_arr{ { 'P', 'C', 'A' } } ) ? 'E' :
					( x == char_3_arr{ { 'P', 'S', 'U' } } ) ? 'U' :
					( x == char_3_arr{ { 'P', 'T', 'R' } } ) ? 'Y' :
					( x == char_3_arr{ { 'S', 'E', 'P' } } ) ? 'S' :
					                                           'X' ;
			}
			char operator()(const amino_acid::dna_variant_t    & ) const {
				return 'X';
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

	/// \brief Get the single letter for the amino_acid if it is an amino acid, or none otherwise
	inline char_opt amino_acid::get_letter_if_amino_acid() const {
		const aa_variant_t * const aa_ptr = boost::get<aa_variant_t>( &data );
		return common::make_optional_if_fn(
			aa_ptr != nullptr,
			[&] { return get_label<char, 0>( *aa_ptr ); }
		);
	}

	/// \brief Get a single letter for the amino_acid, which is the expected letter for
	///        the standard amino acids and a few HETATMs (eg MSE -> M) and is 'X' otherwise
	inline char amino_acid::get_letter_tolerantly() const {
		return boost::apply_visitor( detail::tolerant_aa_letter_getter{}, data );
	}

	/// \brief TODOCUMENT
	inline char_3_arr amino_acid::get_code() const {
		return boost::apply_visitor( detail::aa_code_getter{}, data );
	}

	/// \brief TODOCUMENT
	inline std::string amino_acid::get_name() const {
		return get_label<std::string, 2>( check_is_proper_amino_acid() );
	}

	/// \brief Get a string of the three-letter-code for the specified amino_acid
	///
	/// \relates amino_acid
	inline std::string get_code_string(const amino_acid &arg_amino_acid ///< The amino_acid to query
	                                   ) {
		return common::char_arr_to_string( arg_amino_acid.get_code() );
	}

	namespace detail {

		/// \brief Make a less-than comparator for the specified amino_acid
		inline std::tuple<uint8_t, char, char_3_arr> make_amino_acid_lt_comparator(const amino_acid &arg_aa ///< The amino_acid to query
		                                                                           ) {
			/// \todo Come GCC 6 as oldest compiler, remove this type alias and just return using braces
			using uint8_char_char_3_arr_tpl = std::tuple<uint8_t, char, char_3_arr>;
			switch ( arg_aa.get_type() ) {
				case ( amino_acid_type::AA ) : {
					return uint8_char_char_3_arr_tpl{ 0, *arg_aa.get_letter_if_amino_acid(), char_3_arr{ { 0, 0, 0 } } };
				}
				case ( amino_acid_type::HETATOM ) : {
					return uint8_char_char_3_arr_tpl{ 1, 0,                                  arg_aa.get_hetatm_chars() };
				}
				case ( amino_acid_type::DNA ) : {
					return uint8_char_char_3_arr_tpl{ 2, 0,                                  char_3_arr{ { 0, 0, 0 } } };
				}
			}
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Value of arg_aa.get_type() not recognised whilst in make_amino_acid_lt_comparator()"));
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
	///
	/// \relates amino_acid
	inline bool is_proper_amino_acid(const amino_acid &arg_amino_acid ///< The amino_acid to query
	                                 ) {
		return ( arg_amino_acid.get_type() == amino_acid_type::AA );
	}

	/// \brief Return whether the specified amino_acid is for water (HOH)
	///
	/// \relates amino_acid
	inline bool is_water(const amino_acid &arg_amino_acid ///< The amino_acid to query
	                     ) {
		return ( arg_amino_acid.get_code() == char_3_arr{{ 'H', 'O', 'H' }} );
	}

	char_3_arr get_code_of_amino_acid_letter(const char &);
	std::string get_code_str_of_amino_acid_letter(const char &);

	char get_letter_of_amino_acid_code(const std::string &);

	std::ostream & operator<<(std::ostream &,
	                          const amino_acid &);

	std::istream & operator>>(std::istream &,
	                          amino_acid &);
} // namespace cath

#endif
