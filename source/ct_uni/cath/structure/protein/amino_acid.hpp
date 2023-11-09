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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_PROTEIN_AMINO_ACID_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_PROTEIN_AMINO_ACID_HPP

#include <iosfwd>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include <boost/operators.hpp>
#include <boost/throw_exception.hpp>

#include <fmt/core.h>

#include "cath/common/algorithm/contains.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/char_arr_type_aliases.hpp"
#include "cath/common/cpp20/arrays_are_equal.hpp"
#include "cath/common/cpp20/overload.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"
#include "cath/common/function/ident.hpp"
#include "cath/common/make_type_of_first_n.hpp"
#include "cath/common/optional/make_optional_if.hpp"
#include "cath/common/string/string_view_of_char_arr.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/file/pdb/pdb_record.hpp"
#include "cath/structure/protein/dna_atom.hpp"
#include "cath/structure/structure_type_aliases.hpp"

namespace cath {

	using ::std::literals::string_view_literals::operator""sv;

	/// \brief Represent the type of amino_acid record
	enum class amino_acid_type : char {
		AA,      ///< A standard amino_acid amino acid record (!)
		HETATOM, ///< A HETATM amino acid record
		DNA      ///< A DNA/RNA amino acid record
	};

	::std::string to_string(const amino_acid_type &);



	/// \brief TODOCUMENT
	///
	/// \todo Since, some part of this arose as taking excessive & non-trivial time in some profile
	///       it is worth seeing if it can be made faster, probably by using a ::std::array and constexpr_find()
	///       and constexpr_for_n() to generate lookups with minimal runtime overhead.
	constexpr auto LETTER_CODE_AND_NAME_LIST = ::std::array{
		::std::tuple( 'A', make_char_arr( "ALA" ), "Alanine"sv                            ),
		::std::tuple( 'B', make_char_arr( "ASX" ), "Ambiguous Asparagine/Aspartic Acid"sv ), // eg PDBs 156b, 1kp0, 1pgk, 2atc, 2fmd, 2rxn, 3atc, 3bcl and 3e2o
		::std::tuple( 'C', make_char_arr( "CYS" ), "Cysteine"sv                           ),
		::std::tuple( 'D', make_char_arr( "ASP" ), "Aspartic Acid"sv                      ),
		::std::tuple( 'E', make_char_arr( "GLU" ), "Glutamic Acid"sv                      ),
		::std::tuple( 'F', make_char_arr( "PHE" ), "Phenylalanine"sv                      ),
		::std::tuple( 'G', make_char_arr( "GLY" ), "Glycine"sv                            ),
		::std::tuple( 'H', make_char_arr( "HIS" ), "Histidine"sv                          ),
		::std::tuple( 'I', make_char_arr( "ILE" ), "Isoleucine"sv                         ),
		::std::tuple( 'J', make_char_arr( "XLE" ), "Leucine/Isoleucine"sv                 ), // eg PDBs?
		::std::tuple( 'K', make_char_arr( "LYS" ), "Lysine"sv                             ),
		::std::tuple( 'L', make_char_arr( "LEU" ), "Leucine"sv                            ),
		::std::tuple( 'M', make_char_arr( "MET" ), "Methionine"sv                         ),
		::std::tuple( 'N', make_char_arr( "ASN" ), "Asparagine"sv                         ),
		::std::tuple( 'O', make_char_arr( "PYL" ), "Pyrrolysine"sv                        ),
		::std::tuple( 'P', make_char_arr( "PRO" ), "Proline"sv                            ),
		::std::tuple( 'Q', make_char_arr( "GLN" ), "Glutamine"sv                          ),
		::std::tuple( 'R', make_char_arr( "ARG" ), "Arginine"sv                           ),
		::std::tuple( 'S', make_char_arr( "SER" ), "Serine"sv                             ),
		::std::tuple( 'T', make_char_arr( "THR" ), "Threonine"sv                          ),
		::std::tuple( 'U', make_char_arr( "SEC" ), "Selenocysteine"sv                     ), // eg PDBs 1aa6, 1cc1, 1fdi, 1fdo, 1h0h, 1kqf, 1kqg, 1pae, 1pfp, 2bc7, 2bc8, 2iv2, 2wpn, 2xsk, 3ean, 3eao, 3fwf, 3fwi, 3fwj, 3u5s, 3ze7, 3ze8, 3ze9, 3zea, 4kl8, 4kn9, 4ko1, 4ko2, 4ko3, 4ko4
		::std::tuple( 'V', make_char_arr( "VAL" ), "Valine"sv                             ),
		::std::tuple( 'W', make_char_arr( "TRP" ), "Tryptophan"sv                         ),
		::std::tuple( 'X', make_char_arr( "UNK" ), "Unknown"sv                            ),
		::std::tuple( 'Y', make_char_arr( "TYR" ), "Tyrosine"sv                           ),
		::std::tuple( 'Z', make_char_arr( "GLX" ), "Ambiguous Glutamine/Glutamic Acid"sv  )  // eg PDBs 156b, 1kp0, 1pgk, 2rxn, 2tnc, 3bcl and 4cpa
	};

	/// \brief TODOCUMENT
	///
	/// \param prm_letter The 1 letter code
	constexpr ::std::optional<uint> index_opt_of_letter( const char &prm_letter ) {
		for ( unsigned int x = 0; x < LETTER_CODE_AND_NAME_LIST.size(); ++x ) {
			if ( ::std::get<0>( LETTER_CODE_AND_NAME_LIST.at( x ) ) == prm_letter ) {
				return { x };
			}
		}
		return ::std::nullopt;
	}

	/// \brief TODOCUMENT
	///
	/// \param prm_letter The 1 letter code
	constexpr uint index_of_letter( const char &prm_letter ) {
		const auto index_opt = index_opt_of_letter( prm_letter );
		if ( index_opt ) {
			return *index_opt;
		}
		BOOST_THROW_EXCEPTION(
		  common::invalid_argument_exception( ::fmt::format( "Amino acid letter '{}' is not valid", prm_letter ) ) );
	}

	constexpr ::std::optional<uint> index_opt_of_code( const char_3_arr &prm_code ) {
		for ( unsigned int x = 0; x < LETTER_CODE_AND_NAME_LIST.size(); ++x ) {
			if ( common::arrays_are_equal( ::std::get<1>( LETTER_CODE_AND_NAME_LIST.at( x ) ), prm_code ) ) {
				return { x };
			}
		}
		return ::std::nullopt;
	}

	constexpr uint index_of_code( const char_3_arr &prm_code ) {
		const auto index_opt = index_opt_of_code( prm_code );
		if ( index_opt ) {
			return *index_opt;
		}
		BOOST_THROW_EXCEPTION( common::invalid_argument_exception(
		  ::fmt::format( "Amino acid code '{}' is not valid", common::string_view_of_char_arr( prm_code ) ) ) );
	}

	constexpr ::std::optional<uint> index_opt_of_name( const ::std::string_view &prm_name ) {
		for ( unsigned int x = 0; x < LETTER_CODE_AND_NAME_LIST.size(); ++x ) {
			if ( ::std::get<2>( LETTER_CODE_AND_NAME_LIST.at( x ) ) == prm_name ) {
				return { x };
			}
		}
		return ::std::nullopt;
	}

	constexpr uint index_of_name( const ::std::string_view &prm_name ) {
		const auto index_opt = index_opt_of_name( prm_name );
		if ( index_opt ) {
			return *index_opt;
		}
		BOOST_THROW_EXCEPTION(
		  common::invalid_argument_exception( ::fmt::format( "Amino acid name '{}' is not valid", prm_name ) ) );
	}

	/// \brief TODOCUMENT
	///
	/// \param prm_index TODOCUMENT
	template <typename T>
	constexpr T get_letter_code_or_name( const size_t &prm_index ) {
		return ::std::get<T>( LETTER_CODE_AND_NAME_LIST.at( prm_index ) );
	}

	namespace detail {

		/// \brief The number of chars used to store a HETATM
		constexpr size_t NUM_HETATM_CHARS_IN_AA_REPN = 3;

		/// TODOCUMENT
		///
		/// The typical size of aa_aa_repn is 4 bytes
		///
		/// This is an index into LETTER_CODE_AND_NAME_LIST
		using aa_aa_repn = uint;

		/// TODOCUMENT
		///
		/// The typical size of dna_aa_repn is 1 bytes
		using dna_aa_repn = dna_atom;

		/// TODOCUMENT
		///
		/// The typical size of hetatm_aa_repn is 3 bytes
		using hetatm_aa_repn = ::std::array<char, NUM_HETATM_CHARS_IN_AA_REPN>;

		/// TODOCUMENT
		///
		/// The typical size of aa_data is 8 bytes
		///
		/// If there's strong reason, this could be squeezed into 4 bytes, maybe even 3 at a push
		using aa_data = ::std::variant<aa_aa_repn, hetatm_aa_repn, dna_aa_repn>;


		/// \brief TODOCUMENT
		///
		/// This only matches true (ATOM) amino acids or DNA, not HETATMs (except for "ACE" and "NH2")
		/// so it only returns returns aa_aa_repn or dna_aa_repn, not hetatm_aa_repn (except for "ACE" and "NH2")
		///
		/// \param prm_id The three letter code
		constexpr aa_data atom_aa_or_dna_of_id( const char_3_arr &prm_id ) {
			// If this matches an AA code, return the index for that
			const ::std::optional<unsigned int> index_opt = index_opt_of_code( prm_id );
			if ( index_opt ) {
				return *index_opt;
			}

			// If this matches a DNA/RNA base, return the dna_atom for that
			const ::std::optional<dna_atom> dna_opt = dna_atom_of_code( prm_id );
			if ( dna_opt ) {
				return *dna_opt;
			}

			// Hacks to handle erroneous ATOM AA in versions of 2yjd from before it was fixed in January 2017
			//
			// From http://deposit.rcsb.org/format-faq-v1.html :
			//
			// > Noteworthy exceptions to the above treatment of modified residues are the cases of
			// > acetylation of the N-terminus (residue ACE) and amidation of the C-terminus (residue NH2).
			if ( common::arrays_are_equal( prm_id, make_char_arr( "ACE" ) ) ) {
				return prm_id;
			}
			if ( common::arrays_are_equal( prm_id, make_char_arr( "NH2" ) ) ) {
				return prm_id;
			}
			BOOST_THROW_EXCEPTION( common::invalid_argument_exception( ::fmt::format(
			  R"(Amino acid string "{}" has three characters but is not a recognised code (currently case-sensitive))",
			  common::string_view_of_char_arr( prm_id ) ) ) );
		}

		/// \brief TODOCUMENT
		///
		/// This only matches true (ATOM) amino acids or DNA, not HETATMs (except for "ACE" and "NH2")
		/// so it only returns returns aa_aa_repn or dna_aa_repn, not hetatm_aa_repn (except for "ACE" and "NH2")
		///
		/// \param prm_id The 1 letter code, three letter code or name to which the amino_acid should be set
		constexpr aa_data atom_aa_or_dna_of_id( const ::std::string_view &prm_id ) {
			// If the argument has more than three characters, try to match it against a name
			if ( prm_id.length() > 3 ) {
				return index_of_name( prm_id );
			}

			// If the argument has one character, try to match it against one of the letters
			if ( prm_id.length() == 1 ) {
				return index_of_letter( prm_id.front() );
			}

			// If the argument doesn't now have three characters, it must have two and that's invalid so throw
			if ( prm_id.length() != 3 ) {
				BOOST_THROW_EXCEPTION( common::invalid_argument_exception(
				  ::fmt::format( R"(Two-character amino acid id "{}" not recognised)", prm_id ) ) );
			}

			// Must be three characters
			return atom_aa_or_dna_of_id( common::make_array_of_first_n<3>( prm_id ) );
		}

		/// \brief TODOCUMENT
		///
		/// If pdb_record::ATOM, this acts as if that argument hadn't been passed (try AA code, try DNA code, throw)
		/// If pdb_record::HETATM, this tries for an AA code then falls back on HETATM
		///
		/// \param prm_id         TODOCUMENT
		/// \param prm_pdb_record TODOCUMENT
		constexpr aa_data aa_data_of_id_and_pdb_record( const char_3_arr &prm_id, const file::pdb_record &prm_pdb_record ) {
			switch ( prm_pdb_record ) {
				case ( file::pdb_record::ATOM ): {
					return detail::atom_aa_or_dna_of_id( prm_id );
				}
				case ( file::pdb_record::HETATM ): {
					const ::std::optional<uint> index_opt = index_opt_of_code( prm_id );
					return index_opt ? aa_data{ *index_opt } : aa_data{ prm_id };
				}
			}
			BOOST_THROW_EXCEPTION( common::invalid_argument_exception(
			  "Value of prm_pdb_record not recognised whilst constructing an amino_acid" ) );
		}

		/// \brief TODOCUMENT
		///
		/// If pdb_record::ATOM, this acts as if that argument hadn't been passed (try AA code, try DNA code, throw)
		/// If pdb_record::HETATM, this throws if not a 3-char code, tries for an AA code then falls back on HETATM
		///
		/// \param prm_id         TODOCUMENT
		/// \param prm_pdb_record TODOCUMENT
		constexpr aa_data aa_data_of_id_and_pdb_record( const ::std::string_view &prm_id, const file::pdb_record &prm_pdb_record ) {
			switch ( prm_pdb_record ) {
				case ( file::pdb_record::ATOM ): {
					return detail::atom_aa_or_dna_of_id( prm_id );
				}
				case ( file::pdb_record::HETATM ): {
					if ( prm_id.length() != 3 ) {
						BOOST_THROW_EXCEPTION( common::invalid_argument_exception(
						  "Cannot create a HETATM amino acid from a string that is not 3 characters long" ) );
					}
					const ::std::array str_as_arr = common::make_array_of_first_n<3>( prm_id );
					const ::std::optional<uint> index_opt = index_opt_of_code( str_as_arr );
					return index_opt ? aa_data{ *index_opt } : aa_data{ str_as_arr };
				}
			}
			BOOST_THROW_EXCEPTION( common::invalid_argument_exception(
			  "Value of prm_pdb_record not recognised whilst constructing an amino_acid" ) );
		}

		/// \brief Get the amino_acid_type of this amino_acid
		constexpr amino_acid_type get_type( const aa_data &prm_aa_data ) {
			// clang-format off
			return ::std::visit(
				common::overload(
					[]( aa_aa_repn     ) { return amino_acid_type::AA      ; },
					[]( dna_aa_repn    ) { return amino_acid_type::DNA     ; },
					[]( hetatm_aa_repn ) { return amino_acid_type::HETATOM ; }
				),
				prm_aa_data
			);
			// clang-format on
		}

		/// \brief TODOCUMENT
		constexpr char_3_arr get_code( const aa_data &prm_aa_data ) {
			// clang-format off
			return ::std::visit(
				common::overload(
					[] ( const aa_aa_repn     &x ) { return get_letter_code_or_name<char_3_arr>( x ) ; },
					[] ( const dna_aa_repn    &x ) { return to_three_char_arr( x )                   ; },
					[] ( const hetatm_aa_repn &x ) { return x                                        ; }
				),
				prm_aa_data
			);
			// clang-format on
		}

		/// \brief Get a single letter for the amino_acid, which is the expected letter for
		///        the standard amino acids and a few HETATMs (eg MSE -> M) and is 'X' otherwise
		/// \brief A visitor to get the one-letter code from the variant in amino_acid.
		///        This just returns 'X' for DNA or unrecognised HETATM.
		constexpr char get_letter_tolerantly( const aa_data &prm_aa_data ) {
			// clang-format off
			return ::std::visit(
				common::overload(
					[] ( const aa_aa_repn     &x ) { return get_letter_code_or_name<char>( x ) ; },
					[] ( const dna_aa_repn    &  ) { return 'X'                                ; },
					[] ( const hetatm_aa_repn &x ) {
						// Some of the most common codes, as extracted from PDB dir with command like:
						//
						//     ls -1 | xargs grep -hPB99999 '^TER   ' | grep -P '^HETATM' | awk '{print substr( $0, 18, 3 )}' | sort | uniq -c | sort -g
						//
						// Codes from pages like : http://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/MSE
						//
						// Default to 'X'
						return
							( x == make_char_arr( " IC" ) ) ? 'C' :
							( x == make_char_arr( " IG" ) ) ? 'G' :
							( x == make_char_arr( "HYP" ) ) ? 'P' :
							( x == make_char_arr( "LCG" ) ) ? 'G' :
							( x == make_char_arr( "MLY" ) ) ? 'K' :
							( x == make_char_arr( "MSE" ) ) ? 'M' :
							( x == make_char_arr( "NCX" ) ) ? 'N' :
							( x == make_char_arr( "OMG" ) ) ? 'G' :
							( x == make_char_arr( "PCA" ) ) ? 'E' :
							( x == make_char_arr( "PSU" ) ) ? 'U' :
							( x == make_char_arr( "PTR" ) ) ? 'Y' :
							( x == make_char_arr( "SEP" ) ) ? 'S' :
							                                  'X' ;
					}
				),
				prm_aa_data
			);
			// clang-format on
		}

	} // namespace detail

	/// \brief TODOCUMENT
	///
	/// \todo Make get_code() return char_3_arr and make LETTER_CODE_AND_NAME_LIST store char_3_arr
	class amino_acid final : private boost::equivalent      < amino_acid,
	                                 boost::totally_ordered < amino_acid > > {
	private:
		using char_size_unordered_map   = ::std::unordered_map<char,        uint>;
		using string_size_unordered_map = ::std::unordered_map<::std::string, uint>;

		/// \brief The data for the amino acid, one of three different types in the variant
		detail::aa_data data;

		[[nodiscard]] constexpr const detail::aa_aa_repn &check_is_proper_amino_acid() const;

	  public:
		constexpr amino_acid( const char_3_arr &, const file::pdb_record & );
		constexpr amino_acid( const ::std::string_view &, const file::pdb_record & );
		constexpr explicit amino_acid( const char & );
		constexpr explicit amino_acid( const char_3_arr & );
		constexpr explicit amino_acid( const ::std::string_view & );

		[[nodiscard]] constexpr amino_acid_type get_type() const;

		[[nodiscard]] constexpr const char_3_arr &get_hetatm_chars() const;

		[[nodiscard]] constexpr char_opt           get_letter_if_amino_acid() const;
		[[nodiscard]] constexpr char               get_letter_tolerantly() const;
		[[nodiscard]] constexpr char_3_arr         get_code() const;
		[[nodiscard]] constexpr ::std::string_view get_name() const;
	};

	amino_acid_vec make_amino_acids_of_chars(const char_vec &);

	/// \brief Check that this is a proper amino acid (rather than a HETATM record or DNA/RNA pseudo-amino-acid)
	///        and throw an exception if not, otherwise return the index value
	constexpr auto amino_acid::check_is_proper_amino_acid() const -> const detail::aa_aa_repn & {
		const detail::aa_aa_repn * const aa_ptr = ::std::get_if<detail::aa_aa_repn>( &data );
		if ( aa_ptr == nullptr ) {
			const auto code = get_code();
			BOOST_THROW_EXCEPTION(common::out_of_range_exception(
				R"(Cannot use a generic HETATM amino_acid or DNA value as a proper, ATOM-record amino acid. Problem was a )"
				+ to_string( get_type() )
				+ R"( with code ")"
				+ ::std::string( code.begin(), code.end() )
				+ R"(".)"
			));
		}
		return *aa_ptr;
	}

	/// \brief Ctor for amino_acid
	///
	/// If pdb_record::ATOM, this acts as if that argument hadn't been passed (try AA code, try DNA code, throw)
	/// If pdb_record::HETATM, this tries for an AA code then falls back on HETATM
	///
	/// \param pdb_id         TODOCUMENT
	/// \param prm_pdb_record TODOCUMENT
	constexpr amino_acid::amino_acid( const char_3_arr &pdb_id, const file::pdb_record &prm_pdb_record ) :
	        data{ detail::aa_data_of_id_and_pdb_record( pdb_id, prm_pdb_record ) } {
	}

	/// \brief Ctor for amino_acid
	///
	/// If pdb_record::ATOM, this acts as if that argument hadn't been passed (try AA code, try DNA code, throw)
	/// If pdb_record::HETATM, this throws if not a 3-char code, tries for an AA code then falls back on HETATM
	///
	/// \param prm_string     TODOCUMENT
	/// \param prm_pdb_record TODOCUMENT
	constexpr amino_acid::amino_acid( const ::std::string_view &prm_string, const file::pdb_record &prm_pdb_record ) :
	        data{ detail::aa_data_of_id_and_pdb_record( prm_string, prm_pdb_record ) } {
	}

	/// \brief Ctor for amino_acid
	///
	/// \param prm_letter The 1 letter code
	constexpr amino_acid::amino_acid(const char &prm_letter
	                                 ) : data{ index_of_letter( prm_letter ) } {
	}

	/// \brief Ctor for amino_acid
	constexpr amino_acid::amino_acid( const char_3_arr &prm_code ) : data{ detail::atom_aa_or_dna_of_id( prm_code ) } {
	}

	/// \brief Ctor for amino_acid
	constexpr amino_acid::amino_acid( const ::std::string_view &prm_id ) :
	        data{ detail::atom_aa_or_dna_of_id( prm_id ) } {
	}

	/// \brief Get the amino_acid_type of this amino_acid
	constexpr amino_acid_type amino_acid::get_type() const {
		return detail::get_type( data );
	}

	/// \brief Getter for the raw string for HETATM amino_acids
	constexpr const char_3_arr & amino_acid::get_hetatm_chars() const {
		const detail::hetatm_aa_repn * const raw_string_ptr = ::std::get_if<detail::hetatm_aa_repn>( &data );
		if ( raw_string_ptr == nullptr ) {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Unable to get raw_string from non-HETATM amino_acid"));
		}
		return *raw_string_ptr;
	}

	/// \brief Get the single letter for the amino_acid if it is an amino acid, or nullopt otherwise
	constexpr char_opt amino_acid::get_letter_if_amino_acid() const {
		const detail::aa_aa_repn * const aa_ptr = ::std::get_if<detail::aa_aa_repn>( &data );
		return if_then_optional(
			aa_ptr != nullptr,
			( get_letter_code_or_name<char>( *aa_ptr ) )
		);
	}

	/// \brief Get a single letter for the amino_acid, which is the expected letter for
	///        the standard amino acids and a few HETATMs (eg MSE -> M) and is 'X' otherwise
	constexpr char amino_acid::get_letter_tolerantly() const {
		return detail::get_letter_tolerantly( data );
	}

	constexpr char_3_arr amino_acid::get_code() const {
		return detail::get_code( data );
	}

	/// \brief TODOCUMENT
	constexpr ::std::string_view amino_acid::get_name() const {
		return get_letter_code_or_name<::std::string_view>( check_is_proper_amino_acid() );
	}

	/// \brief Get a string of the three-letter-code for the specified amino_acid
	///
	/// \relates amino_acid
	inline ::std::string get_code_string(const amino_acid &prm_amino_acid ///< The amino_acid to query
	                                     ) {
		return common::string_of_char_arr( prm_amino_acid.get_code() );
	}

	namespace detail {

		/// \brief Make a less-than comparator for the specified amino_acid
		///
		/// TODO: Redo this:
		///         * no need for the char field
		///         * change it to operate on an aa_data argument
		///         * add constexpr
		///         * use in a friend constexpr operator< in amino_acid
		inline ::std::tuple<uint8_t, char, char_3_arr> make_amino_acid_lt_comparator(const amino_acid &prm_aa ///< The amino_acid to query
		                                                                             ) {
			switch ( prm_aa.get_type() ) {
				case ( amino_acid_type::AA ) : {
					return { 0, *prm_aa.get_letter_if_amino_acid(), char_3_arr{ { 0, 0, 0 } } };
				}
				case ( amino_acid_type::HETATOM ) : {
					return { 1, 0,                                  prm_aa.get_hetatm_chars() };
				}
				case ( amino_acid_type::DNA ) : {
					return { 2, 0,                                  char_3_arr{ { 0, 0, 0 } } };
				}
			}
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Value of prm_aa.get_type() not recognised whilst in make_amino_acid_lt_comparator()"));
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
	inline bool operator<(const amino_acid &prm_amino_acid_1, ///< The first  amino acid to compare
	                      const amino_acid &prm_amino_acid_2  ///< The second amino acid to compare
	                      ) {
		return (
			detail::make_amino_acid_lt_comparator( prm_amino_acid_1 )
			<
			detail::make_amino_acid_lt_comparator( prm_amino_acid_2 )
		);
	}

	/// \brief Return whether the specified amino_acid is a proper amino_acid (rather than a HETATM record or DNA/RNA pseudo-amino-acid)
	///
	/// \relates amino_acid
	inline bool is_proper_amino_acid(const amino_acid &prm_amino_acid ///< The amino_acid to query
	                                 ) {
		return ( prm_amino_acid.get_type() == amino_acid_type::AA );
	}

	/// \brief Return whether the specified amino_acid is for water (HOH)
	///
	/// \relates amino_acid
	inline bool is_water(const amino_acid &prm_amino_acid ///< The amino_acid to query
	                     ) {
		return ( prm_amino_acid.get_code() == make_char_arr( "HOH" ) );
	}

	char_3_arr get_code_of_amino_acid_letter( const char & );

	::std::string get_code_str_of_amino_acid_letter( const char & );

	char get_letter_of_amino_acid_code( const ::std::string_view & );

	::std::ostream &operator<<( ::std::ostream &, const amino_acid & );

	::std::istream &operator>>( ::std::istream &, amino_acid & );

} // namespace cath

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_PROTEIN_AMINO_ACID_HPP
