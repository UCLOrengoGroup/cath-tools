/// \file
/// \brief The pdb_atom class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_PDB_ATOM_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_PDB_ATOM_HPP

#include <iosfwd>
#include <limits>
#include <string_view>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast/bad_lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>

#include "cath/biocore/chain_label.hpp"
#include "cath/biocore/residue_id.hpp"
#include "cath/common/algorithm/contains.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/string/string_parse_tools.hpp"
#include "cath/file/file_type_aliases.hpp"
#include "cath/file/pdb/element_type_string.hpp"
#include "cath/file/pdb/pdb_atom_parse_status.hpp"
#include "cath/file/pdb/pdb_record.hpp"
#include "cath/structure/geometry/coord.hpp"
#include "cath/structure/protein/amino_acid.hpp"

// clang-format off
namespace cath::geom { class rotation; }
// clang-format on

namespace cath::util {

	template <typename T>
	[[nodiscard]] constexpr bool constexpr_isfinite( const T &t ) {
		static_assert( ::std::is_same_v<T, double> || ::std::is_same_v<T, float> );
		// clang-format off
		return (
				( t == t )
				&&
				( t != ::std::numeric_limits<T>::infinity() )
				&&
				( t != ( 0.0 -::std::numeric_limits<T>::infinity() ) )
		);
		// clang-format on
	}

} // namespace cath::util

namespace cath::file {

	/// \brief TODOCUMENT
	class pdb_atom final {
	private:
		/// \brief TODOCUMENT
		geom::coord         atom_coord;

		/// \brief TODOCUMENT
		amino_acid          the_amino_acid;

		/// \brief TODOCUMENT
		element_type_string the_element_type;

		/// \brief TODOCUMENT
		pdb_record          record_type;

		/// \brief TODOCUMENT
		char                alt_locn;

		/// \brief TODOCUMENT
		uint                atom_serial;

		/// \brief TODOCUMENT
		float               occupancy;

		/// \brief TODOCUMENT
		float               temp_factor;

		/// \brief TODOCUMENT
		char_2_arr          element_symbol;

		/// \brief TODOCUMENT
		char_2_arr          charge;

	public:
		constexpr pdb_atom(const pdb_record &,
		                   const uint &,
		                   const char_4_arr &,
		                   const char &,
		                   amino_acid,
		                   geom::coord,
		                   const float &,
		                   const float &,
		                   char_2_arr,
		                   char_2_arr);

		[[nodiscard]] constexpr const pdb_record & get_record_type() const;
		[[nodiscard]] constexpr const uint &       get_atom_serial() const;
		[[nodiscard]] constexpr const char_4_arr & get_element_type_untrimmed() const;
		[[nodiscard]] constexpr ::std::string_view get_element_type() const;
		[[nodiscard]] constexpr const char &       get_alt_locn() const;
		[[nodiscard]] constexpr const amino_acid & get_amino_acid() const;
		[[nodiscard]] constexpr const geom::coord &get_coord() const;
		[[nodiscard]] constexpr const float &      get_occupancy() const;
		[[nodiscard]] constexpr const float &      get_temp_factor() const;
		[[nodiscard]] constexpr const char_2_arr & get_element_symbol() const;
		[[nodiscard]] constexpr const char_2_arr & get_charge() const;

		void rotate(const geom::rotation &);
		void operator+=(const geom::coord &);
		void operator-=(const geom::coord &);

		/// \brief Specify how short a line can be before it will be rejected.
		static constexpr size_t MIN_NUM_PDB_COLS =  66;

		/// \brief Specify how long a line can be before it will be rejected.
		///
		/// This is set to 160, twice the standard 80, to avoid rejecting lines that have just got some extra stuff at the end
		static constexpr size_t MAX_NUM_PDB_COLS =  2 * 80;
	};

	constexpr ::std::string_view get_element_type_untrimmed_str_ref(const pdb_atom &);
	constexpr coarse_element_type get_coarse_element_type(const pdb_atom &);
	char get_amino_acid_letter_tolerantly(const pdb_atom &);
	char_3_arr get_amino_acid_code(const pdb_atom &);
	std::string get_amino_acid_code_string(const pdb_atom &);
	std::string_view get_amino_acid_name( const pdb_atom & );
	bool is_pdb_record_of_type(const std::string &,
	                           const pdb_record &);

	std::ostream & write_pdb_file_entry(std::ostream &,
	                                    const residue_id &,
	                                    const pdb_atom &);
	std::string to_pdb_file_entry(const residue_id &,
	                              const pdb_atom &);
	bool alt_locn_is_dssp_accepted(const pdb_atom &);
	std::ostream & operator<<(std::ostream &,
	                          const pdb_atom &);
	::std::string_view get_element_symbol_str_ref(const pdb_atom &);
	::std::string_view get_charge_str_ref(const pdb_atom &);

	/// \brief Ctor for pdb_atom
	constexpr pdb_atom::pdb_atom(const pdb_record &prm_record_type,    ///< TODOCUMENT
	                             const uint       &prm_atom_serial,    ///< TODOCUMENT
	                             const char_4_arr &prm_element_type,   ///< TODOCUMENT
	                             const char       &prm_alt_locn,       ///< TODOCUMENT
	                             amino_acid        prm_amino_acid,     ///< TODOCUMENT
	                             geom::coord       prm_coord,          ///< TODOCUMENT
	                             const float      &prm_occupancy,      ///< TODOCUMENT
	                             const float      &prm_temp_factor,    ///< TODOCUMENT
	                             char_2_arr        prm_element_symbol, ///< TODOCUMENT
	                             char_2_arr        prm_charge          ///< TODOCUMENT
	                             ) : atom_coord       { std::move( prm_coord )          },
	                                 the_amino_acid   { std::move( prm_amino_acid )     },
	                                 the_element_type { prm_element_type                },
	                                 record_type      { prm_record_type                 },
	                                 alt_locn         { prm_alt_locn                    },
	                                 atom_serial      { prm_atom_serial                 },
	                                 occupancy        { prm_occupancy                   },
	                                 temp_factor      { prm_temp_factor                 },
	                                 element_symbol   ( std::move( prm_element_symbol ) ), //< Don't change these brackets to braces - it breaks the build on the older Clang on Travis-CI
	                                 charge           ( std::move( prm_charge )         )  //< Don't change these brackets to braces - it breaks the build on the older Clang on Travis-CI
	                                 {
		if ( ! util::constexpr_isfinite( occupancy ) ) {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Argument occupancy must be a normal, finite floating-point number"));
		}
		if ( ! util::constexpr_isfinite( temp_factor ) ) {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Argument temp_factor must be a normal, finite floating-point number"));
		}
	}

	/// \brief TODOCUMENT
	constexpr const pdb_record & pdb_atom::get_record_type() const {
		return record_type;
	}

	/// \brief TODOCUMENT
	constexpr const uint & pdb_atom::get_atom_serial() const {
		return atom_serial;
	}

	/// \brief TODOCUMENT
	constexpr const char_4_arr & pdb_atom::get_element_type_untrimmed() const {
		return the_element_type.get_element_type_untrimmed();
	}

	/// \brief TODOCUMENT
	constexpr ::std::string_view pdb_atom::get_element_type() const {
		return the_element_type.get_element_type();
	}

	/// \brief TODOCUMENT
	constexpr const char & pdb_atom::get_alt_locn() const {
		return alt_locn;
	}

	/// \brief TODOCUMENT
	constexpr const amino_acid & pdb_atom::get_amino_acid() const {
		return the_amino_acid;
	}

	/// \brief TODOCUMENT
	constexpr const geom::coord & pdb_atom::get_coord() const {
		return atom_coord;
	}

	/// \brief TODOCUMENT
	constexpr const float & pdb_atom::get_occupancy() const {
		return occupancy;
	}

	/// \brief TODOCUMENT
	constexpr const float & pdb_atom::get_temp_factor() const {
		return temp_factor;
	}

	/// \brief Getter for the element symbol characters
	constexpr const char_2_arr & pdb_atom::get_element_symbol() const {
		return element_symbol;
	}

	/// \brief Getter for the charge characters
	constexpr const char_2_arr & pdb_atom::get_charge() const {
		return charge;
	}

	/// \brief Get a string_view for the untrimmed element type of the specified pdb_atom
	constexpr ::std::string_view get_element_type_untrimmed_str_ref(const pdb_atom &prm_pdb_atom ///< The pdb_atom to query
	                                                               ) {
		return common::string_view_of_char_arr( prm_pdb_atom.get_element_type_untrimmed() );
	}

	/// \brief Get the coarse_element_type corresponding to the specified pdb_atom
	constexpr coarse_element_type get_coarse_element_type(const pdb_atom &prm_pdb_atom ///< The pdb_atom to query
	                                                      ) {
		return get_coarse_element_type( prm_pdb_atom.get_element_type() );
	}

	/// \brief Return whether this line purports to be an record of the specified pdb_record
	///
	/// Note: This does NOT check whether the line is a valid record;
	///       Use atom_record_parse_problem() for that.
	inline bool is_pdb_record_of_type(const std::string &prm_pdb_record_string, ///< The string to check
	                                  const pdb_record  &prm_record_type        ///< TODOCUMENT
	                                  ) {
		if ( prm_pdb_record_string.empty() ) {
			return false;
		}
		switch ( prm_record_type ) {
			case ( pdb_record::ATOM   ) : { return boost::algorithm::starts_with( prm_pdb_record_string, "ATOM  " ); }
			case ( pdb_record::HETATM ) : { return boost::algorithm::starts_with( prm_pdb_record_string, "HETATM" ); }
		}
		BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Value of pdb_record not recognised whilst checking is_pdb_record_of_type()"));
	}

	/// \brief TODOCUMENT
	///
	/// \relates pdb_atom
	inline bool alt_locn_is_dssp_accepted(const pdb_atom &prm_pdb_atom ///< TODOCUMENT
	                                      ) {
		return prm_pdb_atom.get_alt_locn() == ' ' || prm_pdb_atom.get_alt_locn() == 'A';;
	}

	/// \brief Convert the specified three-letter string and pdb_record to an amino_acid
	///
	/// This is less strict than the amino_acid ctor because it allows certain combinations
	/// to be accepted for decay to UNK/X. See function body for the accepted combinations.
	///
	/// \relates pdb_atom
	inline amino_acid get_amino_acid_of_string_and_record(const char_3_arr &prm_aa_string,  ///< The three-letter amino acid string (eg "SER")
	                                                      const pdb_record &prm_record_type ///< Whether this is an ATOM or HETATM record
	                                                      ) {
		try {
			return { prm_aa_string, prm_record_type };
		}
		catch (...) {
			constexpr ::std::array PDB_AA_ENTRIES_ALLOWED_TO_DECAY_TO_UNK = {
				::std::pair{ make_char_arr( "MLY" ), pdb_record::HETATM },
				::std::pair{ make_char_arr( "MSE" ), pdb_record::ATOM   },
				::std::pair{ make_char_arr( "MSE" ), pdb_record::HETATM },
			};
			if ( common::contains( PDB_AA_ENTRIES_ALLOWED_TO_DECAY_TO_UNK, ::std::pair( prm_aa_string, prm_record_type ) ) ) {
				return { "UNK", prm_record_type };
			}
			throw;
		}
	}

	/// \brief Parse the amino_acid from the specified PDB atom record string using the
	///        looser pdb_atom criteria that allow a few extra amino acids, which decay to UNK/X.
	///
	/// See the body of get_amino_acid_of_string_and_record() for the accepted combinations.
	///
	/// \relates pdb_atom
	inline amino_acid parse_amino_acid_from_pdb_atom_record(const ::std::string_view &prm_pdb_atom_record_string ///< TODOCUMENT
	                                                        ) {
		assert( prm_pdb_atom_record_string.length() >= 20 );
		const pdb_record         record_type = pdb_rec_of_six_chars_in_string( prm_pdb_atom_record_string                ); //  1 -  6        Record name   "ATOM  "
		const ::std::string_view the_a_a     =                                 prm_pdb_atom_record_string.substr( 17, 3  ); // 18 - 20        Residue name  resName      Residue name.
		return get_amino_acid_of_string_and_record(
			char_3_arr{ the_a_a[0], the_a_a[1], the_a_a[2] },
			record_type
		);
	}

	/// \brief Return a string containing the parse problem with a PDB ATOM/HETATM record string or "" if no problem
	///
	/// \relates pdb_atom
	///
	/// Note: This does NOT check whether the line is a valid record.
	///       Use atom_record_parse_problem() for that.
	inline status_string_aa_tuple pdb_record_parse_problem(const std::string &prm_pdb_atom_record_string ///< The string to check
	                                                       ) {
		amino_acid the_aa{ 'X' };
		if ( ! is_pdb_record_of_type( prm_pdb_atom_record_string, pdb_record::ATOM ) && ! is_pdb_record_of_type( prm_pdb_atom_record_string, pdb_record::HETATM ) ) {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot check for atom record parse problems because string is not an ATOM record"));
		}
		if ( prm_pdb_atom_record_string.length() < pdb_atom::MIN_NUM_PDB_COLS ) {
			return { pdb_atom_parse_status::ABORT, "Is too long", the_aa };
		}
		if ( prm_pdb_atom_record_string.length() > pdb_atom::MAX_NUM_PDB_COLS ) {
			return { pdb_atom_parse_status::ABORT, "Is too long", the_aa };
		}
		if ( prm_pdb_atom_record_string.at( 11 ) != ' '   ) {
			return { pdb_atom_parse_status::ABORT, "Does not contain a space at column 12", the_aa };
		}
		if ( prm_pdb_atom_record_string.at( 20 ) != ' '   ) {
			return { pdb_atom_parse_status::ABORT, "Does not contain a space at column 21", the_aa };
		}
		if ( prm_pdb_atom_record_string.at( 27 ) != ' ' || prm_pdb_atom_record_string.at( 28 ) != ' ' || prm_pdb_atom_record_string.at( 29 ) != ' ' ) {
			return { pdb_atom_parse_status::ABORT, "Does not contain spaces at columns 28-30", the_aa };
		}
	//	if ( prm_pdb_atom_record_string.at( 16 ) != ' ' && prm_pdb_atom_record_string.at( 16 ) != 'A' ) {
	//		return { pdb_atom_parse_status::SKIP, "Has alternate location indicator other than 'A' or' '", the_aa };
	//	}
		try {
			the_aa = parse_amino_acid_from_pdb_atom_record( prm_pdb_atom_record_string );
		}
		catch (const boost::exception &e) {
			return {
				pdb_atom_parse_status::SKIP,
				"Do not recognise amino acid entry: \"" + prm_pdb_atom_record_string.substr( 17, 3 ) + "\" - " + diagnostic_information( e ),
				the_aa
			};
		}
		catch (...) {
			return {
				pdb_atom_parse_status::SKIP,
				"Do not recognise amino acid entry: " + prm_pdb_atom_record_string.substr( 17, 3 ),
				the_aa
			};
		}
		return { pdb_atom_parse_status::OK, "", the_aa };
	}

	/// \brief TODOCUMENT
	///
	/// \relates pdb_atom
	inline resid_atom_pair parse_pdb_atom_record(const std::string &prm_pdb_atom_record_string, ///< TODOCUMENT
	                                             const amino_acid  &prm_amino_acid
	                                             ) {
		if ( isspace( prm_pdb_atom_record_string.at( 25 ) ) != 0 ) {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("PDB ATOM/HETATOM record malformed: space found in column 26 which should contain the end of a right justified residue number"));
		}

		try {
			                                                                                                                       // Comments with PDB format documentation
			                                                                                                                       // (http://www.wwpdb.org/documentation/format33/sect9.html#ATOM)
			const pdb_record       record_type =       pdb_rec_of_six_chars_in_string( prm_pdb_atom_record_string                ); //  1 -  6        Record name   "ATOM  "
			const uint             serial      =    common::parse_uint_from_substring( prm_pdb_atom_record_string,         6, 5  ); //  7 - 11        Integer       serial       Atom  serial number.
			const char_4_arr       element     = common::get_char_arr_of_substring<4>( prm_pdb_atom_record_string,        12     ); // 13 - 16        Atom          name         Atom name.
			const char            &alt_locn    =                                       prm_pdb_atom_record_string.at(     16     ); // 17             Character     altLoc       Alternate location indicator.
	//		const std::string      the_a_a     =                                       prm_pdb_atom_record_string.substr( 17, 3  ); // 18 - 20        Residue name  resName      Residue name.
			const char            &chain_char  =                                       prm_pdb_atom_record_string.at(     21     ); // 22             Character     chainID      Chain identifier.
			const int              res_num     =     common::parse_int_from_substring( prm_pdb_atom_record_string,        22, 4  ); // 23 - 26        Integer       resSeq       Residue sequence number.
			const char            &insert_code =                                       prm_pdb_atom_record_string.at    ( 26     ); // 27             AChar         iCode        Code for insertion of residues.
			const double           coord_x     =  common::parse_double_from_substring( prm_pdb_atom_record_string,        30, 8  ); // 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
			const double           coord_y     =  common::parse_double_from_substring( prm_pdb_atom_record_string,        38, 8  ); // 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
			const double           coord_z     =  common::parse_double_from_substring( prm_pdb_atom_record_string,        46, 8  ); // 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
			const float            occupancy   =  common::parse_float_from_substring ( prm_pdb_atom_record_string,        54, 6  ); // 55 - 60        Real(6.2)     occupancy    Occupancy.
			const float            temp_factor =  common::parse_float_from_substring ( prm_pdb_atom_record_string,        60, 6  ); // 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
			const char_2_arr       element_sym = common::get_char_arr_of_substring<2>( prm_pdb_atom_record_string,        76     ); // 77 - 78        LString(2)    element      Element symbol, right-justified.
			const char_2_arr       charge      = common::get_char_arr_of_substring<2>( prm_pdb_atom_record_string,        78     ); // 79 - 80        LString(2)    charge       Charge  on the atom.

			return {
				residue_id{
					chain_label( chain_char ),
					make_residue_name_with_non_insert_char( res_num, insert_code, ' ' )
				},
				pdb_atom(
					record_type,
					serial,
					element,
					alt_locn,
					prm_amino_acid,
					geom::coord{ coord_x, coord_y, coord_z },
					occupancy,
					temp_factor,
					element_sym,
					charge
				)
			};
		}
		catch (const boost::bad_lexical_cast &) {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Unable to cast a column whilst parsing a PDB ATOM record, which probably means it's malformed.\nRecord was \"" + prm_pdb_atom_record_string + "\""));
		}
		catch (const std::invalid_argument &) {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Unable to cast a column whilst parsing a PDB ATOM record, which probably means it's malformed.\nRecord was \"" + prm_pdb_atom_record_string + "\""));
		}
		catch (const std::out_of_range &) {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Casted column out of range whilst parsing a PDB ATOM record, which probably means it's malformed.\nRecord was \"" + prm_pdb_atom_record_string + "\""));
		}
	}

	/// \brief TODOCUMENT
	///
	/// \relates pdb_atom
	inline resid_atom_pair parse_pdb_atom_record(const std::string &prm_pdb_atom_record_string ///< TODOCUMENT
	                                             ) {
		return parse_pdb_atom_record(
			prm_pdb_atom_record_string,
			parse_amino_acid_from_pdb_atom_record( prm_pdb_atom_record_string )
		);
	}

	/// \brief Get a string_view to (the non-null part of) the specified pdb_atom's element_symbol string
	///
	/// \relates pdb_atom
	inline ::std::string_view get_element_symbol_str_ref(const pdb_atom &prm_pdb_atom ///< The pdb_atom to query
	                                                    ) {
		return common::string_view_of_null_term_char_arr( prm_pdb_atom.get_element_symbol() );
	}

	/// \brief Get a string_view to (the non-null part of) the specified pdb_atom's charge string
	///
	/// \relates pdb_atom
	inline ::std::string_view get_charge_str_ref(const pdb_atom &prm_pdb_atom ///< The pdb_atom to query
	                                            ) {
		return common::string_view_of_null_term_char_arr( prm_pdb_atom.get_charge() );
	}

	/// \brief Return whether the specified pdb_atom's amino_acid is for water (HOH)
	///
	/// \relates pdb_atom
	inline bool is_water(const pdb_atom &prm_pdb_atom ///< The pdb_atom to query
	                     ) {
		return is_water( prm_pdb_atom.get_amino_acid() );
	}

} // namespace cath::file

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_PDB_ATOM_HPP
