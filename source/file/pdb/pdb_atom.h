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

#ifndef _CATH_TOOLS_SOURCE_FILE_PDB_PDB_ATOM_H
#define _CATH_TOOLS_SOURCE_FILE_PDB_PDB_ATOM_H

#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast/bad_lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/utility/string_ref.hpp>

#include "common/algorithm/contains.h"
#include "common/string/string_parse_tools.h"
#include "exception/invalid_argument_exception.h"
#include "file/file_type_aliases.h"
#include "file/pdb/element_type_string.h"
#include "file/pdb/pdb_atom_parse_status.h"
#include "file/pdb/pdb_base.h"
#include "file/pdb/pdb_record.h"
#include "structure/chain_label.h"
#include "structure/geometry/coord.h"
#include "structure/protein/amino_acid.h"
#include "structure/residue_name.h"

#include <iosfwd>
namespace cath { namespace geom { class rotation; } }

namespace cath {
	namespace file {

		/// \brief TODOCUMENT
		class pdb_atom final {
		private:
			/// \brief TODOCUMENT
			pdb_record  record_type;

			/// \brief TODOCUMENT
			size_t      atom_serial;

			/// \brief TODOCUMENT
			element_type_string the_element_type;

			/// \brief TODOCUMENT
			char        alt_locn;

			/// \brief TODOCUMENT
			amino_acid  the_amino_acid;

			/// \brief TODOCUMENT
			geom::coord atom_coord;

			/// \brief TODOCUMENT
			double      occupancy;

			/// \brief TODOCUMENT
			double      temp_factor;

		public:
			pdb_atom(const pdb_record &,
			         const size_t &,
			         const std::string &&,
			         const char &,
			         const amino_acid &,
			         const geom::coord &,
			         const double &,
			         const double &);

			const pdb_record & get_record_type() const;
			const size_t & get_atom_serial() const;
			const std::string & get_element_type_untrimmed() const;
			const boost::string_ref & get_element_type() const;
			const char & get_alt_locn() const;
			const amino_acid & get_amino_acid() const;
			const geom::coord & get_coord() const;
			const double & get_occupancy() const;
			const double & get_temp_factor() const;

			void rotate(const geom::rotation &);
			void operator+=(const geom::coord &);
			void operator-=(const geom::coord &);

			static const std::string PDB_ID_NITROGEN;
			static const std::string PDB_ID_CARBON_ALPHA;
			static const std::string PDB_ID_CARBON;

			static const std::string PDB_ID_CARBON_BETA;
			static const std::string PDB_ID_OXYGEN;
		};

		char get_amino_acid_letter(const pdb_atom &);
		std::string get_amino_acid_code(const pdb_atom &);
		std::string get_amino_acid_name(const pdb_atom &);
		bool is_pdb_record_of_type(const std::string &,
		                           const pdb_record &);

		std::ostream & write_pdb_file_entry(std::ostream &,
		                                    const chain_label &,
		                                    const residue_name &,
		                                    const pdb_atom &);
		std::ostream & operator<<(std::ostream &,
		                          const pdb_atom &);

		/// \brief Ctor for pdb_atom
		inline pdb_atom::pdb_atom(const pdb_record   &arg_record_type,  ///< TODOCUMENT
		                          const size_t       &arg_atom_serial,  ///< TODOCUMENT
		                          const std::string &&arg_element_type, ///< TODOCUMENT
		                          const char         &arg_alt_locn,     ///< TODOCUMENT
		                          const amino_acid   &arg_amino_acid,   ///< TODOCUMENT
		                          const geom::coord  &arg_coord,        ///< TODOCUMENT
		                          const double       &arg_occupancy,    ///< TODOCUMENT
		                          const double       &arg_temp_factor   ///< TODOCUMENT
		                          ) : record_type     ( arg_record_type               ),
		                              atom_serial     ( arg_atom_serial               ),
		                              the_element_type( std::move( arg_element_type ) ),
		                              alt_locn        ( arg_alt_locn                  ),
		                              the_amino_acid  ( arg_amino_acid                ),
		                              atom_coord      ( arg_coord                     ),
		                              occupancy       ( arg_occupancy                 ),
		                              temp_factor     ( arg_temp_factor               ) {
			if ( ! boost::math::isfinite( occupancy ) ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Argument occupancy must be a normal, finite floating-point number"));
			}
			if ( ! boost::math::isfinite( temp_factor ) ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Argument temp_factor must be a normal, finite floating-point number"));
			}
		}


		/// \brief Return whether this line purports to be an record of the specified pdb_record
		///
		/// Note: This does NOT check whether the line is a valid record;
		///       Use atom_record_parse_problem() for that.
		inline bool is_pdb_record_of_type(const std::string &arg_pdb_record_string, ///< The string to check
		                                  const pdb_record  &arg_record_type        ///< TODOCUMENT
		                                  ) {
			if ( arg_pdb_record_string.empty() ) {
				return false;
			}
			switch ( arg_record_type ) {
				case ( pdb_record::ATOM   ) : { return boost::algorithm::starts_with( arg_pdb_record_string, "ATOM  " ); }
				case ( pdb_record::HETATM ) : { return boost::algorithm::starts_with( arg_pdb_record_string, "HETATM" ); }
				default : {
					BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Value of pdb_record not recognised whilst checking is_pdb_record_of_type()"));
					return false; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
				}
			}
		}

		/// \brief Convert the specified three-letter string and pdb_record to an amino_acid
		///
		/// This is less strict than the amino_acid ctor because it allows certain combinations
		/// to be accepted for decay to UNK/X. See function body for the accepted combinations.
		///
		/// \relates pdb_atom
		inline amino_acid get_amino_acid_of_string_and_record(const std::string &arg_aa_string,  ///< The three-letter amino acid string (eg "SER")
		                                                      const pdb_record  &arg_record_type ///< Whether this is an ATOM or HETATM record
		                                                      ) {
			if ( arg_aa_string.length() != 3 ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Amino acid string must contain three characters"));
			}
			try {
				return { arg_aa_string, arg_record_type };
			}
			catch (...) {
				const std::set<std::pair<std::string, pdb_record>> pdb_aa_entries_allowed_to_decay_to_unk = {
					{ "MLY", pdb_record::HETATM },
					{ "MSE", pdb_record::ATOM   },
					{ "MSE", pdb_record::HETATM },
				};
				if ( common::contains( pdb_aa_entries_allowed_to_decay_to_unk, std::make_pair( arg_aa_string, arg_record_type ) ) ) {
					return { "UNK", arg_record_type };
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
		inline amino_acid parse_amino_acid_from_pdb_atom_record(const std::string &arg_pdb_atom_record_string ///< TODOCUMENT
		                                                        ) {
			const pdb_record  record_type = pdb_rec_of_six_chars_in_string( arg_pdb_atom_record_string                ); //  1 -  6        Record name   "ATOM  "
			const std::string the_a_a     =                                 arg_pdb_atom_record_string.substr( 17, 3  ); // 18 - 20        Residue name  resName      Residue name.
			return get_amino_acid_of_string_and_record(
				the_a_a,
				record_type
			);
		}

		/// \brief Return a string containing the parse problem with a PDB ATOM/HETATM record string or "" if no problem
		///
		/// \relates pdb_atom
		///
		/// Note: This does NOT check whether the line is a valid record.
		///       Use atom_record_parse_problem() for that.
		///
		/// \todo Come C++17 (or whenever the explicit on std::tuple's ctor is made dependent on whether any of the tuple's types are explicit),
		///       remove the leading status_string_aa_tuple from all the return expressions.
		inline status_string_aa_tuple pdb_record_parse_problem(const std::string &arg_pdb_atom_record_string ///< The string to check
		                                                       ) {
			amino_acid the_aa{ 'X' };
			if ( ! is_pdb_record_of_type( arg_pdb_atom_record_string, pdb_record::ATOM ) && ! is_pdb_record_of_type( arg_pdb_atom_record_string, pdb_record::HETATM ) ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot check for atom record parse problems because string is not an ATOM record"));
			}
			if ( arg_pdb_atom_record_string.length() < pdb_base::MIN_NUM_PDB_COLS ) {
				return status_string_aa_tuple{ pdb_atom_parse_status::ABORT, "Is too long", the_aa };
			}
			if ( arg_pdb_atom_record_string.length() > pdb_base::MAX_NUM_PDB_COLS ) {
				return status_string_aa_tuple{ pdb_atom_parse_status::ABORT, "Is too long", the_aa };
			}
			if ( arg_pdb_atom_record_string.at( 11 ) != ' '   ) {
				return status_string_aa_tuple{ pdb_atom_parse_status::ABORT, "Does not contain a space at column 12", the_aa };
			}
			if ( arg_pdb_atom_record_string.at( 20 ) != ' '   ) {
				return status_string_aa_tuple{ pdb_atom_parse_status::ABORT, "Does not contain a space at column 21", the_aa };
			}
			if ( arg_pdb_atom_record_string.at( 27 ) != ' ' || arg_pdb_atom_record_string.at( 28 ) != ' ' || arg_pdb_atom_record_string.at( 29 ) != ' ' ) {
				return status_string_aa_tuple{ pdb_atom_parse_status::ABORT, "Does not contain spaces at columns 28-30", the_aa };
			}
		//	if ( arg_pdb_atom_record_string.at( 16 ) != ' ' && arg_pdb_atom_record_string.at( 16 ) != 'A' ) {
		//		return status_string_aa_tuple{ pdb_atom_parse_status::SKIP, "Has alternate location indicator other than 'A' or' '", the_aa };
		//	}
			try {
				the_aa = parse_amino_acid_from_pdb_atom_record( arg_pdb_atom_record_string );
			}
			catch (const boost::exception &e) {
				return status_string_aa_tuple{
					pdb_atom_parse_status::SKIP,
					"Do not recognise amino acid entry: \"" + arg_pdb_atom_record_string.substr( 17, 3 ) + "\" - " + diagnostic_information( e ),
					the_aa
				};
			}
			catch (...) {
				return status_string_aa_tuple{
					pdb_atom_parse_status::SKIP,
					"Do not recognise amino acid entry: " + arg_pdb_atom_record_string.substr( 17, 3 ),
					the_aa
				};
			}
			return status_string_aa_tuple{ pdb_atom_parse_status::OK, "", the_aa };
		}

		/// \brief TODOCUMENT
		///
		/// \relates pdb_atom
		inline chain_resname_atom_tuple parse_pdb_atom_record(const std::string &arg_pdb_atom_record_string, ///< TODOCUMENT
		                                                      const amino_acid  &arg_amino_acid
		                                                      ) {
			if ( isspace( arg_pdb_atom_record_string.at( 25 ) ) ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("PDB ATOM/HETATOM record malformed: space found in column 26 which should contain the end of a right justified residue number"));
			}

			try {
				                                                                                                                       // Comments with PDB format documentation
				                                                                                                                       // (http://www.wwpdb.org/documentation/format33/sect9.html#ATOM)
				const pdb_record       record_type =      pdb_rec_of_six_chars_in_string( arg_pdb_atom_record_string                ); //  1 -  6        Record name   "ATOM  "
				const size_t           serial      =  common::parse_ulong_from_substring( arg_pdb_atom_record_string,         6, 5  ); //  7 - 11        Integer       serial       Atom  serial number.
				      std::string      element     =                                      arg_pdb_atom_record_string.substr( 12, 4  ); // 13 - 16        Atom          name         Atom name.
				const char            &alt_locn    =                                      arg_pdb_atom_record_string.at(     16     ); // 17             Character     altLoc       Alternate location indicator.
		//		const std::string      the_a_a     =                                      arg_pdb_atom_record_string.substr( 17, 3  ); // 18 - 20        Residue name  resName      Residue name.
				const char            &chain_char  =                                      arg_pdb_atom_record_string.at(     21     ); // 22             Character     chainID      Chain identifier.
				const int              res_num     =    common::parse_int_from_substring( arg_pdb_atom_record_string,        22, 4  ); // 23 - 26        Integer       resSeq       Residue sequence number.
				const char            &insert_code =                                      arg_pdb_atom_record_string.at    ( 26     ); // 27             AChar         iCode        Code for insertion of residues.
				const double           coord_x     = common::parse_double_from_substring( arg_pdb_atom_record_string,        30, 8  ); // 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
				const double           coord_y     = common::parse_double_from_substring( arg_pdb_atom_record_string,        38, 8  ); // 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
				const double           coord_z     = common::parse_double_from_substring( arg_pdb_atom_record_string,        46, 8  ); // 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
				const double           occupancy   = common::parse_double_from_substring( arg_pdb_atom_record_string,        54, 6  ); // 55 - 60        Real(6.2)     occupancy    Occupancy.
				const double           temp_factor = common::parse_double_from_substring( arg_pdb_atom_record_string,        60, 6  ); // 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
		//		const std::string      element     =                           trim_copy( arg_pdb_atom_record_string.substr( 76, 3 ));  // 77 - 78        LString(2)    element      Element symbol, right-justified.
		//		const std::string      charge      =                           trim_copy( arg_pdb_atom_record_string.substr( 78, 2 ));  // 79 - 80        LString(2)    charge       Charge  on the atom.

				return std::make_tuple(
					chain_label( chain_char ),
					make_residue_name_with_non_insert_char( res_num, insert_code, ' ' ),
					pdb_atom(
						record_type,
						serial,
						std::move( element ),
						alt_locn,
						arg_amino_acid,
						geom::coord{ coord_x, coord_y, coord_z },
						occupancy,
						temp_factor
					)
				);
			}
			catch (const boost::bad_lexical_cast &) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Unable to cast a column whilst parsing a PDB ATOM record, which probably means it's malformed.\nRecord was \"" + arg_pdb_atom_record_string + "\""));
			}
			catch (const std::invalid_argument &) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Unable to cast a column whilst parsing a PDB ATOM record, which probably means it's malformed.\nRecord was \"" + arg_pdb_atom_record_string + "\""));
			}
			catch (const std::out_of_range &) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Casted column out of range whilst parsing a PDB ATOM record, which probably means it's malformed.\nRecord was \"" + arg_pdb_atom_record_string + "\""));
			}
		}

		/// \brief TODOCUMENT
		///
		/// \relates pdb_atom
		inline chain_resname_atom_tuple parse_pdb_atom_record(const std::string &arg_pdb_atom_record_string ///< TODOCUMENT
		                                                      ) {
			return parse_pdb_atom_record(
				arg_pdb_atom_record_string,
				parse_amino_acid_from_pdb_atom_record( arg_pdb_atom_record_string )
			);
		}

	} // namespace file
} // namespace cath

#endif
