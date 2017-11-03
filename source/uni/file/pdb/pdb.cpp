/// \file
/// \brief The pdb class definitions

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

#include "pdb.hpp"

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/log/trivial.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/binary_search.hpp>
#include <boost/range/algorithm/count_if.hpp>

#include "common/algorithm/copy_build.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/log/log_to_ostream_guard.hpp"
#include "common/boost_addenda/range/adaptor/adjacented.hpp"
#include "common/boost_addenda/range/back.hpp"
#include "common/boost_addenda/range/front.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/exception/runtime_error_exception.hpp"
#include "common/file/open_fstream.hpp"
#include "common/size_t_literal.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_list.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "file/pdb/protein_info.hpp"
#include "structure/geometry/coord.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"

#include <fstream>
#include <iostream>

using namespace boost::log;
using namespace boost::math::constants;
using namespace cath;
using namespace cath::chop;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;

using boost::adaptors::filtered;
using boost::adaptors::transformed;
using boost::algorithm::all;
using boost::algorithm::any_of;
using boost::algorithm::is_space;
using boost::algorithm::join;
using boost::algorithm::starts_with;
using boost::filesystem::path;
using boost::lexical_cast;
using boost::make_optional;
using boost::none;
using boost::numeric_cast;
using boost::range::binary_search;
using boost::range::count_if;
using std::get;
using std::ifstream;
using std::istream;
using std::make_pair;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::pair;
using std::right;
using std::set;
using std::setw;
using std::string;
using std::stringstream;
using std::vector;

const string pdb::PDB_RECORD_STRING_TER ( "TER   " );

/// \brief TODOCUMENT
void pdb::read_file(const path &arg_filename ///< TODOCUMENT
                    ) {
	ifstream pdb_istream;
	open_ifstream(pdb_istream, arg_filename);

	// Try here to catch any I/O exceptions
	try {
		read_pdb_file( pdb_istream, *this );

		// Close the file
		pdb_istream.close();
	}
	// Catch and immediately rethrow any boost::exceptions
	// (so that it won't get caught in the next block if it's a std::exception)
	catch (const boost::exception &ex) {
		throw;
	}
	// Catch any I/O exceptions
	catch (const std::exception &ex) {
		const string error_message("Cannot read PDB file \"" + arg_filename.string() + "\" [" + ex.what() + "] ");
		perror(error_message.c_str());
		BOOST_THROW_EXCEPTION(runtime_error_exception(error_message));
	};
}

/// \brief TODOCUMENT
void pdb::append_to_file(const path &arg_filename ///< TODOCUMENT
                         ) const {
	ofstream pdb_appstream;
	open_ofstream(pdb_appstream, arg_filename);

	// Try here to catch any I/O exceptions
	try {
		write_pdb_file(pdb_appstream, *this);
	
		// Close the file
		pdb_appstream.close();
	}
	// Catch and immediately rethrow any boost::exceptions
	// (so that it won't get caught in the next block if it's a std::exception)
	catch (const boost::exception &ex) {
		throw;
	}
	// Catch any I/O exceptions
	catch (const std::exception &ex) {
		const string error_message("Cannot append to PDB file \"" + arg_filename.string() + "\" [" + ex.what() + "] ");
		perror(error_message.c_str());
		BOOST_THROW_EXCEPTION(runtime_error_exception(error_message));
	};
}

/// \brief TODOCUMENT
pdb & pdb::set_chain_label(const chain_label &arg_chain_label ///< TODOCUMENT
                           ) {
	for (pdb_residue &my_pdb_residue : pdb_residues) {
		my_pdb_residue.set_chain_label( arg_chain_label );
	}
	return *this;
}

/// \brief TODOCUMENT
residue_id_vec pdb::get_residue_ids_of_first_chain__backbone_unchecked() const {
	return get_backbone_complete_residue_ids_of_first_chain( *this, false );
}

/// \brief TODOCUMENT
coord pdb::get_residue_ca_coord_of_index__backbone_unchecked(const size_t &arg_index ///< TODOCUMENT
                                                             ) const {
	return get_carbon_alpha_coord( get_residue_of_index__backbone_unchecked( arg_index ) );
}

/// \brief TODOCUMENT
size_t pdb::get_num_atoms() const {
	auto num_atoms = 0_z;
	for (const pdb_residue &residue : pdb_residues) {
		num_atoms += residue.get_num_atoms();
	}
	return num_atoms;
}

/// \brief TODOCUMENT
pdb & pdb::rotate(const rotation &arg_rotation ///< TODOCUMENT
                  ) {
	for (pdb_residue &my_pdb_residue : pdb_residues) {
		my_pdb_residue.rotate( arg_rotation );
	}
	for (pdb_residue &my_pdb_residue : post_ter_residues) {
		my_pdb_residue.rotate( arg_rotation );
	}
	return *this;
}

/// \brief TODOCUMENT
pdb & pdb::operator+=(const coord &arg_coord ///< TODOCUMENT
                      ) {
	for (pdb_residue &my_pdb_residue : pdb_residues) {
		my_pdb_residue += arg_coord;
	}
	for (pdb_residue &my_pdb_residue : post_ter_residues) {
		my_pdb_residue += arg_coord;
	}
	return *this;
}

/// \brief TODOCUMENT
pdb & pdb::operator-=(const coord &arg_coord ///< TODOCUMENT
                      ) {
	for (pdb_residue &my_pdb_residue : pdb_residues) {
		my_pdb_residue -= arg_coord;
	}
	for (pdb_residue &my_pdb_residue : post_ter_residues) {
		my_pdb_residue -= arg_coord;
	}
	return *this;
}

/// \brief TODOCUMENT
pdb & pdb::set_residues(pdb_residue_vec arg_pdb_residues ///< TODOCUMENT
                       ) {
	pdb_residues = std::move( arg_pdb_residues );
	return *this;
}

/// \brief Rvalue-reference setter for the residues that appear after a TER record in their respective chains
pdb & pdb::set_post_ter_residues(pdb_residue_vec arg_post_ter_residues ///< The residues that appear after a TER record in their respective chains
                                ) {
	post_ter_residues = std::move( arg_post_ter_residues );
	return *this;
}

/// \brief Getter for the residues that appear after a TER record in their respective chains
const pdb_residue_vec & pdb::get_post_ter_residues() const {
	return post_ter_residues;
}

/// \brief TODOCUMENT
///
/// \relates pdb
size_t cath::file::get_num_backbone_complete_residues(const pdb &arg_pdb ///< The pdb to query
                                                      ) {
	return numeric_cast<size_t>( count_if(
		arg_pdb,
		[] (const pdb_residue &x) { return is_backbone_complete( x ); }
	) );
}

/// \brief TODOCUMENT
///
/// \relates pdb
size_t cath::file::get_index_of_backbone_complete_index(const pdb    &arg_pdb,  ///< The pdb to query
                                                        const size_t &arg_index ///< The index of the required residue
                                                        ) {
	using std::to_string;

	if ( arg_index >= arg_pdb.get_num_residues() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to get_residue_ca_coord_of_backbone_complete_index() for index >= number of residues"));
	}
	size_t count = 0;
	for (const size_t &index : indices( arg_pdb.get_num_residues() ) ) {
		const pdb_residue &the_res = arg_pdb.get_residue_of_index__backbone_unchecked( index );
		if ( is_backbone_complete( the_res ) ) {
			if ( count == arg_index ) {
				return index;
			}
			++count;
		}
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception(
		"Cannot find enough backbone_complete residues to reach backbone_complete_index "
		+ to_string( arg_index )
	));
	return 0; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
}

/// \brief TODOCUMENT
///
/// \relates pdb
const pdb_residue & cath::file::get_residue_of_backbone_complete_index(const pdb    &arg_pdb,  ///< The pdb to query
                                                                       const size_t &arg_index ///< The index of the required residue
                                                                       ) {
	return arg_pdb.get_residue_of_index__backbone_unchecked(
		get_index_of_backbone_complete_index(
			arg_pdb,
			arg_index
		)
	);
}

/// \brief TODOCUMENT
///
/// \relates pdb
coord cath::file::get_residue_ca_coord_of_backbone_complete_index(const pdb    &arg_pdb,  ///< The pdb to query
                                                                  const size_t &arg_index ///< The index of the required residue
                                                                  ) {
	return arg_pdb.get_residue_ca_coord_of_index__backbone_unchecked(
		get_index_of_backbone_complete_index(
			arg_pdb,
			arg_index
		)
	);
}
/// \brief TODOCUMENT
///
/// \relates pdb
size_t cath::file::get_num_region_limited_backbone_complete_residues(const pdb             &arg_pdb,    ///< The pdb to query
                                                                     const region_vec_opt  &arg_regions ///< The regions within which the count applies
                                                                     ) {
	regions_limiter limiter{ arg_regions };
	size_t count = 0;
	for (const pdb_residue &the_res : arg_pdb) {
		if ( limiter.update_residue_is_included( the_res.get_residue_id() ) && is_backbone_complete( the_res ) ) {
			++count;
		}
	}
	return count;
}

/// \brief TODOCUMENT
///
/// \relates pdb
size_t cath::file::get_index_of_region_limited_backbone_complete_index(const pdb             &arg_pdb,    ///< The pdb to query
                                                                       const size_t          &arg_index,  ///< The index of the required residue
                                                                       const region_vec_opt  &arg_regions ///< The regions within which the index applies
                                                                       ) {
	using std::to_string;

	if ( arg_index >= arg_pdb.get_num_residues() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to get_index_of_region_limited_backbone_complete_index() for index >= number of residues"));
	}

	regions_limiter limiter{ arg_regions };

	size_t count = 0;
	for (const size_t &index : indices( arg_pdb.get_num_residues() ) ) {
		const pdb_residue &the_res = arg_pdb.get_residue_of_index__backbone_unchecked( index );
		if ( limiter.update_residue_is_included( the_res.get_residue_id() ) && is_backbone_complete( the_res ) ) {
			if ( count == arg_index ) {
				return index;
			}
			++count;
		}
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception(
		"Cannot find enough backbone_complete residues to reach region_limited_backbone_complete_index "
		+ to_string( arg_index )
	));
	return 0; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
}

/// \brief TODOCUMENT
///
/// \relates pdb
const pdb_residue & cath::file::get_residue_of_region_limited_backbone_complete_index(const pdb             &arg_pdb,    ///< The pdb to query
                                                                                      const size_t          &arg_index,  ///< The index of the required residue
                                                                                      const region_vec_opt  &arg_regions ///< The regions within which the index applies
                                                                                      ) {
	return arg_pdb.get_residue_of_index__backbone_unchecked(
		get_index_of_region_limited_backbone_complete_index(
			arg_pdb,
			arg_index,
			arg_regions
		)
	);
}

/// \brief TODOCUMENT
///
/// \relates pdb
coord cath::file::get_residue_ca_coord_of_region_limited_backbone_complete_index(const pdb             &arg_pdb,    ///< The pdb to query
                                                                                 const size_t          &arg_index,  ///< The index of the required residue
                                                                                 const region_vec_opt  &arg_regions ///< The regions within which the index applies
                                                                                 ) {
	return arg_pdb.get_residue_ca_coord_of_index__backbone_unchecked(
		get_index_of_region_limited_backbone_complete_index(
			arg_pdb,
			arg_index,
			arg_regions
		)
	);
}

/// \brief Get the list of residues IDs for the backbone-complete residues on first chain of the specified PDB
///
/// \relates pdb
residue_id_vec cath::file::get_backbone_complete_residue_ids_of_first_chain(const pdb  &arg_pdb,                   ///< The PDB to query
                                                                            const bool &arg_complete_backbone_only ///< Whether to restrict to the backbone-complete residues
                                                                            ) {
	if ( arg_pdb.empty() ) {
		return {};
	}

	const chain_label first_chain_label = get_chain_label( front( arg_pdb ) );
	return transform_build<residue_id_vec>(
		arg_pdb
			| filtered( [&] (const pdb_residue &x) {
				return (
					( get_chain_label( x ) == first_chain_label )
					&&
					( ! arg_complete_backbone_only || is_backbone_complete( x ) )
				);
			} ),
		[] (const pdb_residue &x) { return x.get_residue_id(); }
	);
}

/// \brief Get the list of residues IDs for the backbone-complete residues on all chains of the specified PDB
///
/// \relates pdb
residue_id_vec cath::file::get_backbone_complete_residue_ids(const pdb &arg_pdb ///< The PDB to query
                                                             ) {
	return arg_pdb.empty()
		? residue_id_vec{}
		: transform_build<residue_id_vec>(
			arg_pdb | filtered( is_backbone_complete ),
			[] (const pdb_residue &x) { return x.get_residue_id(); }
		);
}

/// \brief TODOCUMENT
///
/// \relates pdb
pdb cath::file::read_pdb_file(const path &arg_pdb_filename ///< TODOCUMENT
                              ) {
	pdb new_pdb;
	new_pdb.read_file( arg_pdb_filename );
	return new_pdb;
}

/// \brief TODOCUMENT
///
/// \relates pdb
pdb cath::file::read_pdb_file(istream &arg_input_stream ///< TODOCUMENT
                              ) {
	pdb new_pdb;
	read_pdb_file( arg_input_stream, new_pdb );
	return new_pdb;
}

/// \brief TODOCUMENT
///
/// \relates pdb
istream & cath::file::read_pdb_file(istream &input_stream, ///< TODOCUMENT
                                    pdb     &arg_pdb       ///< TODOCUMENT
                                    ) {
	// Variables to store the details of parsed atoms
	pdb_residue_vec  residues;
	pdb_residue_vec  post_ter_residues;

	set<chain_label> terminated_chains;
	string           line_string;

	char_3_arr_opt   prev_amino_acid_3_char_code;
	pdb_atom_vec     prev_atoms;
	residue_id       prev_res_id;
	bool             prev_warned_conflict = false;

	const auto add_atoms_and_reset_fn = [&] (const chain_label &x) {
		// if ( ! prev_atoms.empty() ) {
		// 	prev_atoms.last()
		// }
		pdb_residue_vec &write_residues = contains( terminated_chains, x ) ? post_ter_residues
		                                                                   : residues;
		write_residues.emplace_back(
			prev_res_id,
			std::move( prev_atoms )
		);
		prev_amino_acid_3_char_code = none;
		prev_atoms                  = pdb_atom_vec{};
		prev_warned_conflict        = false;
	};

	// Loop over the lines of the file
	//
	// This code is made a bit more complicated because the aim is to
	// add all of the atoms within a residue at the same time but it isn't
	// clear that the residue has finished until the first line of the next residue
	// (or the end of the file)
	while ( getline( input_stream, line_string ) ) {
		// If this line is an ATOM or HETATM record
		if ( is_pdb_record_of_type( line_string, pdb_record::ATOM ) || is_pdb_record_of_type( line_string, pdb_record::HETATM ) ) {
			const bool is_atom                 = is_pdb_record_of_type( line_string, pdb_record::ATOM );
			const auto parse_status_str_and_aa = pdb_record_parse_problem( line_string );
			const auto &parse_status = get<0>( parse_status_str_and_aa );
			const auto &parse_string = get<1>( parse_status_str_and_aa );
			const auto &parse_aa     = get<2>( parse_status_str_and_aa );
			if ( parse_status == pdb_atom_parse_status::ABORT ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception(
					"ATOM record is malformed : " + parse_string
					+ "\nRecord was \"" + line_string.substr(0, pdb_atom::MAX_NUM_PDB_COLS)
					+ "\""
				));
			}
			else if ( parse_status == pdb_atom_parse_status::SKIP ) {
				BOOST_LOG_TRIVIAL( warning ) << "Skipping PDB atom record \""
					<< line_string
					<< "\" with message: "
					<< parse_string;
				continue;
			}

			// Grab the details from parsing this ATOM record
			const auto        new_entry              = parse_pdb_atom_record( line_string, parse_aa );
			const residue_id &res_id                 = new_entry.first;
			const pdb_atom   &atom                   = new_entry.second;
			const char_3_arr  amino_acid_3_char_code = get_amino_acid_code( atom );

			// If there are previously seen atoms that don't match this chain/res_id,
			// then add those atoms' residue and reset prev_atoms
			if ( ! prev_atoms.empty() && res_id != prev_res_id ) {
				add_atoms_and_reset_fn( prev_res_id.get_chain_label() );
			}

			// Some PDBs (eg 4tsw) may have erroneous consecutive duplicate residues.
			// Though that's a bit rubbish, it shouldn't break the whole comparison
			// so if that's detected, just warn and move on (without appending to new_residues).
			if (
				is_atom
				&&
				! prev_atoms.empty()
				&&
				prev_amino_acid_3_char_code
				&&
				amino_acid_3_char_code != prev_amino_acid_3_char_code
				&&
				atom.get_alt_locn() == ' '
				) {
				if ( ! prev_warned_conflict ) {
					BOOST_LOG_TRIVIAL( warning ) << "Whilst parsing PDB file, found conflicting consecutive entries for residue \""
					                             << res_id
					                             << "\" (with amino acids \""
					                             << char_arr_to_string( *prev_amino_acid_3_char_code )
					                             << "\" and then \""
					                             << char_arr_to_string( amino_acid_3_char_code )
					                             << "\") - ignoring latter entry (and any further entries)";
					prev_warned_conflict = true;
				}
			}
			// Otherwise update the records of previously seen atoms
			else {
				if ( is_atom ) {
					prev_amino_acid_3_char_code = amino_acid_3_char_code;
				}
				prev_res_id = res_id;
				prev_atoms.push_back( atom );
			}
		}
		else if ( boost::algorithm::starts_with( line_string, "ENDMDL" ) ) {
			break;
		}
		else if ( boost::algorithm::starts_with( line_string, "TER" ) ) {
			if ( ! prev_atoms.empty() ) {
				add_atoms_and_reset_fn( prev_res_id.get_chain_label() );
			}
			if ( line_string.length() >= 22 ) {
				terminated_chains.insert( chain_label( line_string.at( 21 ) ) );
			}
			else if ( ! is_null( prev_res_id ) ) {
				terminated_chains.insert( prev_res_id.get_chain_label() );
			}
		}
	};

	// Add any last remaining atoms
	if ( ! prev_atoms.empty() ) {
		add_atoms_and_reset_fn( prev_res_id.get_chain_label() );
	}

	arg_pdb.set_residues         ( std::move( residues          ) );
	arg_pdb.set_post_ter_residues( std::move( post_ter_residues ) );

	return input_stream;
}

/// \brief TODOCUMENT
///
/// \relates pdb
pdb_list cath::file::read_end_separated_pdb_files(istream &arg_in_stream ///< TODOCUMENT
                                                  ) {
	pdb_list pdbs;
	stringstream pdb_file_stream;
	string line_str;
	while (getline(arg_in_stream, line_str)) {
		// If this line begins with end, then add any PDB that's been accumulated in pdb_file_stream
		// and reset the pdb_file_stream
		if (starts_with(line_str, "END")) {
			// If there is something other than whitespace in pdb_file_stream, then
			// try to parse it into a PDB
			if (!all(pdb_file_stream.str(), is_space())) {
				pdbs.push_back(read_pdb_file(pdb_file_stream));
			}

			// Reset pdb_file_stream
			pdb_file_stream.str(string());
			pdb_file_stream.clear();
		}
		// Otherwise, just add this line to pdb_file_stream
		else {
			pdb_file_stream << line_str << "\n";
		}
	}

	// If there's anything other than whitespace left in pdb_file_stream, then
	// try to parse it into a PDB
	if (!all(pdb_file_stream.str(), is_space())) {
		pdbs.push_back(read_pdb_file(pdb_file_stream));
	}

	// Return any PDBs that have been parsed
	return pdbs;
}

/// \brief Write the specified PDB to a pdb file string,
///        restricted to the specified regions and written in the specified mode
///
/// \relates pdb
string cath::file::to_pdb_file_string(const pdb             &arg_pdb,             ///< The pdb to describe
                                      const regions_limiter &arg_regions_limiter, ///< Optional specification of regions to which the written records should be restricted
                                      const pdb_write_mode  &arg_pdb_write_mode   ///< Whether this is the only/last part of the PDB file
                                      ) {
	ostringstream output_ss;
	write_pdb_file(
		output_ss,
		arg_pdb,
		arg_regions_limiter,
		arg_pdb_write_mode
	);
	return output_ss.str();
}

/// \brief Insert a PDB file of the specified pdb into the specified ostream,
///        restricted to the specified regions and written in the specified mode
///
/// \relates pdb
ostream & cath::file::write_pdb_file(ostream               &arg_os,              ///< The ostream into which the PDB file should be inserted
                                     const pdb             &arg_pdb,             ///< The pdb to describe
                                     const regions_limiter &arg_regions_limiter, ///< Optional specification of regions to which the written records should be restricted
                                     const pdb_write_mode  &arg_pdb_write_mode   ///< Whether this is the only/last part of the PDB file
                                     ) {
	regions_limiter the_regions_limiter{ arg_regions_limiter };
	const auto &num_residues = arg_pdb.get_num_residues();
//	size_t atom_ctr = 1;
	for (const size_t &residue_ctr : indices( num_residues ) ) {
		const pdb_residue &the_residue = arg_pdb.get_residue_of_index__backbone_unchecked( residue_ctr );
		if ( the_regions_limiter.update_residue_is_included( the_residue.get_residue_id() ) ) {
			write_pdb_file_entry( arg_os, the_residue );
		}

		const auto &the_chain_label = get_chain_label( the_residue );

		if ( residue_ctr + 1 == num_residues
		     || get_chain_label( arg_pdb.get_residue_of_index__backbone_unchecked( residue_ctr + 1 ) ) != the_chain_label
		     ) {
			const auto last_atom = the_residue.empty() ? none : make_optional( back( the_residue ) );
			const string residue_name_with_insert_or_space = make_residue_name_string_with_insert_or_space(
				get_residue_name( the_residue )
			);
			arg_os << pdb::PDB_RECORD_STRING_TER
				<< right
				<< setw( 5 ) << ( last_atom ? ( last_atom->get_atom_serial() + 1u ) : 0u )
				<< "      "
				<<              ( last_atom ?   last_atom->get_amino_acid() : amino_acid{ 'X' } )
				<< " "
				<< the_chain_label
				<< setw( 5 ) << residue_name_with_insert_or_space
				<< "                                                     \n";
		}
	}

	for (const pdb_residue &the_residue : arg_pdb.get_post_ter_residues() ) {
		write_pdb_file_entry( arg_os, the_residue );
	}

	// If this is the only or last PDB then "END   " the file
	if ( arg_pdb_write_mode == pdb_write_mode::ONLY_OR_LAST_PDB ) {
		arg_os << "END   \n";
	}

	return arg_os;
}

/// \brief Write a PDB file of the specified pdb to the specified file,
///        restricted to the specified regions and written in the specified mode
///
/// \relates pdb
void cath::file::write_pdb_file(const path            &arg_filename,        ///< The file to which the PDB file should be written
                                const pdb             &arg_pdb,             ///< The pdb to describe
                                const regions_limiter &arg_regions_limiter, ///< Optional specification of regions to which the written records should be restricted
                                const pdb_write_mode  &arg_pdb_write_mode   ///< Whether this is the only/last part of the PDB file
                                ) {
	ofstream out_ofstream;
	open_ofstream( out_ofstream, arg_filename );
	write_pdb_file(
		out_ofstream,
		arg_pdb,
		arg_regions_limiter,
		arg_pdb_write_mode
	);
	out_ofstream.close();
}

/// \brief TODOCUMENT
///
/// \relates pdb
amino_acid_vec cath::file::get_amino_acid_list(const pdb &arg_pdb ///< TODOCUMENT
                                               ) {
	amino_acid_vec amino_acids;
	amino_acids.reserve( arg_pdb.get_num_residues() );
	for (const pdb_residue &the_pdb_residue : arg_pdb) {
		amino_acids.push_back( the_pdb_residue.get_amino_acid() );
	}
	return amino_acids;
}

/// \brief Calculate the (ascending) list of indices of the residues that follow a chain break in the specified PDB
///
/// \relates pdb
size_vec cath::file::indices_of_residues_following_chain_breaks(const pdb &arg_pdb /// The PDB to query
                                                                ) {
	// The gap between consecutive residues' carbon and nitrogen atoms respectively,
	// above which the residues are not treated as neighbours
	constexpr double INTER_C_TO_N_DIST_FOR_NEIGHBOURS = 2.5;

	return copy_build<size_vec>(
		indices( arg_pdb.get_num_residues() )
			| filtered( [&] (const size_t &x) {
				if ( x == 0 ) {
					return false;
				}
				const auto &prev_res = arg_pdb.get_residue_of_index__backbone_unchecked( x - 1 );
				const auto &this_res = arg_pdb.get_residue_of_index__backbone_unchecked( x     );
				return (
					( get_chain_label( prev_res ) != get_chain_label( this_res ) )
					||
					(
						distance_between_points( get_carbon_coord( prev_res ), get_nitrogen_coord( this_res ) )
						>
						INTER_C_TO_N_DIST_FOR_NEIGHBOURS
					)
				);
			} )
	);
}

/// \brief TODOCUMENT
///
/// \relates pdb
ostream & cath::file::operator<<(ostream   &arg_os,         ///< TODOCUMENT
                                 const pdb &arg_pdb_residue ///< TODOCUMENT
                                 ) {
	arg_os << "PDB[";

	const size_t num_residues = arg_pdb_residue.get_num_residues();
	for (const size_t &residue_ctr : indices( num_residues ) ) {
		arg_os << arg_pdb_residue.get_residue_of_index__backbone_unchecked(residue_ctr);
	}
	return arg_os;
}

///// \brief TODOCUMENT
/////
///// \relates pdb
//const pdb_residue & cath::get_residue_ref_of_index__offset_1(const pdb          &arg_pdb,  ///< TODOCUMENT
//                                                                   const size_t &arg_index ///< TODOCUMENT
//                                                                   ) {
//	if (arg_index <= 0) {
//		BOOST_THROW_EXCEPTION(invalid_argument_exception("An index that uses offset 1 must be 1 or greater"));
//	}
//	return arg_pdb.get_residue_of_index__backbone_unchecked(arg_index - 1);
//}

/// \brief Get the list of phi/psi angle pairs for each residue, in radians within (0, 2 * pi]
///
/// There are two different possible causes for potentially skipping residues
/// ( and so possibly causing breaks in the phi/psi angles):
///  * the residue may not have a complete set of backbone atoms (N, CA, C, O)
///  * the residue may be complete but its alt-locn specifications may not be to DSSP's liking
///
/// In the first case, we would expect the residues to already have been stripped out before
/// this code is called. But not necessarily in the latter case.
///
/// It turns out that DSSP just uses distance from one residue's C to the next residue's N
/// to check for chain breaks. In some cases, that can mean that it doesn't assign a break
/// even when that means calculating phi/psi angles between residues that just happen to
/// be close enough even though they straddle a residue that DSSP has rejected
/// (see https://github.com/cmbi/xssp/issues/86). That functionality has been duplicated here.
///
/// \todo Depending on the outcome of https://github.com/cmbi/xssp/issues/86, either
///       change this code to do smarter break detection (see commented parts) or
///       remove that commented code and remove the second indices_in_new_of_skips part
///       of the return type of backbone_complete_subset_of_pdb().
///
/// At present, undetermined angles are set to 2 * pi
///
/// \todo Shouldn't the undetermined angles be handled more explicitly?
///       The DSSP/WOLF files have the angles in the more natural range of [-180, 180]
///       and the currently always get shifted into the (0, 360] range
///       if they were left alone, that would leave 360.0 as a special "undetermined" value
///       but this would require changes in residues_have_similar_area_angle_props()
///
/// \relates pdb
///
/// \relates protein
doub_angle_doub_angle_pair_vec cath::file::get_phi_and_psi_angles(const pdb                      &arg_pdb,                ///< TODOCUMENT
                                                                  // const size_vec                 &/*arg_skip_indices*/,   ///< Indices of residues in the pdb that were preceded by residues that have been skipped due to being backbone complete. This may include an index on greater than the index of the last residue in the pdb to indicate that there were residues skipped after the last residue. Phi/psi angles are not set over these skip breaks.
                                                                  const dssp_skip_angle_skipping &arg_dssp_angle_skipping ///< TODOCUMENT
                                                                  ) {
	// The gap between consecutive residues' carbon and nitrogen atoms respectively,
	// above which the residues are not treated as neighbours
	constexpr double INTER_C_TO_N_DIST_FOR_NEIGHBOURS = 2.5;
	
	const     auto   DEFAULT_PHI_PSI              = residue::DEFAULT_PHI_PSI();
	const     auto   DEFAULT_PHI_PSI_PAIR         = make_pair( DEFAULT_PHI_PSI, DEFAULT_PHI_PSI );
	const     size_t num_residues                 = arg_pdb.get_num_residues();

	// Build a range of indices of the residues to be considered
	//
	// It isn't good enough to check pairs of consecutive indices because
	// this should replicate DSSP and that requires checking between pairs of
	// consecutive *non-skipped* residues (straddling any gaps in indices where necessary)
	const auto non_skipped_residues_indices = indices( num_residues )
		| filtered( [&] (const size_t &x) {
			return (
				( arg_dssp_angle_skipping == dssp_skip_angle_skipping::DONT_BREAK_ANGLES )
				||
				! dssp_will_skip_residue( arg_pdb.get_residue_of_index__backbone_unchecked( x ) )
			);
		} );

	// Prepare a data structure of phi/psi angles
	doub_angle_doub_angle_pair_vec phi_and_psi_angles( num_residues, DEFAULT_PHI_PSI_PAIR );

	// Loop over pairs of consecutive *non-skipped* residues
	for (const size_size_pair &adj_indices : non_skipped_residues_indices | adjacented) {

		// Grab this residue and the next non-skipped one
		const pdb_residue &this_pdb_residue = arg_pdb.get_residue_of_index__backbone_unchecked( adj_indices.first  );
		const pdb_residue &next_pdb_residue = arg_pdb.get_residue_of_index__backbone_unchecked( adj_indices.second );

		// If these consecutive residues are on the same chain...
		if ( get_chain_label( this_pdb_residue ) == get_chain_label( next_pdb_residue ) ) {

			// ...and if they're adequately close together...
			const double inter_c_to_n_dist = distance_between_points(
				get_carbon_coord  ( this_pdb_residue ),
				get_nitrogen_coord( next_pdb_residue )
			);

			// If:
			//  * they're close enough and
			//  * there weren't any residues between them that have been skipped and
			//  * they're not to be skipped due to either being residues DSSP would skip
			const bool close_enough           = inter_c_to_n_dist <= INTER_C_TO_N_DIST_FOR_NEIGHBOURS;

			// const bool not_straddling_skipped = ! binary_search( arg_skip_indices, residue_ctr + 1 );
			// const bool not_dssp_skip          = (
			// 	( arg_dssp_angle_skipping == dssp_skip_angle_skipping::DONT_BREAK_ANGLES )
			// 	||
			// 	! ( dssp_will_skip_residue( this_pdb_residue ) || dssp_will_skip_residue( next_pdb_residue ) )
			// );
			// if ( close_enough && not_straddling_skipped && not_dssp_skip ) {
			if ( close_enough ) {
				// Then calculate the two angles of these two residue and store them
				const auto this_psi_and_next_phi = get_psi_of_this_and_phi_of_next( this_pdb_residue, next_pdb_residue );
				phi_and_psi_angles[ adj_indices.first  ].second = this_psi_and_next_phi.first;
				phi_and_psi_angles[ adj_indices.second ].first  = this_psi_and_next_phi.second;
			}
		}
	}

	// Return the calculated angles
	return phi_and_psi_angles;
}

/// \brief TODOCUMENT
///
/// Returns the backbone complete subset of the PDB along with indices of residues in that pdb that
/// were preceded by residues that have been skipped due to being backbone complete.
/// This may include an index on greater than the index of the last residue in the pdb to indicate
/// that there were residues skipped after the last residue.
/// This information is useful to return so it can be used to prevent phi/psi angles being set over these skip breaks.
///
/// \relates pdb
pdb_size_vec_pair cath::file::backbone_complete_subset_of_pdb(const pdb                    &arg_pdb,             ///< TODOCUMENT
                                                              const ostream_ref_opt        &arg_ostream_ref_opt, ///< An optional reference to an ostream to which any logging should be sent
                                                              const dssp_skip_res_skipping &arg_skip_like_dssp   ///< TODOCUMENT
                                                              ) {
	// Grab the number of residues
	const size_t num_residues = arg_pdb.get_num_residues();

	vector<residue_id> seen_residue_ids;

	residue_id_vec backbone_skipped_residues;
	size_vec indices_in_new_of_skips;

	// Prepare a vector of the new residues
	pdb_residue_vec new_pdb_residues;
	new_pdb_residues.reserve(num_residues);

	// Loop over the residues in the input pdb
	for (const size_t &residue_ctr : indices( num_residues ) ) {
		const pdb_residue &the_residue = arg_pdb.get_residue_of_index__backbone_unchecked( residue_ctr );
		const bool seen_res_id = common::contains( seen_residue_ids, the_residue.get_residue_id() );
		if ( ! seen_res_id ) {
			seen_residue_ids.push_back( the_residue.get_residue_id() );
		}

		// If the residue is backbone_complete,then add it to new_pdb_residues
		const bool ok_to_process = ( arg_skip_like_dssp == dssp_skip_res_skipping::SKIP )
			? ! dssp_will_skip_residue( the_residue )
			:   is_backbone_complete  ( the_residue );
		if ( ok_to_process && ! seen_res_id ) {
			new_pdb_residues.push_back( the_residue );
		}
		// Else if this is a proper amino acid (not just a bunch of HETATMs), record it
		else if ( get_letter_if_amino_acid( the_residue ) ) {
			if ( arg_ostream_ref_opt ) {
				backbone_skipped_residues.push_back( the_residue.get_residue_id() );
			}
			if ( indices_in_new_of_skips.empty() || indices_in_new_of_skips.back() != new_pdb_residues.size() ) {
				indices_in_new_of_skips.push_back( new_pdb_residues.size() );
			}
		}
	}

	if ( ! backbone_skipped_residues.empty() && arg_ostream_ref_opt ) {
		const bool multiple_skippeds = ( backbone_skipped_residues.size() > 1 );
		const log_to_ostream_guard ostream_log_guard{ arg_ostream_ref_opt.get().get() };

		BOOST_LOG_TRIVIAL( warning ) << "Ignoring residue"
		                             << ( multiple_skippeds ? "s"s : ""s )
		                             << " "
		                             << join(
		                             	backbone_skipped_residues
		                             		| transformed( [] (const residue_id &x) { return to_string( x ); } ),
		                             	", "
		                             )
		                             << " whilst extracting a protein structure from PDB file data because "
		                             << ( multiple_skippeds ? "they don't"s : "it doesn't"s )
		                             << " have all of N, CA and C atoms (or because of duplicate residue IDs)";
	}

	// Return a new pdb containing these residues
	pdb new_pdb;
	new_pdb.set_residues( new_pdb_residues );
	return { new_pdb, indices_in_new_of_skips };
}

/// \brief TODOCUMENT
///
/// \relates pdb
///
/// \relates protein
pair<protein, protein_info> cath::file::build_protein_of_pdb(const pdb              &arg_pdb,        ///< TODOCUMENT
                                                             const ostream_ref_opt  &arg_ostream,    ///< An optional reference to an ostream to which any logging should be sent
                                                             const dssp_skip_policy &arg_skip_policy ///< TODOCUMENT
                                                             ) {
	constexpr size_t DEFAULT_ACCESSIBILITY = 0;

	const auto     backbone_complete_data       = backbone_complete_subset_of_pdb( arg_pdb, arg_ostream, res_skipping_of_dssp_skip_policy( arg_skip_policy ) );
	const pdb      backbone_complete_pdb_subset = backbone_complete_data.first;
	// const size_vec indices_of_skips             = backbone_complete_data.second;
	const size_t   num_residues                 = backbone_complete_pdb_subset.get_num_residues();
	const auto     phi_and_psi_angles           = get_phi_and_psi_angles(
		backbone_complete_pdb_subset,
		// indices_of_skips,
		angle_skipping_of_dssp_skip_policy( arg_skip_policy )
	);

	return {
		build_protein( transform_build<residue_vec>(
			indices( num_residues ),
			[&] (const size_t &x) {
				const pdb_residue &the_residue = backbone_complete_pdb_subset.get_residue_of_index__backbone_unchecked( x );
				const auto        &phi         = phi_and_psi_angles[ x ].first;
				const auto        &psi         = phi_and_psi_angles[ x ].second;
				return build_residue_of_pdb_residue(
					the_residue,
					phi,
					psi,
					DEFAULT_ACCESSIBILITY
				);
			}
		) ),
		protein_info{
			transform_build<residue_makeup_vec>(
				backbone_complete_pdb_subset,
				[] (const pdb_residue &x) {
					return contains_non_proper_amino_acids( x );
				}
			)
		}
	};
}

/// \brief TODOCUMENT
///
/// \relates pdb
///
/// \relates protein
protein cath::file::build_protein_of_pdb_and_name(const pdb             &arg_pdb,    ///< TODOCUMENT
                                                  const name_set        &arg_name,   ///< TODOCUMENT
                                                  const ostream_ref_opt &arg_ostream ///< An optional reference to an ostream to which any logging should be sent
                                                  ) {
	protein new_protein = build_protein_of_pdb( arg_pdb, arg_ostream ).first;
	new_protein.set_name_set( arg_name );
	return new_protein;
}

/// \brief Generate a list of protein residue indices (corresponding to those returned by build_protein_of_pdb())
///        that DSSP might be expected to skip
///
/// \relates pdb
size_set cath::file::get_protein_res_indices_that_dssp_might_skip(const pdb             &arg_pdb,    ///< The PDB to query
                                                                  const ostream_ref_opt &arg_ostream ///< An ostream to which any status messages might be sent
                                                                  ) {
	const pdb    backbone_complete_pdb_subset = backbone_complete_subset_of_pdb( arg_pdb, arg_ostream ).first;
	const size_t num_residues                 = backbone_complete_pdb_subset.get_num_residues();

	// Return the indices corresponding to residues with any atoms with non-standard alt_locn values
	return copy_build<size_set>(
		indices( num_residues )
			| filtered( [&] (const size_t &x) {
				return dssp_will_skip_residue(
					backbone_complete_pdb_subset.get_residue_of_index__backbone_unchecked( x )
				);
			} )
	);
}

/// \brief Get a copy of the specified PDB restricted to the specified regions
///
/// This does a similar job to get_supn_content_pdb() but that also handles a
/// superposition_content_spec parameter.
///
/// \relates pdb
pdb cath::file::get_regions_limited_pdb(const region_vec_opt &arg_regions, ///< The regions to which the resulting PDB should be restricted
                                        const pdb            &arg_pdb      ///< The source PDB
                                        ) {
	if ( ! arg_regions ) {
		return arg_pdb;
	}

	// *arg_regions is guaranteed to outlive this regions_limiter
	regions_limiter the_limiter{ *arg_regions };
	pdb_residue_vec residues;
	residues.reserve( arg_pdb.get_num_residues() );
	for (const pdb_residue &the_residue : arg_pdb) {
		if ( the_limiter.update_residue_is_included( the_residue.get_residue_id() ) ) {
			residues.push_back( the_residue );
		}
	}

	warn_if_specified_regions_remain_unseen( the_limiter );

	pdb result_pdb;
	result_pdb.set_residues( std::move( residues ) );
	return result_pdb;
}

/// \brief TODOCUMENT
///
/// This could probably be made substantially more efficient by unifying the steps
/// but that would likely make things quite a bit more complicated.
/// If that does turn out to be worthwhile, ensure that the code still allows
/// for regions that start/stop on backbone-incomplete residues.
///
/// \relates pdb
pdb cath::file::backbone_complete_region_limited_subset_of_pdb(const pdb             &arg_pdb,     ///< TODOCUMENT
                                                               const region_vec_opt  &arg_regions, ///< TODOCUMENT
                                                               const ostream_ref_opt &arg_ostream  ///< An optional reference to an ostream to which any logging should be sent
                                                               ) {
	return backbone_complete_subset_of_pdb(
		get_regions_limited_pdb(
			arg_regions,
			arg_pdb
		),
		arg_ostream
	).first;
}
