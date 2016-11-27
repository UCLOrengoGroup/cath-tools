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

#include "pdb.h"

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/log/trivial.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/count_if.hpp>
#include <boost/range/irange.hpp>

#include "common/algorithm/copy_build.h"
#include "common/algorithm/transform_build.h"
#include "common/boost_addenda/log/log_to_ostream_guard.h"
#include "common/cpp14/cbegin_cend.h"
#include "common/file/open_fstream.h"
#include "common/size_t_literal.h"
#include "exception/invalid_argument_exception.h"
#include "exception/not_implemented_exception.h" // ***** TEMPORARY *****
#include "exception/runtime_error_exception.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_list.h"
#include "file/pdb/pdb_residue.h"
#include "structure/geometry/coord.h"
#include "structure/protein/protein.h"
#include "structure/protein/residue.h"
#include "structure/protein/sec_struc.h"
#include "structure/protein/sec_struc_planar_angles.h"

#include <fstream>
#include <iostream>

using namespace boost::filesystem;
using namespace boost::log;
using namespace boost::math::constants;
using namespace cath;
using namespace cath::chop;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace std;

using boost::adaptors::filtered;
using boost::adaptors::transformed;
using boost::algorithm::all;
using boost::algorithm::any_of;
using boost::algorithm::is_space;
using boost::algorithm::join;
using boost::algorithm::starts_with;
using boost::irange;
using boost::lexical_cast;
using boost::numeric_cast;
using boost::range::count_if;

/// \brief TODOCUMENT
void pdb::do_read_file(const path &arg_filename ///< TODOCUMENT
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
void pdb::do_append_to_file(const path &arg_filename ///< TODOCUMENT
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
void pdb::do_set_chain_label(const chain_label &arg_chain_label ///< TODOCUMENT
                             ) {
	for (pdb_residue &my_pdb_residue : pdb_residues) {
		my_pdb_residue.set_chain_label( arg_chain_label );
	}
}

/// \brief TODOCUMENT
residue_name_vec pdb::do_get_residue_names_of_first_chain__backbone_unchecked() const {
	return get_backbone_complete_residue_names_of_first_chain( false );
}

/// \brief TODOCUMENT
coord pdb::do_get_residue_ca_coord_of_index__backbone_unchecked(const size_t &arg_index ///< TODOCUMENT
                                                                ) const {
	return get_carbon_alpha_coord_of_residue( get_residue_cref_of_index__backbone_unchecked( arg_index ) );
}

/// \brief TODOCUMENT
size_t pdb::do_get_num_atoms() const {
	auto num_atoms = 0_z;
	for (const pdb_residue &residue : pdb_residues) {
		num_atoms += residue.get_num_atoms();
	}
	return num_atoms;
}

/// \brief TODOCUMENT
void pdb::do_rotate(const rotation &arg_rotation ///< TODOCUMENT
                    ) {
	for (pdb_residue &my_pdb_residue : pdb_residues) {
		my_pdb_residue.rotate( arg_rotation );
	}
}

/// \brief TODOCUMENT
void pdb::do_add(const coord &arg_coord ///< TODOCUMENT
                 ) {
	for (pdb_residue &my_pdb_residue : pdb_residues) {
		my_pdb_residue += arg_coord;
	}
}

/// \brief TODOCUMENT
void pdb::do_subtract(const coord &arg_coord ///< TODOCUMENT
                      ) {
	for (pdb_residue &my_pdb_residue : pdb_residues) {
		my_pdb_residue -= arg_coord;
	}
}

/// \brief TODOCUMENT
size_t pdb::get_num_backbone_complete_residues() const {
	return numeric_cast<size_t>( count_if(
		pdb_residues,
		[] (const pdb_residue &x) { return is_backbone_complete( x ); }
	) );
}

/// \brief TODOCUMENT
size_t pdb::get_index_of_backbone_complete_index(const size_t &arg_backbone_complete_index ///< TODOCUMENT
                                                 ) const {
	if ( arg_backbone_complete_index >= pdb_residues.size() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to get_residue_ca_coord_of_backbone_complete_index() for index >= number of residues"));
	}
	size_t backbone_complete_ctr = 0;
	for (size_t index = 0; index < pdb_residues.size(); ++index) {
		const pdb_residue &the_pdb_residue = pdb_residues[ index ];
		if ( is_backbone_complete( the_pdb_residue ) ) {
			if ( backbone_complete_ctr == arg_backbone_complete_index ) {
				return index;
			}
			++backbone_complete_ctr;
		}
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot find enough backbone_complete residues to reach backbone_complete_index " + lexical_cast<string>( arg_backbone_complete_index ) ));
	return 0; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
}

/// \brief TODOCUMENT
const pdb_residue & pdb::get_residue_cref_of_backbone_complete_index(const size_t &arg_backbone_complete_index  ///< TODOCUMENT
                                                                     ) const {
	const size_t index = get_index_of_backbone_complete_index( arg_backbone_complete_index );
	return get_residue_cref_of_index__backbone_unchecked( index );
}

/// \brief TODOCUMENT
void pdb::set_residues(const pdb_residue_vec &arg_pdb_residues ///< TODOCUMENT
                       ) {
	pdb_residues = arg_pdb_residues;
}

/// \brief TODOCUMENT
void pdb::set_residues(pdb_residue_vec &&arg_pdb_residues ///< TODOCUMENT
                       ) {
	pdb_residues = std::move( arg_pdb_residues );
}

/// \brief TODOCUMENT
coord pdb::get_residue_ca_coord_of_backbone_complete_index(const size_t &arg_backbone_complete_index ///< TODOCUMENT
                                                           ) const {
	const size_t index = get_index_of_backbone_complete_index( arg_backbone_complete_index );
	return get_residue_ca_coord_of_index__backbone_unchecked( index );
}

/// \brief TODOCUMENT
residue_name_vec pdb::get_backbone_complete_residue_names_of_first_chain(const bool &arg_complete_backbone_only  ///< TODOCUMENT
                                                                         ) const {
	residue_name_vec residue_names;
	residue_names.reserve( pdb_residues.size() );
	if ( ! pdb_residues.empty() ) {
		const chain_label first_chain_label = pdb_residues.front().get_chain_label();
		for (const pdb_residue &my_pdb_residue : pdb_residues) {
			if ( first_chain_label == my_pdb_residue.get_chain_label() ) {
				if ( is_backbone_complete( my_pdb_residue ) || ! arg_complete_backbone_only ) {
					residue_names.push_back( my_pdb_residue.get_residue_name() );
				}
			}
		}
	}
	return residue_names;
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
///
/// \relatesalso domain
pdb cath::file::read_domain_from_pdb_file(const path   &arg_pdb_filename, ///< TODOCUMENT
                                          const domain &/*arg_domain*/    ///< TODOCUMENT
                                          ) {
	BOOST_THROW_EXCEPTION(not_implemented_exception("This is not implemented this here"));
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
	pdb_residue_vec residues;

	string          line_string;

	chain_label     prev_chain( ' ' );
	string          prev_amino_acid_3_char_code;
	residue_name    prev_res_name;
	pdb_atom_vec    prev_atoms;
	bool            prev_warned_conflict = false;

	// Loop over the lines of the file
	//
	// This code is made a bit more complicated because the aim is to
	// add all of the atoms within a residue at the same time but it isn't
	// clear that the residue has finished until the first line of the next residue
	// (or the end of the file)
	while ( getline( input_stream, line_string ) ) {
		// If this line is an ATOM or HETATM record
		if ( is_pdb_record_of_type( line_string, pdb_record::ATOM ) || is_pdb_record_of_type( line_string, pdb_record::HETATM ) ) {
			const auto parse_status_str_and_aa = pdb_record_parse_problem( line_string );
			const auto &parse_status = get<0>( parse_status_str_and_aa );
			const auto &parse_string = get<1>( parse_status_str_and_aa );
			const auto &parse_aa     = get<2>( parse_status_str_and_aa );
			if ( parse_status == pdb_atom_parse_status::ABORT ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception(
					"ATOM record is malformed : " + parse_string
					+ "\nRecord was \"" + line_string.substr(0, pdb_base::MAX_NUM_PDB_COLS)
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
			const chain_resname_atom_tuple  new_entry              = parse_pdb_atom_record( line_string, parse_aa );
			const chain_label              &chain                  = get<0>( new_entry );
			const residue_name             &res_name               = get<1>( new_entry );
			const pdb_atom                 &atom                   = get<2>( new_entry );
			const string                    amino_acid_3_char_code = get_amino_acid_code( atom );

			// If there are previously seen atoms that don't match this chain/res_name,
			// then add those atoms' residue and reset prev_atoms
			if ( ! prev_atoms.empty() && ( chain != prev_chain || res_name != prev_res_name ) ) {
				residues.emplace_back(
					prev_chain,
					prev_res_name,
					std::move( prev_atoms )
				);
				prev_atoms = pdb_atom_vec{};
				prev_warned_conflict = false;
			}

			// Some PDBs (eg 4tsw) may have erroneous consecutive duplicate residues.
			// Though that's a bit rubbish, it shouldn't break the whole comparison
			// so if that's detected, just warn and move on (without appending to new_residues).
			if ( ! prev_atoms.empty() && amino_acid_3_char_code != prev_amino_acid_3_char_code ) {
				if ( ! prev_warned_conflict ) {
					BOOST_LOG_TRIVIAL( warning ) << "Whilst parsing PDB file, found conflicting consecutive entries for residue \""
					                             << res_name
					                             << "\" on chain '"
					                             << chain
					                             << "' (with amino acids \""
					                             << prev_amino_acid_3_char_code
					                             << "\" and then \""
					                             << amino_acid_3_char_code
					                             << "\") - ignoring latter entry (and any further entries)";
					prev_warned_conflict = true;
				}
			}
			// Otherwise update the records of previously seen atoms
			else {
				prev_chain                  = chain;
				prev_amino_acid_3_char_code = amino_acid_3_char_code;
				prev_res_name               = res_name;
				prev_atoms.push_back( atom );
			}
		}
	};

	// Add any last remaining atoms
	if ( ! prev_atoms.empty() ) {
		residues.emplace_back(
			prev_chain,
			prev_res_name,
			std::move( prev_atoms )
		);
	}

	arg_pdb.set_residues( std::move( residues ) );
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

/// \brief TODOCUMENT
///
/// \relates pdb
ostream & cath::file::write_pdb_file(ostream   &arg_os,     ///< TODOCUMENT
                                     const pdb &arg_pdb ///< TODOCUMENT
                                     ) {
	const size_t num_residues = arg_pdb.get_num_residues();
//	size_t atom_ctr = 1;
	for (size_t residue_ctr = 0; residue_ctr < num_residues; ++residue_ctr) {
		const pdb_residue &residue = arg_pdb.get_residue_cref_of_index__backbone_unchecked( residue_ctr );
		write_pdb_file_entry( arg_os, residue );
//		write_pdb_file_entry( arg_os, residue, atom_ctr );
//		atom_ctr += residue.get_num_atoms();
	}
	arg_os << pdb_base::PDB_RECORD_STRING_TER << "\n";
	return arg_os;
}

/// \brief TODOCUMENT
///
/// \relates pdb
amino_acid_vec cath::file::get_amino_acid_list(const pdb &arg_pdb ///< TODOCUMENT
                                               ) {
	amino_acid_vec amino_acids;
	amino_acids.reserve( arg_pdb.get_num_residues() );
	for (const pdb_residue &the_pdb_residue : arg_pdb) {
		amino_acids.push_back( get_amino_acid( the_pdb_residue ) );
	}
	return amino_acids;
}

/// \brief TODOCUMENT
///
/// \relates pdb
ostream & cath::file::operator<<(ostream   &arg_os,         ///< TODOCUMENT
                                 const pdb &arg_pdb_residue ///< TODOCUMENT
                                 ) {
	arg_os << "PDB[";

	const size_t num_residues = arg_pdb_residue.get_num_residues();
	for (size_t residue_ctr = 0; residue_ctr < num_residues; ++residue_ctr) {
		arg_os << arg_pdb_residue.get_residue_cref_of_index__backbone_unchecked(residue_ctr);
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
//	return arg_pdb.get_residue_cref_of_index__backbone_unchecked(arg_index - 1);
//}

/// \brief Get the list of phi/psi angle pairs for each residue, in radians within (0, 2 * pi]
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
doub_angle_doub_angle_pair_vec cath::file::get_phi_and_psi_angles(const pdb &arg_pdb ///< TODOCUMENT
                                                                  ) {
	// The gap between consecutive residues' carbon-alpha atoms, above which the residues are not treated as neighbours
	//
	// To replicate DSSP's behaviour, this must be at least >= 5.01602 so that residues 161 and 162 on chain A of 139l are treated as connected
	// To replicate DSSP's behaviour, this must be at least >= 6.60324 so that residues  63 and  64 on chain A of 1c8c are treated as connected
	const double INTER_CA_DIST_FOR_NEIGHBOURS = 6.625; //< 6 + 5/8
	const auto   DEFAULT_PHI_PSI              = residue::DEFAULT_PHI_PSI();
	const auto   DEFAULT_PHI_PSI_PAIR         = make_pair( DEFAULT_PHI_PSI, DEFAULT_PHI_PSI );

	const size_t num_residues = arg_pdb.get_num_residues();
	doub_angle_doub_angle_pair_vec phi_and_psi_angles( num_residues, DEFAULT_PHI_PSI_PAIR );
	for (size_t residue_ctr = 0; residue_ctr + 1 < num_residues; ++residue_ctr) {

		// Grab this residue and the next one and calculate this psi and the next phi
		const pdb_residue    &this_pdb_residue      = arg_pdb.get_residue_cref_of_index__backbone_unchecked( residue_ctr     );
		const pdb_residue    &next_pdb_residue      = arg_pdb.get_residue_cref_of_index__backbone_unchecked( residue_ctr + 1 );
		const double          inter_ca_dist         = distance_between_points(
			get_carbon_alpha_coord_of_residue( this_pdb_residue ),
			get_carbon_alpha_coord_of_residue( next_pdb_residue )
		);

		// If these consecutive residues are a sensible distance apart...
		if ( inter_ca_dist < INTER_CA_DIST_FOR_NEIGHBOURS ) {
			// Calculate the two angles of these two residue and store them
			const auto this_psi_and_next_phi = get_psi_of_this_and_phi_of_next( this_pdb_residue, next_pdb_residue );
			phi_and_psi_angles[ residue_ctr     ].second = this_psi_and_next_phi.first;
			phi_and_psi_angles[ residue_ctr + 1 ].first  = this_psi_and_next_phi.second;
		}
	}

	// Return the calculated angles
	return phi_and_psi_angles;
}

/// \brief TODOCUMENT
///
/// \relates pdb
pdb cath::file::backbone_complete_subset_of_pdb(const pdb             &arg_pdb,             ///< TODOCUMENT
                                                const ostream_ref_opt &arg_ostream_ref_opt, ///< An optional reference to an ostream to which any logging should be sent
                                                const bool            &arg_skip_like_dssp   ///< TODOCUMENT
                                                ) {
	// Grab the number of residues
	const size_t num_residues = arg_pdb.get_num_residues();

	residue_name_vec backbone_skipped_residues;

	// Prepare a vector of the new residues
	pdb_residue_vec new_pdb_residues;
	new_pdb_residues.reserve(num_residues);

	// Loop over the residues in the input pdb
	for (size_t residue_ctr = 0; residue_ctr < num_residues; ++residue_ctr) {
		const pdb_residue &the_residue = arg_pdb.get_residue_cref_of_index__backbone_unchecked( residue_ctr );

		// If the residue is backbone_complete,then add it to new_pdb_residues
		const bool     skip_like_dssp_and_ok =   arg_skip_like_dssp && ! dssp_will_skip_residue( the_residue );
		const bool not_skip_like_dssp_and_ok = ! arg_skip_like_dssp && is_backbone_complete( the_residue );
		if ( skip_like_dssp_and_ok || not_skip_like_dssp_and_ok ) {
			new_pdb_residues.push_back( the_residue );
		}
		// Else if this is a proper amino acid (not just a bunch of HETATMs), warn
		else if ( get_amino_acid_letter( the_residue ) && arg_ostream_ref_opt ) {
			backbone_skipped_residues.push_back( the_residue.get_residue_name() );
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
		                             		| transformed( [] (const residue_name &x) { return to_string( x ); } ),
		                             	", "
		                             )
		                             << " whilst extracting a protein structure from PDB file data because "
		                             << ( multiple_skippeds ? "they don't"s : "it doesn't"s )
		                             << " have all of N, CA and C atoms";
	}

	// Return a new pdb containing these residues
	pdb new_pdb;
	new_pdb.set_residues( new_pdb_residues );
	return new_pdb;
}

/// \brief TODOCUMENT
///
/// \relates pdb
///
/// \relates protein
protein cath::file::build_protein_of_pdb(const pdb             &arg_pdb,    ///< TODOCUMENT
                                         const ostream_ref_opt &arg_ostream ///< An optional reference to an ostream to which any logging should be sent
                                         ) {
	constexpr size_t DEFAULT_ACCESSIBILITY = 0;

	const pdb    backbone_complete_pdb_subset = backbone_complete_subset_of_pdb( arg_pdb, arg_ostream );
	const size_t num_residues                 = backbone_complete_pdb_subset.get_num_residues();
	const auto   phi_and_psi_angles           = get_phi_and_psi_angles( backbone_complete_pdb_subset );

	return build_protein( transform_build<residue_vec>(
		irange( 0_z, num_residues ),
		[&] (const size_t &x) {
			const pdb_residue &the_residue = backbone_complete_pdb_subset.get_residue_cref_of_index__backbone_unchecked( x );
			const auto        &phi         = phi_and_psi_angles[ x ].first;
			const auto        &psi         = phi_and_psi_angles[ x ].second;
			return build_residue_of_pdb_residue(
				the_residue,
				phi,
				psi,
				DEFAULT_ACCESSIBILITY
			);
		}
	) );
}

/// \brief TODOCUMENT
///
/// \relates pdb
///
/// \relates protein
protein cath::file::build_protein_of_pdb_and_name(const pdb             &arg_pdb,    ///< TODOCUMENT
                                                  const string          &arg_name,   ///< TODOCUMENT
                                                  const ostream_ref_opt &arg_ostream ///< An optional reference to an ostream to which any logging should be sent
                                                  ) {
	protein new_protein = build_protein_of_pdb( arg_pdb, arg_ostream );
	new_protein.set_title( arg_name );
	return new_protein;
}

/// \brief Generate a list of protein residue indices (corresponding to those returned by build_protein_of_pdb())
///        that DSSP might be expected to skip
///
/// \relates pdb
size_set cath::file::get_protein_res_indices_that_dssp_might_skip(const pdb             &arg_pdb,    ///< The PDB to query
                                                                  const ostream_ref_opt &arg_ostream ///< An ostream to which any status messages might be sent
                                                                  ) {
	const pdb    backbone_complete_pdb_subset = backbone_complete_subset_of_pdb( arg_pdb, arg_ostream );
	const size_t num_residues                 = backbone_complete_pdb_subset.get_num_residues();

	// Return the indices corresponding to residues with any atoms with non-standard alt_locn values
	return copy_build<size_set>(
		irange( 0_z, num_residues )
			| filtered( [&] (const size_t &x) {
				return dssp_might_skip_residue(
					backbone_complete_pdb_subset.get_residue_cref_of_index__backbone_unchecked( x )
				);
			} )
	);
}
