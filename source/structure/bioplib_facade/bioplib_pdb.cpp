/// \file
/// \brief The bioplib_pdb class definitions

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "bioplib_pdb.h"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "exception/invalid_argument_exception.h"
#include "exception/out_of_range_exception.h"
#include "exception/runtime_error_exception.h"
#include "structure/bioplib_facade/bioplib_interface.h"
#include "structure/geometry/coord.h"
#include "structure/geometry/rotation.h"
#include "structure/residue_name.h"

#include <iostream> // **** TEMPORARY ****

using namespace boost::algorithm;
using namespace boost::filesystem;
using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace std;

using boost::algorithm::starts_with;
using boost::algorithm::trim_copy;
using boost::lexical_cast;
using boost::numeric_cast;

/// \brief An private member for clearing up the internals of the PDB
///
/// This is called by the dtor if the cleanup hasn't already been done
void bioplib_pdb::free_pdb_internals() {
	check_ptr();

	// Index the PDB structure
	PDB **indx;
	int local_natoms;
	indx = IndexPDB(get_ptr(), &local_natoms);
	assert(indx != nullptr);

	// If in debug mode, check that the number of residues is as previously recorded
	const size_t previous_natoms = get_natoms();
	assert(previous_natoms == numeric_cast<size_t>(local_natoms));

	// Loop over the array of pointers returned by IndexPDB, freeing them one by one
	for (size_t ptr_ctr = 0; ptr_ctr < previous_natoms; ++ptr_ctr) {
		free(indx[ptr_ctr]);
	}

	// Free the array of pointers itself
	free(indx);

	// Set the PDB pointer to nullptr
	pdb_ptr = nullptr;
}

/// \brief TODOCUMENT
void bioplib_pdb::check_ptr() const {
	if (pdb_ptr == nullptr) {
		BOOST_THROW_EXCEPTION(out_of_range_exception("Attempt to use an uninitialised bioplib_pdb object"));
	}
}

/// \brief TODOCUMENT
PDB * bioplib_pdb::get_ptr() {
	check_ptr();
	return pdb_ptr;
}

/// \brief TODOCUMENT
const PDB * bioplib_pdb::get_ptr() const {
	check_ptr();
	return pdb_ptr;
}

/// \brief TODOCUMENT
size_t bioplib_pdb::get_natoms() const {
	check_ptr();
	return numeric_cast<size_t>(natoms);
}

/// This is the canonical copy-assignment operator using the copy-and-swap idiom.
/// There may be a nicer way to do it in C++11.
//bioplib_pdb & bioplib_pdb::operator=(const bioplib_pdb &arg_bioplib_pdb
//                                     ) {
//	bioplib_pdb temp(arg_bioplib_pdb); // do all the work off to the side
//
//	// "Commit" the work using nonthrowing operations only
//	swap(temp, *this);
//	return *this;
//}

/// \brief TODOCUMENT
void bioplib_pdb::do_read_file(const path &arg_filename
                               ) {
	if (pdb_ptr != nullptr) {
		free_pdb_internals();
	}

	FILE * in_file = fopen(arg_filename.string().c_str(), "r");
	if (in_file == nullptr) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to open PDB file \"" + arg_filename.string() + "\" for reading"));
	}
//	pdb_ptr = ReadPDB( in_file, &natoms);
	pdb_ptr = doReadPDB( in_file, &natoms, TRUE, 1, 1 );
	fclose(in_file);
}

/// \brief TODOCUMENT
void bioplib_pdb::do_append_to_file(const path &arg_filename
                                    ) const {
	check_ptr();

	FILE * out_file = fopen(arg_filename.string().c_str(), "a");
	if (out_file == nullptr) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to open PDB/superposition file \"" + arg_filename.string() + "\" for reading"));
	}
	WritePDB(out_file, pdb_ptr);
	fclose(out_file);
}

/// \brief Label structure with arg_chain_label
void bioplib_pdb::do_set_chain_label(const chain_label &arg_chain_label ///< TODOCUMENT
                                     ) {
	const string chain_label_str = arg_chain_label.to_string();
	if ( chain_label_str.length() != 1 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to set multi-character chain label from bioplib_pdb"));
	}
	check_ptr();
	const size_t my_natoms = get_natoms();
	PDB *current = get_ptr();
	for (size_t atom_ctr = 0; atom_ctr < my_natoms; ++atom_ctr) {
		snprintf(current->chain, 2, "%c", chain_label_str.at( 0 ) );
		current = current->next;
	}
}

/// \brief Get a list of the uniq-ed residue names in the original order
///
/// Any whitespace is trimmed off the insert code.
///
/// This is tested in get_residue_names_test.cpp
residue_name_vec bioplib_pdb::do_get_residue_names_of_first_chain__backbone_unchecked() const {
	check_ptr();
	residue_name_vec residue_names;

	// Loop over the atom records
	const size_t my_natoms = get_natoms();
	const PDB *current = get_ptr();
	for (size_t atom_ctr = 0; atom_ctr < my_natoms; ++atom_ctr) {
		// Construct a residue name from a lexical_cast of the resnum and a trim of a string of the insert
		const string       residue_number   = lexical_cast<string>( current->resnum );
		const string       residue_insert   = trim_copy( string( current->insert ) );
		const string       residue_name_str = residue_number + residue_insert;
		const residue_name res_name         = make_residue_name( residue_name_str  );

		// If this is the first atom or of the residue_name is different from the most recent one,
		// add it to the list
		if ( residue_names.empty() || residue_names.back() != res_name ) {
			residue_names.push_back( res_name );
		}

		// Step forward
		current = current->next;
	}

	// Return the results
	return residue_names;
}

/// \brief Grab the carbon alpha coordinates of the residue with the specified index
///
/// Indexing starts from 0.
///
/// This method of stepping through a linked list is a pretty inefficient way of finding a
/// particular entry. I've made it reasonably efficient, which is good enough for what I need now.
///
/// It'd be more efficient if it weren't necessary to step through a linked list.

coord bioplib_pdb::do_get_residue_ca_coord_of_index__backbone_unchecked(const size_t &arg_index ///< The index of the residue to query (starts from 0)
                                                                        ) const {
	check_ptr();

	// Initialise the number of atoms and the pointer to the start of the linked list
	const size_t my_natoms = get_natoms();
	const PDB *current = get_ptr();
	if (my_natoms <= 0) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to search for CA atom in bioplib_pdb with no entries"));
	}

	// Initialise the records of the first number and insert
//	const size_t NUM_INSERT_CHARS(sizeof(PDB::insert) / sizeof(PDB::insert[0]));
	// The previous line is quite sensible but doesn't compile on orengobuild64's g++ 4.1.2, so instead we do this:
	const size_t NUM_INSERT_CHARS(sizeof(((PDB*)nullptr)->insert) / sizeof(sizeof(((PDB*)nullptr)->insert[0])));

	char_vec insert_code(NUM_INSERT_CHARS+2, 0);
	strncpy(&insert_code.front(), current->insert, NUM_INSERT_CHARS+1);
	int resnum(current->resnum);
	size_t num_previous_residues(0);

	// Loop over the atom records
	for (size_t atom_ctr = 0; atom_ctr < my_natoms; ++atom_ctr) {
		// If this atom's number or insert is different from the previous atom's then increment num_previous_residues
		// and update the records of the most recent number and insert
		if (atom_ctr > 0 && (resnum != current->resnum || strcmp(&insert_code.front(), current->insert) != 0)) {
			++num_previous_residues;
			resnum = current->resnum;
			strncpy(&insert_code.front(), current->insert, NUM_INSERT_CHARS+1);
		}

		// If this line represents the carbon alpha of the residue of interest, return the coordinates
		if (arg_index == num_previous_residues && starts_with(current->atnam, "CA")) {
			return coord(current->x, current->y, current->z);
		}

		// Step forward
		current = current->next;
	}

	BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to find CA atom for residue with specified index in bioplib_pdb::get_residue_ca_coord_of_index__backbone_unchecked()"));
}

/// \brief Get the number of atoms (implemented simply by returning natoms)
size_t bioplib_pdb::do_get_num_atoms() const {
	check_ptr();
	return numeric_cast<size_t>( natoms );
}

/// \brief Rotate the whole structure according to the specified rotation
void bioplib_pdb::do_rotate(const rotation &arg_rotation ///< The specification of the rotation to perform
                            ) {
	check_ptr();
	double raw_rotation_matrix[3][3];
	for (size_t ctr1 = 0; ctr1 < coord::NUM_DIMS; ++ctr1) {
		for (size_t ctr2 = 0; ctr2 < coord::NUM_DIMS; ++ctr2) {
			// NOTE: Switching indices here because bioplib appears to index in the opposite way
			raw_rotation_matrix[ ctr1 ][ ctr2 ] = arg_rotation.get_value( ctr2, ctr1 );
		}
	}
	ApplyMatrixPDB(get_ptr(), raw_rotation_matrix);
}

/// \brief Translate the whole structure by adding a specified coord to all the structure's coordinates
///
/// This class uses boost::additive<> to automatically add other operators implemented in terms of this
void bioplib_pdb::do_add(const coord &arg_coord ///< The coord (vector) by which this structure should be translated
                         ) {
	check_ptr();
	TranslatePDB(get_ptr(), COOR_of_coord(arg_coord));
}

/// \brief Translate the whole structure by subtracting a specified coord from all the structure's coordinates
///
/// This class uses boost::additive<> to automatically add other operators implemented in terms of this
void bioplib_pdb::do_subtract(const coord &arg_coord ///< The coord (vector) by which this structure should be translated
                              ) {
	const coord neg_coord = -arg_coord;
	(*this) += neg_coord;
}

/// \brief TODOCUMENT
bioplib_pdb::bioplib_pdb(const bioplib_pdb &arg_bioplib_pdb
                         ) : pdb_base(),
                             pdb_ptr( DupePDB( arg_bioplib_pdb.pdb_ptr ) ),
                             natoms ( arg_bioplib_pdb.natoms             ) {
}

/// \brief TODOCUMENT
bioplib_pdb & bioplib_pdb::operator=(const bioplib_pdb &arg_bioplib_pdb
                                     ) {
	pdb_ptr = DupePDB( arg_bioplib_pdb.pdb_ptr );
	natoms  = arg_bioplib_pdb.natoms;
	return *this;
}

/// \brief TODOCUMENT
bioplib_pdb::~bioplib_pdb() noexcept {
	try {
		if (pdb_ptr != nullptr) {
			free_pdb_internals();
		}
	}
	catch (...) {
	}
}
