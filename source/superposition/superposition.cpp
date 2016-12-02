/// \file
/// \brief The superposition class definitions

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

#include "superposition.hpp"

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>

#include "common/algorithm/constexpr_is_uniq.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/string_algorithm/split_build.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "exception/not_implemented_exception.hpp"
#include "exception/runtime_error_exception.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "structure/bioplib_facade/bioplib_interface.hpp"
#include "structure/bioplib_facade/bioplib_pdb.hpp"
#include "structure/geometry/coord_list.hpp"
#include "structure/geometry/rotation.hpp"

#include <deque>
#include <fstream>
#include <numeric>
#include <set>

using namespace cath;
using namespace cath::common;
using namespace cath::geom;
using namespace cath::sup;
using namespace std;

using boost::algorithm::any_of;
using boost::algorithm::is_any_of;
using boost::algorithm::token_compress_on;
using boost::lexical_cast;
using std::make_tuple;

//const double superposition::INVALID_RMSD(-1.0);
constexpr size_t          superposition::NUM_ENTRIES_IN_PAIRWISE_SUPERPOSITION;
constexpr size_t          superposition::INDEX_OF_FIRST_IN_PAIRWISE_SUPERPOSITION;
constexpr size_t          superposition::INDEX_OF_SECOND_IN_PAIRWISE_SUPERPOSITION;
constexpr array<char, 52> superposition::SUPERPOSITION_CHAIN_CHARS;

/// \brief The list of standard chain labels to use if relabelling the chains of the structures in a superposition
///
/// \todo Can chain_label be made constexpr and this be converted to an array<chain_label, 52>?
const chain_label_vec     superposition::SUPERPOSITION_CHAIN_LABELS = transform_build<chain_label_vec>(
	SUPERPOSITION_CHAIN_CHARS,
	[] (const char &x) { return chain_label{ x }; }
);
static_assert( constexpr_is_uniq( superposition::SUPERPOSITION_CHAIN_CHARS ), "The list of superposition chain characters should not have any repeats" ) ;

/// \brief Find the translation and rotation that will make the second coord_list best fit the first in the first's current position
///
/// Note: this previously said "[...]make the first coord_list best fit the second in the second's current position[...]"
///
coord_rot_pair superposition::fit_second_to_first(const coord_list &arg_coords_a, ///< TODOCUMENT
                                                  const coord_list &arg_coords_b  ///< TODOCUMENT
                                                  ) {
	// Find the centres of gravity of each and generate versions that are each centred
	const coord      centre_of_gravity_a = centre_of_gravity( arg_coords_a );
	const coord      centre_of_gravity_b = centre_of_gravity( arg_coords_b );
	const coord_list centred_coords_a    = arg_coords_a - centre_of_gravity_a;
	const coord_list centred_coords_b    = arg_coords_b - centre_of_gravity_b;

	// Calculate rotation matrix
	//
	// Occasionally, if the original pair of coordinates are exactly 180 degrees from optimum, one
	// call to bioblip_fit can be insufficient to get a good answer.
	//
	// For example, this problem can occur with some subsets of the residues in chains B and C of 1iph
	//
	// Presumably this is because the algorithm struggles to find a first step of improvement when the
	// correct rotation is exactly 180 degrees, and it takes many iterations to build up a tiny perturbation.
	//
	// So here, a second call is made and the results combined with the first call to get a more robust answer.
	//
	// To provide extra protection against it failing to perturb the original at all,
	// the first rotation is calculated against the centred_coords_b twisted x->y->z
	// and then the second rotation is calculated against the centred_coords_b twisted
	// by the first rotation
	const coord_list twisted_centred_coords_b = rotate_copy( rotation::ROTATE_X_TO_Y_TO_Z_TO_X(), centred_coords_b         );
//	const coord_list twisted_centred_coords_b = centred_coords_b;
	const rotation   rotation_1               = bioplib_fit(centred_coords_a,                     twisted_centred_coords_b );
	const coord_list centred_coords_b_rot     = rotate_copy( rotation_1,                          centred_coords_b         );
	const rotation   rotation_2               = bioplib_fit( centred_coords_a,                    centred_coords_b_rot     );
	const rotation   match_rotation           = rotation_2 * rotation_1;

	// Calculate translation vector
	const coord      match_translation   = rotate_copy(transpose_copy(match_rotation), centre_of_gravity_a) - centre_of_gravity_b;

	// Return the results
	return make_pair(match_translation, match_rotation);
}

/// \brief Constructs a (potentially multiple) superposition from a bunch of pairs of coord_lists.
///
/// Each pair of coord_lists represents the common coords of two entries in the superposition
///
/// These pairs must form a spanning-tree on the list of entries.
///
/// It is useful to allow a superposition for a single entry so that cath-superpose doesn't
/// have to treat that as a special case.
superposition::superposition(const vector<indices_and_coord_lists_type> &arg_indices_and_coord_lists, ///< TODOCUMENT
                             const size_t                               &arg_base_index,       ///< Optionally specify which entry in the resulting superposition is the "base" (default: 0)
                             const coord                                &arg_base_translation, ///< Optionally specify the translation to apply to the base structure
                             const rotation                             &arg_base_rotation     ///< Optionally specify the rotation    to apply to the base structure
                             ) {
	// Sanity check the inputs: check that each entry in the input has indices within range
	const size_t num_entries = arg_indices_and_coord_lists.size() + 1; // A tree
	// Sanity check the inputs: check that there are enough entries in input
	if (num_entries < 1) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot generate a superposition for zero entries"));
	}
	for (const indices_and_coord_lists_type &indices_and_coord_lists_entry : arg_indices_and_coord_lists) {
		if ( get<0>( indices_and_coord_lists_entry ) >= num_entries || get<2>( indices_and_coord_lists_entry ) >= num_entries ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Superposition input data index exceeds or equals the number of entries"));
		}
		if ( get<0>( indices_and_coord_lists_entry ) == get<2>( indices_and_coord_lists_entry ) ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Superposition input data indices are equal"));
		}
	}
	// Sanity check the inputs: check that the arg_index_to_use_as_base is within range
	if (arg_base_index >= num_entries) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Superposition arg_index_to_use_as_base exceeds or equals the number of entries"));
	}

//	rmsds.assign(num_entries - 1, INVALID_RMSD);

	translations.assign( num_entries, coord::ORIGIN_COORD           );
	rotations.assign   ( num_entries, rotation::IDENTITY_ROTATION() );
	translations[ arg_base_index ] = arg_base_translation;
	rotations   [ arg_base_index ] = arg_base_rotation;

	size_set entries_loaded = { arg_base_index };
	bool_deq inputs_processed( arg_indices_and_coord_lists.size(), false );

	// While any of the inputs_processed are false
	while ( any_of( inputs_processed, logical_not<bool>() ) ) {
		bool made_progress = false;

		for (size_t input_ctr = 0; input_ctr < arg_indices_and_coord_lists.size(); ++input_ctr) {
			if ( ! inputs_processed[input_ctr]) {
				const indices_and_coord_lists_type &indices_and_coord_lists_entry = arg_indices_and_coord_lists[input_ctr];
				const size_t     &index_a      = get<0>( indices_and_coord_lists_entry );
				const coord_list &coord_list_a = get<1>( indices_and_coord_lists_entry );
				const size_t     &index_b      = get<2>( indices_and_coord_lists_entry );
				const coord_list &coord_list_b = get<3>( indices_and_coord_lists_entry );

				if (entries_loaded.count(index_a) && entries_loaded.count(index_b)) {
					BOOST_THROW_EXCEPTION(invalid_argument_exception(
						"Superposition cannot be created because there are loops (entry number "
						+ lexical_cast<string>(input_ctr)
						+ " joins the already-joined "
						+ lexical_cast<string>(index_a)
						+ " and "
						+ lexical_cast<string>(index_b)
						+ ")"
					));
				}
				else if (entries_loaded.count(index_a) || entries_loaded.count(index_b)) {
					const bool        b_is_to_be_matched_to_a = entries_loaded.count(index_a);
					const size_t     &source_index           = b_is_to_be_matched_to_a ? index_a      : index_b;
					const size_t     &target_index           = b_is_to_be_matched_to_a ? index_b      : index_a;
					const coord_list &source_coord_list      = b_is_to_be_matched_to_a ? coord_list_a : coord_list_b;
					const coord_list &target_coord_list      = b_is_to_be_matched_to_a ? coord_list_b : coord_list_a;
					const coord      &source_translation     = translations[ source_index ];
					const rotation   &source_rotation        = rotations[    source_index ];

					const coord_rot_pair new_trans_and_rotn = fit_second_to_first(
						rotate_copy(source_rotation, source_coord_list + source_translation),
						target_coord_list
					);

					translations[ target_index ] = new_trans_and_rotn.first;
					rotations   [ target_index ] = new_trans_and_rotn.second;
					inputs_processed[ input_ctr ]  = true;
					made_progress                  = true;
					entries_loaded.insert( target_index );
				}
			}
		}

		if ( ! made_progress ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Superposition cannot be created because input data does not contain links to connect all entries"));
		}
	}
}

/// \brief TODOCUMENT
superposition::superposition(const coord_vec    &arg_translations, ///< TODOCUMENT
                             const rotation_vec &rotations         ///< TODOCUMENT
                             ) : translations( arg_translations ),
                                 rotations   ( rotations        ) {
}

/// \brief TODOCUMENT
size_t superposition::get_num_entries() const {
	return translations.size();
}

/// \brief TODOCUMENT
const coord & superposition::get_translation_of_index(const size_t &arg_index ///< TODOCUMENT
                                                      ) const {
	return translations[arg_index];
}

/// \brief TODOCUMENT
const rotation & superposition::get_rotation_of_index(const size_t &arg_index ///< TODOCUMENT
                                                      ) const {
	return rotations[arg_index];
}

/// \brief TODOCUMENT
///
/// \relates superposition
void cath::sup::transform(const superposition &arg_superposition, ///< TODOCUMENT
                          const size_t        &arg_index,         ///< TODOCUMENT
                          coord               &arg_coord          ///< TODOCUMENT
                          ) {
	// Grab the translation and rotation
	const coord    &the_translation = arg_superposition.get_translation_of_index( arg_index );
	const rotation &the_rotation    = arg_superposition.get_rotation_of_index   ( arg_index );
	
	// Apply them to the specified coord
	arg_coord += the_translation;
	cath::geom::rotate( the_rotation, arg_coord );
}

/// \brief TODOCUMENT
///
/// \relates superposition
coord cath::sup::transform_copy(const superposition &arg_superposition, ///< TODOCUMENT
                                const size_t        &arg_index,         ///< TODOCUMENT
                                coord                arg_coord          ///< TODOCUMENT
                                ) {
	transform( arg_superposition, arg_index, arg_coord );
	return arg_coord;
}

/// \brief TODOCUMENT
///
/// \relates superposition
void cath::sup::transform(const superposition &arg_superposition, ///< TODOCUMENT
                          const size_t        &arg_index,         ///< TODOCUMENT
                          coord_list          &arg_coord_list     ///< TODOCUMENT
                          ) {
	// Grab the translation and rotation
	const coord    &the_translation = arg_superposition.get_translation_of_index( arg_index );
	const rotation &the_rotation    = arg_superposition.get_rotation_of_index   ( arg_index );
	
	// Apply them to the specified coord
	arg_coord_list += the_translation;
	cath::geom::rotate( the_rotation, arg_coord_list );
}

/// \brief TODOCUMENT
///
/// \relates superposition
coord_list cath::sup::transform_copy(const superposition &arg_superposition, ///< TODOCUMENT
                                     const size_t        &arg_index,         ///< TODOCUMENT
                                     coord_list           arg_coord_list     ///< TODOCUMENT
                                     ) {
	transform( arg_superposition, arg_index, arg_coord_list );
	return arg_coord_list;
}

/// \brief TODOCUMENT
///
/// \relates superposition
void cath::sup::transform(const superposition &arg_superposition, ///< TODOCUMENT
                          const size_t        &arg_index,         ///< TODOCUMENT
                          coord_list_vec      &arg_coord_lists    ///< TODOCUMENT
                          ) {
	for_each(
		arg_coord_lists,
		[&] (coord_list &x) { transform( arg_superposition, arg_index, x ); }
	);
}

/// \brief TODOCUMENT
///
/// \relates superposition
geom::coord_list_vec cath::sup::transform_copy(const superposition &arg_superposition, ///< TODOCUMENT
                                               const size_t        &arg_index,         ///< TODOCUMENT
                                               coord_list_vec       arg_coord_lists    ///< TODOCUMENT
                                               ) {
	transform( arg_superposition, arg_index, arg_coord_lists );
	return arg_coord_lists;
}

/// \brief TODOCUMENT
///
/// \relates superposition
double cath::sup::superposed_distance(const superposition &arg_superposition, ///< TODOCUMENT
                                      const size_t        &arg_index_a,       ///< TODOCUMENT
                                      coord                arg_coord_a,       ///< TODOCUMENT
                                      const size_t        &arg_index_b,       ///< TODOCUMENT
                                      coord                arg_coord_b        ///< TODOCUMENT
                                      ) {
	transform( arg_superposition, arg_index_a, arg_coord_a );
	transform( arg_superposition, arg_index_b, arg_coord_b );
	return distance_between_points( arg_coord_a, arg_coord_b );
}

/// \brief TODOCUMENT
///
/// \relates superposition
double cath::sup::calc_rmsd_between_superposed_entries(const superposition &arg_superposition, ///< TODOCUMENT
                                                       const size_t        &arg_index_a,       ///< TODOCUMENT
                                                       coord_list           arg_coord_list_a,  ///< TODOCUMENT
                                                       const size_t        &arg_index_b,       ///< TODOCUMENT
                                                       coord_list           arg_coord_list_b   ///< TODOCUMENT
                                                       ) {
	transform( arg_superposition, arg_index_a, arg_coord_list_a );
	transform( arg_superposition, arg_index_b, arg_coord_list_b );

	return calc_rmsd( arg_coord_list_a, arg_coord_list_b );
}

/// \brief Non-member equality operator for the superposition class
///
/// \relates superposition
bool cath::sup::operator==(const superposition &arg_sup_a, ///< TODOCUMENT
                           const superposition &arg_sup_b  ///< TODOCUMENT
                           ) {
	// Return false if the number of entries differ
	if (arg_sup_a.get_num_entries() != arg_sup_b.get_num_entries()) {
		return false;
	}

	// Otherwise, check each of the positions match
	const size_t num_entries = arg_sup_a.get_num_entries();
	for (size_t entry_ctr = 0; entry_ctr < num_entries; ++entry_ctr) {
		if (arg_sup_a.get_translation_of_index(entry_ctr) != arg_sup_b.get_translation_of_index(entry_ctr) ) {
			return false;
		}
		const rotation arg_rotation_a = arg_sup_a.get_rotation_of_index(entry_ctr);
		const rotation arg_rotation_b = arg_sup_b.get_rotation_of_index(entry_ctr);
		if ( arg_rotation_a != arg_rotation_b ) {
			return false;
		}
	}

	return true;
}

/// \brief Basic insertion operator to output a rough summary of an superposition to an ostream
///
/// \relates superposition
ostream & cath::sup::operator<<(ostream             &arg_ostream,      ///< The stream to which to output
                                const superposition &arg_superposition ///< The superposition to summarise
                                ) {
	const size_t num_entries = arg_superposition.get_num_entries();
	arg_ostream << "superposition[";
	arg_ostream << num_entries;
	arg_ostream << " entries:";
	for (size_t entry_ctr = 0; entry_ctr < num_entries; ++entry_ctr) {
		arg_ostream << (entry_ctr > 0 ? "; " : " ");
		arg_ostream << arg_superposition.get_translation_of_index(entry_ctr);
		arg_ostream << ", ";
		arg_ostream << arg_superposition.get_rotation_of_index(entry_ctr);
	}
	arg_ostream << "]";
	return arg_ostream;
}

/// \brief TODOCUMENT
///
/// \relates superposition
bool cath::sup::are_close(const superposition &arg_sup_1,  ///< TODOCUMENT
                          const superposition &arg_sup_2   ///< TODOCUMENT
                          ) {
	const size_t num_entries = arg_sup_1.get_num_entries();

	if (num_entries != arg_sup_2.get_num_entries()) {
		return false;
	}
	for (size_t entry_ctr = 0; entry_ctr < num_entries; ++entry_ctr) {
		const bool translations_are_close = (
			arg_sup_1.get_translation_of_index(entry_ctr) == arg_sup_2.get_translation_of_index(entry_ctr)
		);
		if (!translations_are_close) {
			return false;
		}
		const bool rotations_are_close = are_close(
			arg_sup_1.get_rotation_of_index(entry_ctr),      arg_sup_2.get_rotation_of_index(entry_ctr)
		);
		if (!rotations_are_close) {
			return false;
		}
	}
	return true;
}

/// \brief TODOCUMENT
///
/// \relates superposition
void cath::sup::write_superposition(ostream             &arg_os,           ///< TODOCUMENT
                                    const superposition &arg_superposition ///< TODOCUMENT
                                    ) {
	const streamsize old_precision = arg_os.precision();
	arg_os.precision(30);
	const size_t num_entries = arg_superposition.get_num_entries();
	for (size_t entry_ctr = 0; entry_ctr < num_entries; ++entry_ctr) {
		const coord    entry_translation = arg_superposition.get_translation_of_index( entry_ctr );
		const rotation entry_rotation    = arg_superposition.get_rotation_of_index(    entry_ctr );
		arg_os << entry_translation.get_x() << " " << entry_translation.get_y() << " " << entry_translation.get_z();
		arg_os << " "  << entry_rotation.get_value<0, 0>();
		arg_os << " "  << entry_rotation.get_value<0, 1>();
		arg_os << " "  << entry_rotation.get_value<0, 2>();
		arg_os << " "  << entry_rotation.get_value<1, 0>();
		arg_os << " "  << entry_rotation.get_value<1, 1>();
		arg_os << " "  << entry_rotation.get_value<1, 2>();
		arg_os << " "  << entry_rotation.get_value<2, 0>();
		arg_os << " "  << entry_rotation.get_value<2, 1>();
		arg_os << " "  << entry_rotation.get_value<2, 2>();
		arg_os << "\n";
	}
	arg_os.precision(old_precision);
}
/// \brief TODOCUMENT
superposition cath::sup::read_superposition(istream &arg_is ///< TODOCUMENT
		                                    ) {
	coord_vec    translations;
	rotation_vec rotations;
	const size_t CORRECT_NUM_PARTS(12);
	string line_string;
	while ( getline( arg_is, line_string ) ) {
		const str_vec line_parts = split_build<str_vec>( line_string, is_any_of( " " ), token_compress_on );
		if (line_parts.size() != CORRECT_NUM_PARTS) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to parse superposition because line does not contain " + lexical_cast<string>(CORRECT_NUM_PARTS) + " parts."));
		}

		translations.push_back(coord(
			stod( line_parts[  0 ] ),
			stod( line_parts[  1 ] ),
			stod( line_parts[  2 ] )
		));
		rotations.push_back(rotation(
			stod( line_parts[  3 ] ),
			stod( line_parts[  4 ] ),
			stod( line_parts[  5 ] ),
			stod( line_parts[  6 ] ),
			stod( line_parts[  7 ] ),
			stod( line_parts[  8 ] ),
			stod( line_parts[  9 ] ),
			stod( line_parts[ 10 ] ),
			stod( line_parts[ 11 ] )
		));
	}
	return superposition(translations, rotations);
}

/// \brief A convenience factory to create a pairwise superposition from two coord lists
///
/// \relates superposition
///
/// Rather than needing to use the more complicated superposition constructor,
/// this can be used to provide and interface to it.
superposition cath::sup::create_pairwise_superposition(const coord_list &arg_coord_list_1,     ///< TODOCUMENT
                                                       const coord_list &arg_coord_list_2,     ///< TODOCUMENT
                                                       const bool       &arg_first_as_base,    ///< TODOCUMENT
                                                       const coord      &arg_base_translation, ///< TODOCUMENT
                                                       const rotation   &arg_base_rotation     ///< TODOCUMENT
                                                       ) {
	const superposition new_superposition(
		vector<superposition::indices_and_coord_lists_type>( {
			std::make_tuple(
				superposition::INDEX_OF_FIRST_IN_PAIRWISE_SUPERPOSITION,
				arg_coord_list_1,
				superposition::INDEX_OF_SECOND_IN_PAIRWISE_SUPERPOSITION,
				arg_coord_list_2
			)
		} ),
		( arg_first_as_base ? superposition::INDEX_OF_FIRST_IN_PAIRWISE_SUPERPOSITION : superposition::INDEX_OF_SECOND_IN_PAIRWISE_SUPERPOSITION ),
		arg_base_translation,
		arg_base_rotation
	);
	return new_superposition;
}

/// \brief A convenience function to superpose one set of coordinates to another
///
/// \relates superposition
void cath::sup::superpose_second_coords_to_first(const coord_list &arg_coord_list_1, ///< The first set of coordinates (to which the others should be superposed)
                                                 coord_list       &arg_coord_list_2  ///< The second set of coordinates (to superpose to the others)
                                                 ) {
	const auto the_sup = create_pairwise_superposition( arg_coord_list_1, arg_coord_list_2 );
	transform( the_sup, superposition::INDEX_OF_SECOND_IN_PAIRWISE_SUPERPOSITION, arg_coord_list_2 );
}

/// \brief A convenience function to create a copy of one set of coordinates superposed to another
///
/// \relates superposition
coord_list cath::sup::superpose_copy_second_coords_to_first(const coord_list &arg_coord_list_1, ///< TODOCUMENT
                                                            coord_list        arg_coord_list_2  ///< TODOCUMENT
                                                            ) {
	superpose_second_coords_to_first( arg_coord_list_1, arg_coord_list_2 );
	return arg_coord_list_2;
}


/// \brief A convenience function to superpose one set of coordinates to another
///
/// \relates superposition
void cath::sup::superpose_second_coords_to_first(const coord_list_vec &arg_coord_list_1, ///< The first set of coordinates (to which the others should be superposed)
                                                 coord_list_vec       &arg_coord_list_2  ///< The second set of coordinates (to superpose to the others)
                                                 ) {
	const auto the_sup = create_pairwise_superposition(
		flatten_coord_lists( arg_coord_list_1 ),
		flatten_coord_lists( arg_coord_list_2 )
	);
	transform( the_sup, superposition::INDEX_OF_SECOND_IN_PAIRWISE_SUPERPOSITION, arg_coord_list_2 );
}

/// \brief A convenience function to create a copy of one set of coordinates superposed to another
///
/// \relates superposition
coord_list_vec cath::sup::superpose_copy_second_coords_to_first(const coord_list_vec &arg_coord_list_1, ///< TODOCUMENT
                                                                coord_list_vec        arg_coord_list_2  ///< TODOCUMENT
                                                                ) {
	superpose_second_coords_to_first( arg_coord_list_1, arg_coord_list_2 );
	return arg_coord_list_2;
}


/// \brief A convenience factory to superpose two coord lists and return the RMSD
///
/// \relates superposition
///
/// Rather than needing to use the more complicated
double cath::sup::calc_pairwise_superposition_rmsd(const coord_list &arg_coord_list_1, ///< TODOCUMENT
                                                   const coord_list &arg_coord_list_2  ///< TODOCUMENT
                                                   ) {
	const superposition new_sup = create_pairwise_superposition(
		arg_coord_list_1,
		arg_coord_list_2
	);
	return calc_rmsd_between_superposed_entries(
		new_sup,
		superposition::INDEX_OF_FIRST_IN_PAIRWISE_SUPERPOSITION,
		arg_coord_list_1,
		superposition::INDEX_OF_SECOND_IN_PAIRWISE_SUPERPOSITION,
		arg_coord_list_2
	);
}

/// \brief TODOCUMENT
///
/// \relates superposition
void cath::sup::check_superposition_is_pairwise(const superposition &arg_superposition ///< TODOCUMENT
                                                ) {
	if (arg_superposition.get_num_entries() != superposition::NUM_ENTRIES_IN_PAIRWISE_SUPERPOSITION) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Superposition is not pairwise"));
	}
}

