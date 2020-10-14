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
#include <boost/range/algorithm/for_each.hpp>
#include <boost/range/combine.hpp>

#include "cath/common/algorithm/constexpr_is_uniq.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/boost_addenda/string_algorithm/split_build.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/exception/not_implemented_exception.hpp"
#include "cath/common/exception/runtime_error_exception.hpp"
#include "cath/file/pdb/pdb.hpp"
#include "cath/file/pdb/pdb_atom.hpp"
#include "cath/file/pdb/pdb_residue.hpp"
#include "cath/structure/geometry/coord_list.hpp"
#include "cath/structure/geometry/rotation.hpp"
#include "cath/structure/geometry/superpose_fit.hpp"

#include <deque>
#include <fstream>
#include <numeric>
#include <set>

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::geom;
using namespace ::cath::sup;
using namespace ::std;

using ::boost::algorithm::any_of;
using ::boost::algorithm::is_any_of;
using ::boost::algorithm::token_compress_on;
using ::boost::lexical_cast;
using ::boost::range::combine;
using ::boost::range::for_each;
using ::std::make_tuple;

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
coord_rot_pair superposition::fit_second_to_first(const coord_list &prm_coords_a, ///< TODOCUMENT
                                                  const coord_list &prm_coords_b  ///< TODOCUMENT
                                                  ) {
	// Find the centres of gravity of each and generate versions that are each centred
	const coord      centre_of_gravity_a = centre_of_gravity( prm_coords_a );
	const coord      centre_of_gravity_b = centre_of_gravity( prm_coords_b );
	const coord_list centred_coords_a    = prm_coords_a - centre_of_gravity_a;
	const coord_list centred_coords_b    = prm_coords_b - centre_of_gravity_b;

	// Calculate rotation matrix
	//
	// This used to be complicated by the old fitting routine's difficulty with two sets of coordinates
	// that are exactly 180-degree rotations of each other. I think that because that routine was iterative,
	// it had some difficulties with situations in which it could find a direction to start the improvements.
	// Cf a ball-bearing that's perfectly balanced at the top of a smooth mound and doesn't roll in any direction.
	//
	// For example, this problem occurred with chains B and C of 1iph.
	//
	// To tackle this, the code previously calculated a first rotation against the centred_coords_b
	// twisted x->y->z and then a second rotation against the centred_coords_b twisted
	// by the first rotation.
	//
	// Now that the non-iterative Kabsch algorithm is being used, 1iphB/1iphC appear to superpose
	// correctly and that extra code has been removed.
	const rotation   match_rotation      = superpose_fit_2nd_to_1st( centred_coords_a, centred_coords_b );

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
superposition::superposition(const vector<indices_and_coord_lists_type> &prm_indices_and_coord_lists, ///< TODOCUMENT
                             const size_t                               &prm_base_index,       ///< Optionally specify which entry in the resulting superposition is the "base" (default: 0)
                             const coord                                &prm_base_translation, ///< Optionally specify the translation to apply to the base structure
                             const rotation                             &prm_base_rotation     ///< Optionally specify the rotation    to apply to the base structure
                             ) {
	// Sanity check the inputs: check that each entry in the input has indices within range
	const size_t num_entries = prm_indices_and_coord_lists.size() + 1; // A tree
	// Sanity check the inputs: check that there are enough entries in input
	if (num_entries < 1) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot generate a superposition for zero entries"));
	}
	for (const indices_and_coord_lists_type &indices_and_coord_lists_entry : prm_indices_and_coord_lists) {
		if ( get<0>( indices_and_coord_lists_entry ) >= num_entries || get<2>( indices_and_coord_lists_entry ) >= num_entries ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Superposition input data index exceeds or equals the number of entries"));
		}
		if ( get<0>( indices_and_coord_lists_entry ) == get<2>( indices_and_coord_lists_entry ) ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Superposition input data indices are equal"));
		}
	}
	// Sanity check the inputs: check that the prm_index_to_use_as_base is within range
	if (prm_base_index >= num_entries) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Superposition prm_index_to_use_as_base exceeds or equals the number of entries"));
	}

//	rmsds.assign(num_entries - 1, INVALID_RMSD);

	translations.assign( num_entries, ORIGIN_COORD                  );
	rotations.assign   ( num_entries, rotation::IDENTITY_ROTATION() );
	translations[ prm_base_index ] = prm_base_translation;
	rotations   [ prm_base_index ] = prm_base_rotation;

	size_set entries_loaded = { prm_base_index };
	bool_deq inputs_processed( prm_indices_and_coord_lists.size(), false );

	// While any of the inputs_processed are false
	while ( any_of( inputs_processed, logical_not<bool>() ) ) {
		bool made_progress = false;

		for (const size_t &input_ctr : indices( prm_indices_and_coord_lists.size() ) ) {
			if ( ! inputs_processed[input_ctr]) {
				const indices_and_coord_lists_type &indices_and_coord_lists_entry = prm_indices_and_coord_lists[input_ctr];
				const size_t     &index_a      = get<0>( indices_and_coord_lists_entry );
				const coord_list &coord_list_a = get<1>( indices_and_coord_lists_entry );
				const size_t     &index_b      = get<2>( indices_and_coord_lists_entry );
				const coord_list &coord_list_b = get<3>( indices_and_coord_lists_entry );

				if ( contains( entries_loaded, index_a ) && contains( entries_loaded, index_b ) ) {
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
				else if ( contains( entries_loaded, index_a ) || contains( entries_loaded, index_b ) ) {
					const bool        b_is_to_be_matched_to_a = contains( entries_loaded, index_a );
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
superposition::superposition(coord_vec    prm_translations, ///< TODOCUMENT
                             rotation_vec prm_rotations     ///< TODOCUMENT
                             ) : translations{ std::move( prm_translations ) },
                                 rotations   { std::move( prm_rotations    ) } {
}

/// \brief TODOCUMENT
size_t superposition::get_num_entries() const {
	return translations.size();
}

/// \brief TODOCUMENT
const coord & superposition::get_translation_of_index(const size_t &prm_index ///< TODOCUMENT
                                                      ) const {
	return translations[ prm_index ];
}

/// \brief TODOCUMENT
const rotation & superposition::get_rotation_of_index(const size_t &prm_index ///< TODOCUMENT
                                                      ) const {
	return rotations[ prm_index ];
}

/// \brief Modify the superposition so that it has the same effect as applying the original superposition and then the specified translation
superposition & superposition::post_translate(const coord &prm_translation ///< The post-superposition translation to add into the superposition
                                              ) {
	// \TODO Come C++17 and structured bindings, use here
	for (const boost::tuple<coord &, const rotation &> &x : combine( translations, rotations ) ) {
		auto       &the_trans = x.get<0>();
		const auto &the_rotn  = x.get<1>();
		the_trans += rotate_copy( transpose_copy( the_rotn ), prm_translation );
	}
	return *this;
}

/// \brief Modify the superposition so that it has the same effect as applying the original superposition and then the specified rotation
superposition & superposition::post_rotate(const rotation &prm_rotation ///< The post-superposition rotation to add into the superposition
                                           ) {
	for (rotation &the_rotn : rotations) {
		// \todo Consider adding a `rotation & left_multiply_by(const rotation &);` method to rotation and using it here
		the_rotn = prm_rotation * the_rotn;
	}
	return *this;
}

/// \brief Modify the specified superposition so that it has the same effect as
///        applying the original superposition and then the specified translation and then rotation
///
/// \relates superposition
void cath::sup::post_translate_and_rotate(superposition  &prm_superposition, ///< The superposition to modify
                                          const coord    &prm_translation,   ///< The translation to apply after the superposition
                                          const rotation &prm_rotation       ///< The rotation to apply last, after the translation
                                          ) {
	prm_superposition.post_translate( prm_translation );
	prm_superposition.post_rotate   ( prm_rotation    );
}

/// \brief Modify and return a copy of the specified superposition so that it has the same effect as
///        applying the original superposition and then the specified translation and then rotation
///
/// \relates superposition
superposition cath::sup::post_translate_and_rotate_copy(superposition   prm_superposition, ///< The superposition from which a copy should be taken, so it can be modified and returned
                                                        const coord    &prm_translation,   ///< The translation to apply after the superposition
                                                        const rotation &prm_rotation       ///< The rotation to apply last, after the translation
                                                        ) {
	post_translate_and_rotate( prm_superposition, prm_translation, prm_rotation );
	return prm_superposition;
}

/// \brief TODOCUMENT
///
/// \relates superposition
void cath::sup::transform(const superposition &prm_superposition, ///< TODOCUMENT
                          const size_t        &prm_index,         ///< TODOCUMENT
                          coord               &prm_coord          ///< TODOCUMENT
                          ) {
	// Grab the translation and rotation
	const coord    &the_translation = prm_superposition.get_translation_of_index( prm_index );
	const rotation &the_rotation    = prm_superposition.get_rotation_of_index   ( prm_index );
	
	// Apply them to the specified coord
	prm_coord += the_translation;
	cath::geom::rotate( the_rotation, prm_coord );
}

/// \brief TODOCUMENT
///
/// \relates superposition
coord cath::sup::transform_copy(const superposition &prm_superposition, ///< TODOCUMENT
                                const size_t        &prm_index,         ///< TODOCUMENT
                                coord                prm_coord          ///< TODOCUMENT
                                ) {
	transform( prm_superposition, prm_index, prm_coord );
	return prm_coord;
}

/// \brief TODOCUMENT
///
/// \relates superposition
void cath::sup::transform(const superposition &prm_superposition, ///< TODOCUMENT
                          const size_t        &prm_index,         ///< TODOCUMENT
                          coord_list          &prm_coord_list     ///< TODOCUMENT
                          ) {
	// Grab the translation and rotation
	const coord    &the_translation = prm_superposition.get_translation_of_index( prm_index );
	const rotation &the_rotation    = prm_superposition.get_rotation_of_index   ( prm_index );
	
	// Apply them to the specified coord
	prm_coord_list += the_translation;
	cath::geom::rotate( the_rotation, prm_coord_list );
}

/// \brief TODOCUMENT
///
/// \relates superposition
coord_list cath::sup::transform_copy(const superposition &prm_superposition, ///< TODOCUMENT
                                     const size_t        &prm_index,         ///< TODOCUMENT
                                     coord_list           prm_coord_list     ///< TODOCUMENT
                                     ) {
	transform( prm_superposition, prm_index, prm_coord_list );
	return prm_coord_list;
}

/// \brief TODOCUMENT
///
/// \relates superposition
void cath::sup::transform(const superposition &prm_superposition, ///< TODOCUMENT
                          const size_t        &prm_index,         ///< TODOCUMENT
                          coord_list_vec      &prm_coord_lists    ///< TODOCUMENT
                          ) {
	for_each(
		prm_coord_lists,
		[&] (coord_list &x) { transform( prm_superposition, prm_index, x ); }
	);
}

/// \brief TODOCUMENT
///
/// \relates superposition
geom::coord_list_vec cath::sup::transform_copy(const superposition &prm_superposition, ///< TODOCUMENT
                                               const size_t        &prm_index,         ///< TODOCUMENT
                                               coord_list_vec       prm_coord_lists    ///< TODOCUMENT
                                               ) {
	transform( prm_superposition, prm_index, prm_coord_lists );
	return prm_coord_lists;
}

/// \brief TODOCUMENT
///
/// \relates superposition
double cath::sup::superposed_distance(const superposition &prm_superposition, ///< TODOCUMENT
                                      const size_t        &prm_index_a,       ///< TODOCUMENT
                                      coord                prm_coord_a,       ///< TODOCUMENT
                                      const size_t        &prm_index_b,       ///< TODOCUMENT
                                      coord                prm_coord_b        ///< TODOCUMENT
                                      ) {
	transform( prm_superposition, prm_index_a, prm_coord_a );
	transform( prm_superposition, prm_index_b, prm_coord_b );
	return distance_between_points( prm_coord_a, prm_coord_b );
}

/// \brief TODOCUMENT
///
/// \relates superposition
double cath::sup::calc_rmsd_between_superposed_entries(const superposition &prm_superposition, ///< TODOCUMENT
                                                       const size_t        &prm_index_a,       ///< TODOCUMENT
                                                       coord_list           prm_coord_list_a,  ///< TODOCUMENT
                                                       const size_t        &prm_index_b,       ///< TODOCUMENT
                                                       coord_list           prm_coord_list_b   ///< TODOCUMENT
                                                       ) {
	transform( prm_superposition, prm_index_a, prm_coord_list_a );
	transform( prm_superposition, prm_index_b, prm_coord_list_b );

	return calc_rmsd( prm_coord_list_a, prm_coord_list_b );
}

/// \brief Non-member equality operator for the superposition class
///
/// \relates superposition
bool cath::sup::operator==(const superposition &prm_sup_a, ///< TODOCUMENT
                           const superposition &prm_sup_b  ///< TODOCUMENT
                           ) {
	// Return false if the number of entries differ
	if (prm_sup_a.get_num_entries() != prm_sup_b.get_num_entries()) {
		return false;
	}

	// Otherwise, check each of the positions match
	const size_t num_entries = prm_sup_a.get_num_entries();
	for (const size_t &entry_ctr : indices( num_entries ) ) {
		if (prm_sup_a.get_translation_of_index(entry_ctr) != prm_sup_b.get_translation_of_index(entry_ctr) ) {
			return false;
		}
		const rotation prm_rotation_a = prm_sup_a.get_rotation_of_index(entry_ctr);
		const rotation prm_rotation_b = prm_sup_b.get_rotation_of_index(entry_ctr);
		if ( prm_rotation_a != prm_rotation_b ) {
			return false;
		}
	}

	return true;
}

/// \brief Basic insertion operator to output a rough summary of an superposition to an ostream
///
/// \relates superposition
ostream & cath::sup::operator<<(ostream             &prm_ostream,      ///< The stream to which to output
                                const superposition &prm_superposition ///< The superposition to summarise
                                ) {
	const size_t num_entries = prm_superposition.get_num_entries();
	prm_ostream << "superposition[";
	prm_ostream << num_entries;
	prm_ostream << " entries:";
	for (const size_t &entry_ctr : indices( num_entries ) ) {
		prm_ostream << (entry_ctr > 0 ? "; " : " ");
		prm_ostream << prm_superposition.get_translation_of_index(entry_ctr);
		prm_ostream << ", ";
		prm_ostream << prm_superposition.get_rotation_of_index(entry_ctr);
	}
	prm_ostream << "]";
	return prm_ostream;
}

/// \brief TODOCUMENT
///
/// \relates superposition
bool cath::sup::are_close(const superposition &prm_sup_1,  ///< TODOCUMENT
                          const superposition &prm_sup_2   ///< TODOCUMENT
                          ) {
	const size_t num_entries = prm_sup_1.get_num_entries();

	if (num_entries != prm_sup_2.get_num_entries()) {
		return false;
	}
	for (const size_t &entry_ctr : indices( num_entries ) ) {
		const bool translations_are_close = (
			prm_sup_1.get_translation_of_index(entry_ctr) == prm_sup_2.get_translation_of_index(entry_ctr)
		);
		if (!translations_are_close) {
			return false;
		}
		const bool rotations_are_close = are_close(
			prm_sup_1.get_rotation_of_index(entry_ctr),      prm_sup_2.get_rotation_of_index(entry_ctr)
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
void cath::sup::write_superposition(ostream             &prm_os,           ///< TODOCUMENT
                                    const superposition &prm_superposition ///< TODOCUMENT
                                    ) {
	const streamsize old_precision = prm_os.precision();
	prm_os.precision(30);
	const size_t num_entries = prm_superposition.get_num_entries();
	for (const size_t &entry_ctr : indices( num_entries ) ) {
		const coord    entry_translation = prm_superposition.get_translation_of_index( entry_ctr );
		const rotation entry_rotation    = prm_superposition.get_rotation_of_index(    entry_ctr );
		prm_os << entry_translation.get_x() << " " << entry_translation.get_y() << " " << entry_translation.get_z();
		prm_os << " "  << entry_rotation.get_value<0, 0>();
		prm_os << " "  << entry_rotation.get_value<0, 1>();
		prm_os << " "  << entry_rotation.get_value<0, 2>();
		prm_os << " "  << entry_rotation.get_value<1, 0>();
		prm_os << " "  << entry_rotation.get_value<1, 1>();
		prm_os << " "  << entry_rotation.get_value<1, 2>();
		prm_os << " "  << entry_rotation.get_value<2, 0>();
		prm_os << " "  << entry_rotation.get_value<2, 1>();
		prm_os << " "  << entry_rotation.get_value<2, 2>();
		prm_os << "\n";
	}
	prm_os.precision(old_precision);
}
/// \brief TODOCUMENT
superposition cath::sup::read_superposition(istream &prm_is ///< TODOCUMENT
		                                    ) {
	constexpr size_t CORRECT_NUM_PARTS = 12;
	coord_vec    translations;
	rotation_vec rotations;
	string       line_string;
	while ( getline( prm_is, line_string ) ) {
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
superposition cath::sup::create_pairwise_superposition(const coord_list &prm_coord_list_1,     ///< TODOCUMENT
                                                       const coord_list &prm_coord_list_2,     ///< TODOCUMENT
                                                       const bool       &prm_first_as_base,    ///< TODOCUMENT
                                                       const coord      &prm_base_translation, ///< TODOCUMENT
                                                       const rotation   &prm_base_rotation     ///< TODOCUMENT
                                                       ) {
	const superposition new_superposition(
		vector<superposition::indices_and_coord_lists_type>( {
			std::make_tuple(
				superposition::INDEX_OF_FIRST_IN_PAIRWISE_SUPERPOSITION,
				prm_coord_list_1,
				superposition::INDEX_OF_SECOND_IN_PAIRWISE_SUPERPOSITION,
				prm_coord_list_2
			)
		} ),
		( prm_first_as_base ? superposition::INDEX_OF_FIRST_IN_PAIRWISE_SUPERPOSITION : superposition::INDEX_OF_SECOND_IN_PAIRWISE_SUPERPOSITION ),
		prm_base_translation,
		prm_base_rotation
	);
	return new_superposition;
}

/// \brief A convenience function to superpose one set of coordinates to another
///
/// \relates superposition
void cath::sup::superpose_second_coords_to_first(const coord_list &prm_coord_list_1, ///< The first set of coordinates (to which the others should be superposed)
                                                 coord_list       &prm_coord_list_2  ///< The second set of coordinates (to superpose to the others)
                                                 ) {
	const auto the_sup = create_pairwise_superposition( prm_coord_list_1, prm_coord_list_2 );
	transform( the_sup, superposition::INDEX_OF_SECOND_IN_PAIRWISE_SUPERPOSITION, prm_coord_list_2 );
}

/// \brief A convenience function to create a copy of one set of coordinates superposed to another
///
/// \relates superposition
coord_list cath::sup::superpose_copy_second_coords_to_first(const coord_list &prm_coord_list_1, ///< TODOCUMENT
                                                            coord_list        prm_coord_list_2  ///< TODOCUMENT
                                                            ) {
	superpose_second_coords_to_first( prm_coord_list_1, prm_coord_list_2 );
	return prm_coord_list_2;
}


/// \brief A convenience function to superpose one set of coordinates to another
///
/// \relates superposition
void cath::sup::superpose_second_coords_to_first(const coord_list_vec &prm_coord_list_1, ///< The first set of coordinates (to which the others should be superposed)
                                                 coord_list_vec       &prm_coord_list_2  ///< The second set of coordinates (to superpose to the others)
                                                 ) {
	const auto the_sup = create_pairwise_superposition(
		flatten_coord_lists( prm_coord_list_1 ),
		flatten_coord_lists( prm_coord_list_2 )
	);
	transform( the_sup, superposition::INDEX_OF_SECOND_IN_PAIRWISE_SUPERPOSITION, prm_coord_list_2 );
}

/// \brief A convenience function to create a copy of one set of coordinates superposed to another
///
/// \relates superposition
coord_list_vec cath::sup::superpose_copy_second_coords_to_first(const coord_list_vec &prm_coord_list_1, ///< TODOCUMENT
                                                                coord_list_vec        prm_coord_list_2  ///< TODOCUMENT
                                                                ) {
	superpose_second_coords_to_first( prm_coord_list_1, prm_coord_list_2 );
	return prm_coord_list_2;
}


/// \brief A convenience factory to superpose two coord lists and return the RMSD
///
/// \relates superposition
///
/// Rather than needing to use the more complicated
double cath::sup::calc_pairwise_superposition_rmsd(const coord_list &prm_coord_list_1, ///< TODOCUMENT
                                                   const coord_list &prm_coord_list_2  ///< TODOCUMENT
                                                   ) {
	const superposition new_sup = create_pairwise_superposition(
		prm_coord_list_1,
		prm_coord_list_2
	);
	return calc_rmsd_between_superposed_entries(
		new_sup,
		superposition::INDEX_OF_FIRST_IN_PAIRWISE_SUPERPOSITION,
		prm_coord_list_1,
		superposition::INDEX_OF_SECOND_IN_PAIRWISE_SUPERPOSITION,
		prm_coord_list_2
	);
}

/// \brief TODOCUMENT
///
/// \relates superposition
void cath::sup::check_superposition_is_pairwise(const superposition &prm_superposition ///< TODOCUMENT
                                                ) {
	if (prm_superposition.get_num_entries() != superposition::NUM_ENTRIES_IN_PAIRWISE_SUPERPOSITION) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Superposition is not pairwise"));
	}
}

/// \brief Make an identity superposition for the specified number of entries
///
/// \relates superposition
superposition cath::sup::make_identity_superposition(const size_t &prm_num_entries ///< The number of entries to be superposed
                                                     ) {
	return {
		coord_vec   ( prm_num_entries, ORIGIN_COORD           ),
		rotation_vec( prm_num_entries, rotation::IDENTITY_ROTATION() )
	};
}
