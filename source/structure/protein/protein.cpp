/// \file
/// \brief The protein class definitions

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

#include "protein.hpp"

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/log/trivial.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/algorithm_ext/for_each.hpp>
#include <boost/range/irange.hpp>

#include "biocore/residue_id.hpp"
#include "chopping/region/regions_limiter.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "exception/not_implemented_exception.hpp"
#include "ssap/context_res.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"

#include <algorithm> // for max, min
#include <iterator>  // for end, begin, etc
#include <sstream>   // for string, etc

using namespace boost::log;
using namespace cath;
using namespace cath::chop;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace std;

using boost::algorithm::any_of;
using boost::irange;
using boost::lexical_cast;
using boost::numeric_cast;
using boost::range::find_if;
using boost::range::for_each;

/// \brief TODOCUMENT
protein::protein(name_set    arg_name_set, ///< TODOCUMENT
                 residue_vec arg_residues  ///< TODOCUMENT
                 ) : the_name_set { std::move( arg_name_set ) },
                     residues     { std::move( arg_residues ) } {
}

/// \brief TODOCUMENT
protein & protein::set_name_set(name_set arg_name_set ///< TODOCUMENT
                                ) {
	the_name_set = std::move( arg_name_set );
	return *this;
}

/// \brief TODOCUMENT
protein & protein::set_residues(residue_vec arg_residues ///< TODOCUMENT
                                ) {
	residues = std::move( arg_residues );
	return *this;
}

/// \brief TODOCUMENT
protein & protein::set_sec_strucs(sec_struc_vec arg_sec_strucs ///< TODOCUMENT
                                  ) {
	sec_strucs = std::move( arg_sec_strucs );
	return *this;
}

/// \brief TODOCUMENT
name_set & protein::get_name_set() {
	return the_name_set;
}

/// \brief TODOCUMENT
const name_set & protein::get_name_set() const {
	return the_name_set;
}

/// \brief TODOCUMENT
sec_struc & protein::get_sec_struc_ref_of_index(const size_t &arg_index ///< TODOCUMENT
                                                ) {
	check_sec_struc_is_valid(arg_index);
	return sec_strucs[arg_index];
}

/// \brief TODOCUMENT
const sec_struc & protein::get_sec_struc_ref_of_index(const size_t &arg_index ///< TODOCUMENT
                                                      ) const {
	check_sec_struc_is_valid(arg_index);
	return sec_strucs[arg_index];
}

/// \brief TODOCUMENT
size_t protein::get_num_sec_strucs() const {
	return sec_strucs.size();
}

/// \brief TODOCUMENT
protein::iterator protein::begin() {
	return std::begin( residues );
}

/// \brief TODOCUMENT
protein::iterator protein::end() {
	return std::end( residues );
}

/// \brief TODOCUMENT
protein::const_iterator protein::begin() const {
	return common::cbegin( residues );
}

/// \brief TODOCUMENT
protein::const_iterator protein::end() const {
	return common::cend( residues );
}
/// \brief TODOCUMENT
protein::sec_struc_crange protein::get_sec_strucs() const {
	return {
		common::cbegin( sec_strucs ),
		common::cend  ( sec_strucs )
	};
}

/// \brief TODOCUMENT
protein cath::build_protein(residue_vec arg_residues ///< TODOCUMENT
                            ) {
	protein new_protein;
	new_protein.set_residues( std::move( arg_residues ) );
	return new_protein;
}

/// \brief TODOCUMENT
protein cath::build_protein(residue_vec   arg_residues,  ///< TODOCUMENT
                            sec_struc_vec arg_sec_strucs ///< TODOCUMENT
                            ) {
	protein new_protein;
	new_protein.set_residues  ( std::move( arg_residues   ) );
	new_protein.set_sec_strucs( std::move( arg_sec_strucs ) );
	return new_protein;
}

/// \brief Retrieve the PDB residue name for the residue with the specified index in the specified protein
///
/// \relates protein
residue_id cath::get_pdb_residue_id_of_index(const protein &arg_protein,      ///< The protein containing the residue whose name should be returned
                                             const size_t  &arg_residue_index ///< The index of the residue within the protein
                                             ) {
	return arg_protein.get_residue_ref_of_index( arg_residue_index ).get_pdb_residue_id();
}

/// \brief Retrieve the PDB residue name for the residue with the specified index in the specified protein
///
/// \relates protein
string cath::get_pdb_residue_id_string_of_index(const protein &arg_protein,      ///< The protein containing the residue whose name should be returned
                                                const size_t  &arg_residue_index ///< The index of the residue within the protein
                                                ) {
	return to_string( get_pdb_residue_id_of_index( arg_protein, arg_residue_index ) );
}

/// \brief TODOCUMENT
///
/// \relates protein
size_t cath::get_index_of_pdb_residue_id(const protein    &arg_protein,   ///< TODOUCMENT
                                         const residue_id &arg_residue_id ///< TODOUCMENT
                                         ) {
	// Return the distance from the start of the protein to the result of finding a residue
	// that matches the specified residue name

	const size_t result_index = numeric_cast<size_t>( distance(
		common::cbegin( arg_protein ),
		find_if(
			arg_protein,
			[&] (const residue &x) { return residue_matches_residue_id( x, arg_residue_id ); }
		)
	) );
	if ( result_index >= arg_protein.get_length() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Unable to find residue with residue name \""
			+ to_string( arg_residue_id )
			+ "\" in protein"
		));
	}
	return result_index;
}

/// \brief TODOCUMENT
amino_acid_vec cath::get_amino_acid_list(const protein &arg_protein ///< TODOCUMENT
                                         ) {
	amino_acid_vec amino_acids;
	amino_acids.reserve( arg_protein.get_length() );
	for (const residue &the_residue : arg_protein) {
		amino_acids.push_back( the_residue.get_amino_acid() );
	}
	return amino_acids;
}

/// \brief TODOCUMENT
///
/// \relates protein
amino_acid cath::get_amino_acid_of_index(const protein &arg_protein,      ///< TODOCUMENT
                                         const size_t  &arg_residue_index ///< TODOCUMENT
                                         ) {
	return arg_protein.get_residue_ref_of_index( arg_residue_index ).get_amino_acid();
}

/// \brief TODOCUMENT
///
/// \relates protein
char cath::get_amino_acid_letter_of_index_tolerantly(const protein &arg_protein,      ///< TODOCUMENT
                                                     const size_t  &arg_residue_index ///< TODOCUMENT
                                                     ) {
	return get_amino_acid_of_index( arg_protein, arg_residue_index ).get_letter_tolerantly();
}

/// \brief Retrieve a list of all the PDB residue names of the residues in the specified protein
///
/// \relates protein
residue_id_vec cath::get_residue_ids(const protein &arg_protein ///< The protein containing the residues whose names should be returned
                                     ) {
	return transform_build<residue_id_vec>(
		arg_protein,
		[] (const residue &x) { return x.get_pdb_residue_id(); }
	);
}

/// \brief Set the specified accessibility values in the specified protein
///
/// \relates protein
void cath::set_accessibilities(protein        &arg_protein,        ///< The protein to modify
                               const size_vec &arg_accessibilities ///< The accessibilities to set
                               ) {
	if ( arg_accessibilities.size() != arg_protein.get_length() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Cannot set "
			+ std::to_string( arg_accessibilities.size() )
			+ " accessibility values for a protein of length "
			+ std::to_string( arg_protein.get_length() )
		));
	}
	for_each(
		arg_protein,
		arg_accessibilities,
		[] (residue &the_residue, const size_t &accessibility) {
			the_residue.set_access( accessibility );
		}
	);
}

/// \brief Set the specified accessibility values in the specified protein,
///        converting them to unsigned integer values first
///
/// \relates protein
void cath::set_accessibilities(protein        &arg_protein,        ///< The protein to modify
                               const doub_vec &arg_accessibilities ///< The accessibilities to set
                               ) {
	if ( any_of( arg_accessibilities, [] (const double &x) { return ! boost::math::isfinite( x ) || x < 0; } ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot set accessibilities given negative/infinite/NaN values"));
	}
	set_accessibilities(
		arg_protein,
		transform_build<size_vec>(
			arg_accessibilities,
			[] (const double &x) {
				return numeric_cast<size_t>( round( x ) );
			}
		)
	);
}

/// \brief Set the specified sec_struc_type values in the specified protein
///
/// \relates protein
void cath::set_sec_struc_types(protein                  &arg_protein,        ///< The protein to modify
                               const sec_struc_type_vec &arg_sec_struc_types ///< The sec_struc_type values to set
                               ) {
	if ( arg_sec_struc_types.size() != arg_protein.get_length() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Cannot set "
			+ std::to_string( arg_sec_struc_types.size() )
			+ " sec_struc types for a protein of length "
			+ std::to_string( arg_protein.get_length() )
		));
	}
	for_each(
		arg_protein,
		arg_sec_struc_types,
		[] (residue &the_residue, const sec_struc_type &the_sec_struc) {
			the_residue.set_sec_struc_type( the_sec_struc );
		}
	);
}

/// \brief For all residues, wipe secondary structure numbers and set the label to sec_struc_type::COIL
///
/// \relatesalso protein
void cath::wipe_sec_strucs_of_residues(protein &arg_protein ///< The protein object to modify
                                       ) {
	for (residue &the_residue : arg_protein) {
		wipe_secondary_structure( the_residue );
	}
}

/// \brief Label residues with appropriate secondary structure numbers
///
/// \relatesalso protein
void cath::label_residues_with_sec_strucs(protein               &arg_protein,   ///< The protein object to modify
                                          const ostream_ref_opt &/*arg_stderr*/ ///< An optional reference to an ostream to which any logging should be performed
                                          ) {
	/// \brief TODOCUMENT
	const size_t num_residues = arg_protein.get_length();

	// Label all residues with secondary structure label of sec_struc_type::COIL
	wipe_sec_strucs_of_residues( arg_protein );

	// Loop over all the secondary structures
	const size_t num_sec_strucs = arg_protein.get_num_sec_strucs();
	for (const size_t &sec_struc_ctr : indices( num_sec_strucs ) ) {

		// Grab data from the secondary structure
		const auto &my_sec_struc    = arg_protein.get_sec_struc_ref_of_index(sec_struc_ctr);
		const auto  sec_struc_start = my_sec_struc.get_start_residue_num();
		const auto  sec_struc_stop  = my_sec_struc.get_stop_residue_num();

		// Sanity check the start and stop of the secondary structure
		if ( sec_struc_start > sec_struc_stop ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception(
				"Secondary structure has a start of " + lexical_cast<string>(sec_struc_start) +
				", which is greater than its stop "   + lexical_cast<string>(sec_struc_stop) +
				" for protein " + get_domain_or_specified_or_name_from_acq( arg_protein )
			));
		}
		if ( sec_struc_stop > num_residues ) {
			BOOST_LOG_TRIVIAL( warning ) << "Ignoring (part of) secondary structure because it has a stop of " << sec_struc_stop
			                             << ", which is greater than the number of residues " << num_residues
			                             << " for protein " << arg_protein.get_name_set()
			                             << " (this is probably because the secondary structure has been read from a sec file"
			                             << ", which refers to residues by their sequential order"
			                             << ", whereas the residues have been read from a dssp or wolf file"
			                             << ", which drops some residues"
			                             << ", due to atom records removed by splitchains and for other reasons)";
		}

		// Get details for any preceding secondary structure
		const bool   there_is_a_preceding_ss = (sec_struc_ctr > 0);
		const size_t prev_start = ( there_is_a_preceding_ss ? arg_protein.get_sec_struc_ref_of_index( sec_struc_ctr - 1 ).get_start_residue_num() : 0 );
		const size_t prev_stop  = ( there_is_a_preceding_ss ? arg_protein.get_sec_struc_ref_of_index( sec_struc_ctr - 1 ).get_stop_residue_num()  : 0 );

		// If there is a previous secondary structure and its start is no earlier than this one, then this
		// is a serious error so throw an exception
		if (there_is_a_preceding_ss && sec_struc_start <= prev_start) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception(
				"Secondary structure has a start of " + lexical_cast<string>(sec_struc_start) +
				", which is not greater than the preceding secondary structure's start of " + lexical_cast<string>(prev_start) +
				" for protein " + get_domain_or_specified_or_name_from_acq( arg_protein )
			));
		}
		if ( there_is_a_preceding_ss && sec_struc_start <= prev_stop ) {
			BOOST_LOG_TRIVIAL( warning ) << "Secondary structure starts at residue number " << lexical_cast<string>(sec_struc_start)
			                             << ", which overlaps with the end of the previous secondary structure at residue number " << lexical_cast<string>(prev_stop)
			                             << " for protein " << arg_protein.get_name_set()
			                             << " - will use the previous secondary structure to label residue(s) within overlapping region.";
		}

		// If there is an overlap with the previous, then start updating from the from the first residue after the end of the previous
		const size_t update_start__os1_uncapped = there_is_a_preceding_ss ? max( prev_stop + 1, sec_struc_start ) : sec_struc_start;
		const size_t update_start__os1          = min( num_residues + 1, update_start__os1_uncapped );
		const size_t update_stop__os1           = min( num_residues + 1, sec_struc_stop + 1         );

		// Loop over the residues that are covered by this secondary structure, labelling each
		for (const size_t &residue_ctr__os1 : irange( update_start__os1, update_stop__os1 ) ) {
			residue &my_residue = get_residue_ref_of_index__offset_1( arg_protein, residue_ctr__os1 );
			my_residue.set_residue_sec_struc_number( sec_struc_ctr + 1       );
			my_residue.set_sec_struc_type          ( my_sec_struc.get_type() );
		}
	}
}

/// \brief TODOCUMENT
///
/// \relates protein
///
/// The way that this selects the coordinate of the anchor point is as follows:
///  * vector from this secondary structure's midpoint to the next secondary structure's midpoint
///  * unless this is the last secondary structure, in which case the previous secondary structure's used instead
coord cath::calculate_inter_sec_struc_vector(const protein            &arg_protein,             ///< The protein to be queried
                                             const size_t &arg_src_sec_struc_index, ///< The index of the source secondary structure in the protein
                                             const size_t &arg_dest_sec_struc_index ///< The index of the destination secondary structure in the protein
                                             ) {
	// Sanity check the inputs
	const size_t num_sec_strucs = arg_protein.get_num_sec_strucs();
//	if (arg_src_sec_struc_index == arg_dest_sec_struc_index) {
//		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to calculate inter-sec_struc vector between a sec_struc and itself"));
//	}
	if (arg_src_sec_struc_index >= num_sec_strucs) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Unable to calculate inter-sec_struc vector because source index "
			+ lexical_cast<string>(arg_src_sec_struc_index)
			+ " is out of range in a protein with "
			+ lexical_cast<string>( num_sec_strucs )
			+ " secondary structures"
		));
	}
	if (arg_dest_sec_struc_index >= num_sec_strucs) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Unable to calculate inter-sec_struc vector because source index "
			+ lexical_cast<string>(arg_src_sec_struc_index)
			+ " is out of range in a protein with "
			+ lexical_cast<string>( num_sec_strucs )
			+ " secondary structures"
		));
	}

	// If the source and destination sec_struc indices are equal, then just return vector (0, 0, 0)
	if (arg_src_sec_struc_index == arg_dest_sec_struc_index) {
		return coord::ORIGIN_COORD;
	}

	// Grab the two sec_strucs along with a third, anchor sec_struc
	const ptrdiff_t  anchor_sec_struc_offset = (arg_src_sec_struc_index + 1 != arg_protein.get_num_sec_strucs()) ? 1 : -1;
	const sec_struc &src_sec_struc           = arg_protein.get_sec_struc_ref_of_index( arg_src_sec_struc_index  );
	const sec_struc &dest_sec_struc          = arg_protein.get_sec_struc_ref_of_index( arg_dest_sec_struc_index );
	const size_t     sec_struc_index         = numeric_cast<size_t>( numeric_cast<ptrdiff_t>( arg_src_sec_struc_index )  + anchor_sec_struc_offset );
	const sec_struc &anchor_sec_struc        = arg_protein.get_sec_struc_ref_of_index( sec_struc_index );

	// Use the two sec_strucs and the anchor sec_struc to calculate the inter-sec_struc vector
	return calculate_inter_sec_struc_vector( src_sec_struc, dest_sec_struc, anchor_sec_struc );
}

/// \brief TODOCUMENT
///
/// \relates protein
coord cath::view_vector(const protein &arg_protein,    ///< TODOCUMENT
                        const size_t  &arg_from_index, ///< TODOCUMENT
                        const size_t  &arg_to_index    ///< TODOCUMENT
                        ) {
	return view_vector_of_residue_pair(
		arg_protein.get_residue_ref_of_index( arg_from_index ),
		arg_protein.get_residue_ref_of_index( arg_to_index   )
	);
}

/// \brief Get the indices of the residues in the specified protein within the specified regions
///
/// \relates protein
size_vec cath::get_indices_of_residues_within_regions(const protein        &arg_protein, ///< The protein to query
                                                      const region_vec_opt &arg_regions  ///< The regions to query
                                                      ) {
	size_vec result;
	regions_limiter the_limiter{ arg_regions };
	for (const size_t &res_idx : indices( arg_protein.get_length() ) ) {
		if ( the_limiter.update_residue_is_included( get_pdb_residue_id_of_index( arg_protein, res_idx ) ) ) {
			result.push_back( res_idx );
		}
	}
	return result;
}

/// \brief Restrict the specified, existing protein to the specified regions
///
/// \relates protein
void cath::restrict_to_regions(protein              &/*arg_protein*/, ///< The initial protein to restrict
                               const region_vec_opt &/*arg_regions*/  ///< The regions to which the resulting protein should be restricted
                               ) {

	BOOST_THROW_EXCEPTION(not_implemented_exception("Cannot yet restrict to regions with this combination of input files"));

// 	residue_vec residues_to_keep;
// 	residues_to_keep.reserve( arg_protein.get_length() );

// 	regions_limiter the_limiter{ arg_regions };

// 	for (const residue &res : arg_protein) {
// 		if ( the_limiter.update_residue_is_included( res.get_pdb_residue_id() ) ) {
// 			residues_to_keep.push_back( res );
// 		}
// 	}

// 	// A vector to populate with the sec_strucs that are to be kept
// 	const size_t num_sec_strucs = arg_protein.get_num_sec_strucs();
// 	sec_struc_vec sec_strucs_to_keep;
// 	sec_strucs_to_keep.reserve( num_sec_strucs );

// 	// A vector of the old numbers of each of the kept sec_strucs indexed by their new numbers.
// 	size_vec sec_struc_index_conv;
// 	sec_struc_index_conv.reserve( num_sec_strucs );

// 	for (const size_t &sec_struc_ctr : indices( num_sec_strucs ) ) {
// 		const sec_struc &my_sec_struc = arg_protein.get_sec_struc_ref_of_index(sec_struc_ctr);

// 		// Look to see if we want to exclude this sec_struc
// 		const bool found = any_of(
// 			arg_clique_starts_and_ends,
// 			[&] (const size_size_pair &x) {
// 				const size_t &a_start = x.first;
// 				const size_t &a_end   = x.second;

// 				// If the secondary structure's start or end is within the range, then mark it at as one to exclude
// 				// (Should this code also exclude secondary structures that straddle the range?
// 				//  ie (my_sec_struc.from <= a_end) || (my_sec_struc.to >= a_start) ? )
// 				return (
// 					( my_sec_struc.get_start_residue_num() >= a_start && my_sec_struc.get_start_residue_num() <= a_end )
// 					||
// 					( my_sec_struc.get_stop_residue_num()  >= a_start && my_sec_struc.get_stop_residue_num()  <= a_end ) );
// 			}
// 		);

// 		// If it's not designated to be removed, then keep
// 		if ( ! found ) {
// 			// Record the old sec_struc number according to the new sec_struc number
// 			// (the new sec_struc number is implied because sec_struc_index_conv is populated synchronously with sec_strucs_to_keep)
// 			sec_struc_index_conv.push_back(sec_struc_ctr);

// 			// Add the new sec_struc to sec_strucs_to_keep, first updating its sec_struc number to the new value
// 			sec_struc new_sec_struc(my_sec_struc);
// 			sec_strucs_to_keep.push_back(my_sec_struc);
// 		}
// 	}

// 	// Update pair data
// 	const size_t num_sec_strucs_to_keep = sec_strucs_to_keep.size();
// 	for (const size_t &new_sec_struc_ctr_i : indices( num_sec_strucs_to_keep ) ) {
// 		const size_t    &old_sec_struc_index_i = sec_struc_index_conv[new_sec_struc_ctr_i];
// 		const sec_struc &old_sec_struc         = arg_protein.get_sec_struc_ref_of_index( old_sec_struc_index_i );

// 		sec_struc_planar_angles_vec new_planar_angles;
// 		new_planar_angles.reserve(num_sec_strucs_to_keep);

// 		for (const size_t &old_sec_struc_index_j : sec_struc_index_conv) {
// 			const sec_struc_planar_angles old_planar_angles = old_sec_struc.get_planar_angles_of_index( old_sec_struc_index_j );

// 			// The previous code appeared to be wiping all planar_angles here!
// 			// Not sure what was going on there.
// 			//                                        Tony Lewis, 30th September 2012
// 			new_planar_angles.push_back(old_planar_angles);
// 		}

// 		sec_struc &new_sec_struc = sec_strucs_to_keep[new_sec_struc_ctr_i];
// 		new_sec_struc.set_planar_angles(new_planar_angles);
// 	}
// 	arg_protein.set_sec_strucs(sec_strucs_to_keep);



// 	/////
// 	// 2. Process the residues that are to be removed
// 	/////

// 	const size_t protein_length = arg_protein.get_length();
// 	residue_vec residues_to_keep;
// 	residues_to_keep.reserve(protein_length);

// 	// Look for residues to exclude
// 	for (const size_t &residue_ctr : indices( protein_length ) ) {
// 		const residue &my_residue = arg_protein.get_residue_ref_of_index(residue_ctr);
// 		bool           found      = false;

// 		// Look to see if we want to exclude this residue
// 		for (const size_size_pair &clique_start_and_end : arg_clique_starts_and_ends) {
// 			const size_t &a_start = clique_start_and_end.first;
// 			const size_t &a_end   = clique_start_and_end.second;

// 			/// Problem in which residues names are compared numerically
// 			/// --------------------------------------------------------
// 			///
// 			/// \todo Fix this! this is comparing residue names without handling insert-codes and without any consideration
// 			/// that the numbers in residue names can be completely out of order.
// 			///
// 			/// Scanning a chain with lots of insert-codes, eg:
// 			///
// 			///     > cathedral_scan -d 2qriB -c /cath/data/v3_5_0/release_data/Class1-4.s35.cathedral.library --clique --dir /cath/data/v3_5_0/grath
// 			///
// 			/// ...produces files like this:
// 			///
// 			///     > cat 2qriB1nepA00.clique
// 			///     5
// 			///      2   21   30  1    5    8
// 			///      3   36   41  7   62   65
// 			///      8   78   83  6   52   59
// 			///     29  233  236  5   35   44
// 			///     30  241  250  3   16   21
// 			///
// 			/// ...but then CATHEDRAL can hardly be expected to get this right because the grath file is wrong:
// 			///
// 			///     > head -n 6 /cath/data/v3_5_0/grath/2qriB.gth
// 			///     2qriB
// 			///     32
// 			///     6 11
// 			///     21 30
// 			///     36 41
// 			///     43 46
// 			///
// 			/// Perhaps it'd be better if everything used sequential numbers rather than insert codes.
// 			///
// 			/// \todo Is the initial `pdb_number( my_residue ) != 0` check just an error?
// 			if ( ( pdb_number( my_residue ) != 0 ) && pdb_number( my_residue ) >= numeric_cast<int>( a_start ) && pdb_number( my_residue ) <= numeric_cast<int>( a_end ) ) {
// 				found = true;
// 				break;
// 			}
// 		}

// 		// If it's not designated to be removed, then keep
// 		if (!found) {
// 			residues_to_keep.push_back(my_residue);
// 		}
// 	}

// //	const size_t num_residues_to_keep = residues_to_keep.size();
// //
// //	for (const size_t &i : indices( num_residues_to_keep ) ) {
// //		residue &newp = residues_to_keep[ i ];
// //
// //		// Disulfide information (not sure it's important)
// //		if ( islower(newp.get_amino_acid()) ) {
// //			cerr << "IMPORTANT NOTICE : ******** THIS BIT OF CODE IS NOT COMPLETELY USELESS - HURRAH - PLEASE TELL SOMEONE *******" << endl;
// //			newp.set_amino_acid('C');
// //		}
// //	}

// 	arg_protein.set_residues( residues_to_keep );

// 	label_residues_with_sec_strucs( arg_protein, arg_stderr );
}

/// \brief TODOCUMENT
///
/// \relates protein
protein cath::restrict_to_regions_copy(protein               arg_protein, ///< TODOCUMENT
                                       const region_vec_opt &arg_regions  ///< TODOCUMENT
                                       ) {
	restrict_to_regions( arg_protein, arg_regions );
	return arg_protein;
}

/// \brief Get the domain or specified name from the specified protein
///
/// \relates protein
string cath::get_domain_or_specified_or_name_from_acq(const protein &arg_protein ///< The protein to query
                                                      ) {
	return get_domain_or_specified_or_name_from_acq( arg_protein.get_name_set() );
}
