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

#include <boost/lexical_cast.hpp>
#include <boost/log/trivial.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/irange.hpp>

#include "common/algorithm/transform_build.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "common/size_t_literal.hpp"
//#include "common/temp_check_offset_1.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "ssap/context_res.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"
#include "structure/residue_name.hpp"

#include <algorithm> // for max, min
#include <iterator>  // for end, begin, etc
#include <sstream>   // for string, etc

using namespace boost::log;
using namespace cath;
using namespace cath::common;
using namespace cath::geom;
using namespace std;

using boost::irange;
using boost::numeric_cast;
using boost::range::find_if;
using boost::lexical_cast;

/// \brief TODOCUMENT
protein::protein(const string      &arg_title,   ///< TODOCUMENT
                 const residue_vec &arg_residues ///< TODOCUMENT
                 ) : title    ( arg_title    ),
                     residues ( arg_residues ) {
}

/// \brief TODOCUMENT
void protein::set_title(const string &arg_title ///< TODOCUMENT
                        ) {
	title = arg_title;
}

/// \brief TODOCUMENT
void protein::set_residues(const residue_vec &arg_residues ///< TODOCUMENT
                           ) {
	residues = arg_residues;
}

/// \brief TODOCUMENT
void protein::set_sec_strucs(const sec_struc_vec &arg_sec_strucs ///< TODOCUMENT
                             ) {
	sec_strucs = arg_sec_strucs;
}

/// \brief TODOCUMENT
string protein::get_title() const {
	return title;
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
protein cath::build_protein(const residue_vec &arg_residues ///< TODOCUMENT
                            ) {
	protein new_protein;
	new_protein.set_residues(arg_residues);
	return new_protein;
}

/// \brief TODOCUMENT
protein cath::build_protein(const residue_vec       &arg_residues,  ///< TODOCUMENT
                            const sec_struc_vec &arg_sec_strucs ///< TODOCUMENT
                            ) {
	protein new_protein;
	new_protein.set_residues(arg_residues);
	new_protein.set_sec_strucs(arg_sec_strucs);
	return new_protein;
}

/// \brief Retrieve the PDB residue name for the residue with the specified index in the specified protein
///
/// \relates protein
residue_name cath::get_pdb_residue_name_of_index(const protein &arg_protein,      ///< The protein containing the residue whose name should be returned
                                                 const size_t  &arg_residue_index ///< The index of the residue within the protein
                                                 ) {
	return arg_protein.get_residue_ref_of_index( arg_residue_index ).get_pdb_residue_name();
}

/// \brief Retrieve the PDB residue name for the residue with the specified index in the specified protein
///
/// \relates protein
string cath::get_pdb_residue_name_string_of_index(const protein &arg_protein,      ///< The protein containing the residue whose name should be returned
                                                  const size_t  &arg_residue_index ///< The index of the residue within the protein
                                                  ) {
	return lexical_cast<string>( get_pdb_residue_name_of_index( arg_protein, arg_residue_index ) );
}

/// \brief TODOCUMENT
///
/// \relates protein
size_t cath::get_index_of_pdb_residue_name(const protein      &arg_protein,     ///< TODOUCMENT
                                           const residue_name &arg_residue_name ///< TODOUCMENT
                                           ) {
	// Return the distance from the start of the protein to the result of finding a residue
	// that matches the specified residue name

	const size_t result_index = numeric_cast<size_t>( distance(
		common::cbegin( arg_protein ),
		find_if(
			arg_protein,
			[&] (const residue &x) { return residue_matches_residue_name( x, arg_residue_name ); }
		)
	) );
	if ( result_index >= arg_protein.get_length() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Unable to find residue with residue name \""
			+ lexical_cast<string>( arg_residue_name )
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
char cath::get_amino_acid_letter_of_index(const protein &arg_protein,      ///< TODOCUMENT
                                          const size_t  &arg_residue_index ///< TODOCUMENT
                                          ) {
	return get_amino_acid_of_index( arg_protein, arg_residue_index ).get_letter();
}

/// \brief Retrieve a list of all the PDB residue names of the residues in the specified protein
///
/// \relates protein
residue_name_vec cath::get_residue_names(const protein &arg_protein ///< The protein containing the residues whose names should be returned
                                         ) {
	return transform_build<residue_name_vec>(
		arg_protein,
		[] (const residue &x) { return x.get_pdb_residue_name(); }
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
void cath::label_residues_with_sec_strucs(protein &arg_protein,   ///< The protein object to modify
                                          ostream &/*arg_stderr*/ ///< TODOCUMENT
                                          ) {
	/// \brief TODOCUMENT
	const size_t num_residues = arg_protein.get_length();

	// Label all residues with secondary structure label of sec_struc_type::COIL
	wipe_sec_strucs_of_residues( arg_protein );

	// Loop over all the secondary structures
	const size_t num_sec_strucs = arg_protein.get_num_sec_strucs();
	for (const size_t &sec_struc_ctr : irange( 0_z, num_sec_strucs ) ) {

		// Grab data from the secondary structure
		const auto &my_sec_struc    = arg_protein.get_sec_struc_ref_of_index(sec_struc_ctr);
		const auto  sec_struc_start = my_sec_struc.get_start_residue_num();
		const auto  sec_struc_stop  = my_sec_struc.get_stop_residue_num();

		// Sanity check the start and stop of the secondary structure
		if ( sec_struc_start > sec_struc_stop ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception(
				"Secondary structure has a start of " + lexical_cast<string>(sec_struc_start) +
				", which is greater than its stop "   + lexical_cast<string>(sec_struc_stop) +
				" for protein " + arg_protein.get_title()
			));
		}
		if ( sec_struc_stop > num_residues ) {
			BOOST_LOG_TRIVIAL( warning ) << "Ignoring (part of) secondary structure because it has a stop of " << sec_struc_stop
			                             << ", which is greater than the number of residues " << num_residues
			                             << " for protein " << arg_protein.get_title()
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
				" for protein " + arg_protein.get_title()
			));
		}
		if ( there_is_a_preceding_ss && sec_struc_start <= prev_stop ) {
			BOOST_LOG_TRIVIAL( warning ) << "Secondary structure starts at residue number " << lexical_cast<string>(sec_struc_start)
			                             << ", which overlaps with the end of the previous secondary structure at residue number " << lexical_cast<string>(prev_stop)
			                             << " for protein " << arg_protein.get_title()
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
			+ lexical_cast<string>(num_sec_strucs)
			+ " secondary structures"
		));
	}
	if (arg_dest_sec_struc_index >= num_sec_strucs) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Unable to calculate inter-sec_struc vector because source index "
			+ lexical_cast<string>(arg_src_sec_struc_index)
			+ " is out of range in a protein with "
			+ lexical_cast<string>(num_sec_strucs)
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

