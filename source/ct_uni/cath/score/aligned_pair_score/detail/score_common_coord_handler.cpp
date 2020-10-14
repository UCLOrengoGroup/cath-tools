/// \file
/// \brief The score_common_coord_handler class definitions

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

#include "score_common_coord_handler.hpp"

#include <boost/lexical_cast.hpp>

#include "cath/alignment/alignment_coord_extractor.hpp"
#include "cath/alignment/common_residue_selection_policy/common_residue_select_all_policy.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"
#include "cath/common/less_than_helper.hpp"
#include "cath/structure/geometry/coord.hpp"
#include "cath/structure/geometry/coord_list.hpp"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::geom;
using namespace cath::score::detail;
using namespace std;

using boost::lexical_cast;

/// \brief Private getter for comm_res_seln_pol_ptr
const common_residue_selection_policy & score_common_coord_handler::get_comm_res_seln_pol() const {
	return *comm_res_seln_pol_ptr;
}

/// \brief Private getter for comm_atom_seln_pol_ptr
const common_atom_selection_policy & score_common_coord_handler::get_comm_atom_seln_pol() const {
	return *comm_atom_seln_pol_ptr;
}

/// \brief Get the descriptions strings from the common_residue_selection_policy and common_atom_selection_policy,
///        replacing each with an empty string if it is the default policy
///
/// \todo Change the callers of this to use short_name_suffixes() instead and then remove this
str_str_pair score_common_coord_handler::get_policy_description_strings() const {
	const string res_long_suffix  = is_default_policy( get_comm_res_seln_pol()  ) ? "" : get_comm_res_seln_pol().get_descriptive_name();
	const string atom_long_suffix = is_default_policy( get_comm_atom_seln_pol() ) ? "" : get_comm_atom_seln_pol().get_descriptive_name();
	return make_pair(res_long_suffix, atom_long_suffix);
}

/// \brief Ctor for score_common_coord_handler that allows the caller to specify the common_residue_selection_policy to be used
score_common_coord_handler::score_common_coord_handler(const common_residue_selection_policy &prm_comm_res_seln_pol, ///< The policy to use for selecting common residues
                                                       const common_atom_selection_policy    &prm_atom_seln_pol      ///< The policy to use for selecting common atoms
                                                       ) : comm_res_seln_pol_ptr  ( prm_comm_res_seln_pol.clone() ),
                                                           comm_atom_seln_pol_ptr ( prm_atom_seln_pol.clone()     ) {
}

/// \brief Return a string to append to the end of short names to provide more detail on how the common
///        coordinates are generated (or an empty string if this is just the default)
string score_common_coord_handler::short_suffix_string() const {
	const str_str_pair policy_strings = get_policy_description_strings();
	const string &res_pol_string  = policy_strings.first;
	const string &atom_pol_string = policy_strings.second;
	return ( res_pol_string.empty()  ? "" : "." + res_pol_string  )
	     + ( atom_pol_string.empty() ? "" : "." + atom_pol_string );
}

/// \brief Return a string to append to the end of long names to provide more detail on how the common
///        coordinates are generated (or an empty string if this is just the default)
string score_common_coord_handler::long_suffix_string() const {
	const str_str_pair policy_strings = get_policy_description_strings();
	const string &res_pol_string  = policy_strings.first;
	const string &atom_pol_string = policy_strings.second;
	return   (   res_pol_string.empty() &&   atom_pol_string.empty() ) ? ""
	       : ( ! res_pol_string.empty() && ! atom_pol_string.empty() ) ? " (" + res_pol_string + "; " + atom_pol_string + ")"
	                                                                   : " (" + res_pol_string        + atom_pol_string + ")";
}

/// \brief Return a free text description in brackets to insert into the descriptions to provide more detail
///        on how the common coordinates are generated (or an empty string if this is just the default)
string score_common_coord_handler::description_brackets_string() const {
	const str_str_pair policy_strings = get_policy_description_strings();
	const string &res_pol_string  = policy_strings.first;
	const string &atom_pol_string = policy_strings.second;
	return   (   res_pol_string.empty() &&   atom_pol_string.empty() ) ? ""
	       : ( ! res_pol_string.empty() && ! atom_pol_string.empty() ) ? " (as filtered through the " + res_pol_string + " and " + atom_pol_string + " selection policies)"
	                                                                   : " (as filtered through the " + res_pol_string           + atom_pol_string + " selection policy)";
}

/// \brief TODOCUMENT
str_bool_pair_vec score_common_coord_handler::short_name_suffixes() const {
	return {
		{ get_comm_res_seln_pol().get_descriptive_name(),  ! is_default_policy( get_comm_res_seln_pol()  ) },
		{ get_comm_atom_seln_pol().get_descriptive_name(), ! is_default_policy( get_comm_atom_seln_pol() ) }
	};
}

/// \brief Get the common coordinates between residues, as specified by the common_residue_selection_policy
///
/// \post The two coord_lists will have matching sizes,
///       else out_of_range_exception will be thrown
///
/// Note that this may contain more than one coordinate pair per residue, depending on the common atom policy
coord_list_coord_list_pair score_common_coord_handler::get_common_coords(const alignment &prm_alignment, ///< The pair alignment to be scored
                                                                         const protein   &prm_protein_a, ///< The protein associated with the first  half of the alignment
                                                                         const protein   &prm_protein_b  ///< The protein associated with the second half of the alignment
                                                                         ) const {
	// Extract the common coordinates to be chosen
	const pair<coord_list, coord_list> common_coords = alignment_coord_extractor::get_common_coords(
		prm_alignment,
		prm_protein_a,
		prm_protein_b,
		get_comm_res_seln_pol(),
		get_comm_atom_seln_pol()
	);

	// Check that the two lists' sizes match
	if ( common_coords.first.size() != common_coords.second.size() ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception("Lists of common coordinates have mismatching sizes"));
	}

	// Return the list of common coordinates
	return common_coords;
}

/// \brief Get the common coordinates between residues, as specified by the common_residue_selection_policy
///
/// \post The two coord_lists will have matching sizes,
///       else out_of_range_exception will be thrown
///
/// Note that this may contain more than one coordinate pair per residue, depending on the common atom policy
pair<coord_list_vec, coord_list_vec> score_common_coord_handler::get_common_coords_by_residue(const alignment &prm_alignment, ///< The pair alignment to be scored
                                                                                              const protein   &prm_protein_a, ///< The protein associated with the first  half of the alignment
                                                                                              const protein   &prm_protein_b  ///< The protein associated with the second half of the alignment
                                                                                              ) const {
	// Extract the common coordinates to be chosen
	const pair<coord_list_vec, coord_list_vec> common_coords = alignment_coord_extractor::get_common_coords_by_residue(
		prm_alignment,
		prm_protein_a,
		prm_protein_b,
		get_comm_res_seln_pol(),
		get_comm_atom_seln_pol()
	);

	// Check that the two lists' sizes match
	if ( common_coords.first.size() != common_coords.second.size() ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception("Lists of common coordinates have mismatching numbers of residues"));
	}
	const size_t num_residues = common_coords.first.size();
	for (const size_t &res_ctr : indices( num_residues ) ) {
		if ( common_coords.first[res_ctr].size() != common_coords.second[res_ctr].size() ) {
			BOOST_THROW_EXCEPTION(out_of_range_exception(
				"Lists of common coordinates have mismatching numbers of atoms for common residue "
				+ lexical_cast<string>(res_ctr)
			));
		}
		;
	}


	// Return the list of common coordinates
	return common_coords;
}

/// \brief Less-than operator for score_common_coord_handler
///
/// \relates score_common_coord_handler
score_common_coord_handler_vec cath::score::detail::get_all_score_common_coord_handlers() {
	score_common_coord_handler_vec score_common_coord_handlers;
	for (const common_residue_selection_policy &res_pol : get_all_common_residue_selection_policies() ) {
		for (const common_atom_selection_policy &atom_pol : get_all_common_atom_selection_policies() ) {
			score_common_coord_handlers.emplace_back( res_pol, atom_pol );
		}
	}
	return score_common_coord_handlers;
}

/// \brief Less-than operator for score_common_coord_handler
///
/// \relates score_common_coord_handler
bool cath::score::detail::operator<(const score_common_coord_handler &prm_score_common_coord_handler_a, ///< The first  score_common_coord_handler to compare
                                    const score_common_coord_handler &prm_score_common_coord_handler_b  ///< The second score_common_coord_handler to compare
                                    ) {
	auto the_helper = make_less_than_helper( prm_score_common_coord_handler_a, prm_score_common_coord_handler_b );
	the_helper.register_comparison_field( &score_common_coord_handler::get_comm_res_seln_pol  );
	the_helper.register_comparison_field( &score_common_coord_handler::get_comm_atom_seln_pol );
	return final_less_than_result( the_helper );
}

/// \brief Simple insertion operator for score_common_coord_handler
///
/// \relates score_common_coord_handler
ostream & cath::score::detail::operator<<(ostream                          &prm_os,                        ///< The ostream to which the score_common_coord_handler should be output
                                          const score_common_coord_handler &prm_score_common_coord_handler ///< The score_common_coord_handler to output
                                          ) {
	prm_os << ( "score_common_coord_handler[" + prm_score_common_coord_handler.short_suffix_string() + "]" );
	return prm_os;
}

