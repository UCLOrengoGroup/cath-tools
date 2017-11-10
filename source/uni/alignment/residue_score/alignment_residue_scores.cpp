/// \file
/// \brief The alignment_residue_scores class definitions

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

#include "alignment_residue_scores.hpp"

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "alignment/alignment.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/type_aliases.hpp"

#include <iostream> // ***** TEMPORARY *****

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace std;

using boost::algorithm::any_of;
using boost::numeric_cast;

/// \brief Ctor for alignment_residue_scores
alignment_residue_scores::alignment_residue_scores(const size_t     &arg_num_entries,                    ///< TODOCUMENT
                                                   size_vec          arg_num_present_entries_by_index,   ///< TODOCUMENT
                                                   score_opt_vec_vec arg_scores_to_other_present_entries ///< TODOCUMENT
                                                   ) : num_entries                    { arg_num_entries                                  },
                                                       num_present_entries_by_index   { std::move( arg_num_present_entries_by_index    ) },
                                                       scores_to_other_present_entries{ std::move( arg_scores_to_other_present_entries ) } {
	sanity_check();
}

/// \brief TODOCUMENT
void alignment_residue_scores::sanity_check() const {
	for (const size_t &num_present_entries : num_present_entries_by_index) {
		if ( num_present_entries > num_entries ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Number of present entries in alignment_residue_scores exceeds num_entries"));
		}
	}
}

/// \brief TODOCUMENT
size_t alignment_residue_scores::get_num_entries() const {
	return num_entries;
}

/// \brief TODOCUMENT
size_t alignment_residue_scores::get_length() const {
	return num_present_entries_by_index.size();
}

/// \brief TODOCUMENT
size_t alignment_residue_scores::get_num_present_entries_of_index(const size_t &arg_index ///< TODOCUMENT
                                                                  ) const {
	return num_present_entries_by_index[ arg_index ];
}

/// \brief TODOCUMENT
score_opt alignment_residue_scores::get_opt_score_to_other_present_entries(const size_t &arg_entry, ///< TODOCUMENT
                                                                           const size_t &arg_index  ///< TODOCUMENT
                                                                           ) const {
	return scores_to_other_present_entries[ arg_entry ][ arg_index ];
}

/// \brief TODOCUMENT
///
/// \relates alignment_residue_scores
bool cath::align::has_score(const alignment_residue_scores &arg_scores, ///< TODOCUMENT
                            const size_t                   &arg_entry,  ///< TODOCUMENT
                            const size_t                   &arg_index   ///< TODOCUMENT
                            ) {
	return static_cast<bool>( arg_scores.get_opt_score_to_other_present_entries( arg_entry, arg_index ) );
}

/// \brief TODOCUMENT
///
/// \relates alignment_residue_scores
float_score_type cath::align::get_unnormalised_score(const alignment_residue_scores &arg_scores,              ///< TODOCUMENT
                                                     const size_t                   &arg_entry,               ///< TODOCUMENT
                                                     const size_t                   &arg_index,               ///< TODOCUMENT
                                                     const bool                     &arg_to_all_other_entries ///< TODOCUMENT to all other entries rather than just to other present entries (ie penalise for positions where few entries are aligned)
                                                     ) {
	if ( ! has_score( arg_scores, arg_entry, arg_index ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot get absent score from alignment_residue_scores"));
	}
	const size_t num_entries         = arg_scores.get_num_entries();
	const size_t num_present_entries = arg_scores.get_num_present_entries_of_index( arg_index );

//	if ( num_present_entries <= 1 ) {
//		return 0.0;
////		BOOST_THROW_EXCEPTION(invalid_argument_exception(
////			"Cannot get score from alignment_residue_scores at entry "
////			+ lexical_cast<string>( arg_entry )
////			+ " and index "
////			+ lexical_cast<string>( arg_index )
////			+ " because it only has "
////			+ lexical_cast<string>( num_present_entries )
////			+ " present entries"
////		));
//	}
	const score_opt &score_to_other_present_entries = arg_scores.get_opt_score_to_other_present_entries( arg_entry,  arg_index );
	if ( ! arg_to_all_other_entries ) {
		return *score_to_other_present_entries;
	}
	const float_score_type num_other_present_entries = numeric_cast<float_score_type>( num_present_entries - 1 );
	const float_score_type num_all_other_entries     = numeric_cast<float_score_type>( num_entries         - 1 );
	return ( *score_to_other_present_entries * num_other_present_entries / num_all_other_entries );
}

/// \brief TODOCUMENT
///
/// \relates alignment_residue_scores
float_score_type cath::align::get_max_score(const alignment_residue_scores &arg_scores,              ///< TODOCUMENT
                                            const bool                     &arg_to_all_other_entries ///< TODOCUMENT
                                            ) {
	float_score_type max_score = 0.0;
	const size_t num_entries = arg_scores.get_num_entries();
	const size_t length      = arg_scores.get_length();
	for (const size_t &entry : indices( num_entries ) ) {
		for (const size_t &index : indices( length ) ) {
			if ( has_score( arg_scores, entry, index ) ) {
				max_score = max( max_score, get_unnormalised_score( arg_scores, entry, index, arg_to_all_other_entries ) );
			}
		}
	}
	return max_score;
}

/// \brief TODOCUMENT
///
/// \relates alignment_residue_scores
float_score_type cath::align::get_score(const alignment_residue_scores &arg_scores,               ///< TODOCUMENT
                                        const size_t                   &arg_entry,                ///< TODOCUMENT
                                        const size_t                   &arg_index,                ///< TODOCUMENT
                                        const bool                     &arg_to_all_other_entries, ///< TODOCUMENT
                                        const bool                     &arg_normalise             ///< TODOCUMENT
                                        ) {
	const float_score_type raw_score = get_unnormalised_score( arg_scores, arg_entry, arg_index, arg_to_all_other_entries );
	return arg_normalise ? raw_score / get_max_score( arg_scores, arg_to_all_other_entries )
	                     : raw_score;
}

/// \brief TODOCUMENT
///
/// \relates alignment_residue_scores
float_score_type cath::align::get_score_to_present_entries_of_index(const alignment_residue_scores &arg_scores, ///< TODOCUMENT
                                                                    const size_t                   &arg_entry,  ///< TODOCUMENT
                                                                    const size_t                   &arg_index   ///< TODOCUMENT
                                                                    ) {
	return get_unnormalised_score( arg_scores, arg_entry, arg_index, false );
}

/// \brief TODOCUMENT
///
/// \relates alignment_residue_scores
float_score_type cath::align::get_score_to_all_entries(const alignment_residue_scores &arg_scores, ///< TODOCUMENT
                                                       const size_t                   &arg_entry,  ///< TODOCUMENT
                                                       const size_t                   &arg_index   ///< TODOCUMENT
                                                       ) {
	return get_unnormalised_score( arg_scores, arg_entry, arg_index, true );
}

/// \brief TODOCUMENT
///
/// \relates alignment_residue_scores
float_score_type cath::align::get_normalised_score_to_present_entries_of_index(const alignment_residue_scores &arg_scores, ///< TODOCUMENT
                                                                               const size_t                   &arg_entry,  ///< TODOCUMENT
                                                                               const size_t                   &arg_index   ///< TODOCUMENT
                                                                               ) {
	return get_score( arg_scores, arg_entry, arg_index, true, true );
}

/// \brief TODOCUMENT
///
/// \relates alignment_residue_scores
float_score_type cath::align::get_normalised_score_to_all_entries(const alignment_residue_scores &arg_scores, ///< TODOCUMENT
                                                                  const size_t                   &arg_entry,  ///< TODOCUMENT
                                                                  const size_t                   &arg_index   ///< TODOCUMENT
                                                                  ) {
	return get_score( arg_scores, arg_entry, arg_index, false, true );
}

/// \brief TODOCUMENT
///
/// \relates alignment_residue_scores
alignment_residue_scores cath::align::make_alignment_residue_scores(const alignment         &arg_alignment, ///< TODOCUMENT
                                                                    const score_opt_vec_vec &arg_scores     ///< TODOCUMENT
                                                                    ) {
	const size_t num_entries = arg_alignment.num_entries();
	const size_t length      = arg_alignment.length();

	if ( arg_scores.size() != num_entries ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Number of entries in scores doesn't match number of entries in alignment"));
	}
	if ( any_of( arg_scores, [&] (const score_opt_vec &x) { return x.size() != length; } ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Length of scores for entry doesn't match length of alignment"));
	}

	return alignment_residue_scores(
		arg_alignment.num_entries(),
		num_present_positions_by_index( arg_alignment ),
		arg_scores
	);
}



