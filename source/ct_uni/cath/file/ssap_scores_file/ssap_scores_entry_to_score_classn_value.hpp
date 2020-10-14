/// \file
/// \brief The ssap_scores_entry_to_score_classn_value header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_FILE_SSAP_SCORES_FILE_SSAP_SCORES_ENTRY_TO_SCORE_CLASSN_VALUE_HPP
#define _CATH_TOOLS_SOURCE_UNI_FILE_SSAP_SCORES_FILE_SSAP_SCORES_ENTRY_TO_SCORE_CLASSN_VALUE_HPP

#include <boost/tuple/tuple.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "cath/common/algorithm/transform_build.hpp"
#include "cath/file/file_type_aliases.hpp"
#include "cath/score/score_classification/label_pair_is_positive/label_pair_is_positive.hpp"
#include "cath/score/score_classification/score_classn_value.hpp"
#include "cath/score/score_classification/score_classn_value_list.hpp"

namespace cath {
	namespace file {
		namespace detail {

		/// \brief Helper class for building score_classn_value entries from ssap_scores_entry objects using the
		///        contained ssap_scores_entry getter and label_pair_is_positive to determine which pairs are positive
		template <typename GETR_FN, typename NAME_FN>
			class score_classn_value_maker final {
			private:
				/// \brief A label_pair_is_positive to determine which pairs are positive
				score::label_pair_is_positive the_is_positive;

				/// \brief The ssap_scores_entry getter member function with which to extract the score value
				GETR_FN f_score_getter;

				/// \brief The ssap_scores_entry getter member function with which to extract the score value
				NAME_FN f_name_getter;

			public:
				/// \brief Ctor from the ssap_scores_entry getter and label_pair_is_positive
				score_classn_value_maker(score::label_pair_is_positive  prm_is_positive,  ///< A label_pair_is_positive to determine which pairs are positive
				                         const GETR_FN                 &prm_score_getter, ///< The ssap_scores_entry getter member function with which to extract the score value
				                         const NAME_FN                 &prm_name_getter   ///< TODOCUMENT
				                         ) : the_is_positive { std::move( prm_is_positive ) },
				                             f_score_getter  { prm_score_getter             },
				                             f_name_getter   { prm_name_getter              } {
				}

				/// \brief Build a score_classn_value from the specified ssap_scores_entry using the stored field getter and label_pair_is_positive
				template <typename T>
				score::score_classn_value operator()(const T &prm_entry ///< The ssap_scores_entry containing the data from which the score_classn_value should be built
				                                     ) {
					const auto names = f_name_getter( prm_entry );
					const std::string &name_1 = names.first;
					const std::string &name_2 = names.second;
					return {
						boost::numeric_cast<double>( f_score_getter( prm_entry ) ),
						the_is_positive.is_positive( name_1, name_2 ),
						name_1 + " " + name_2
					};
				}
			};
		} // namespace detail

		/// \brief Build a score_classn_value_list for the specified ssap_scores_entry field from the specified ssap_scores_entry_vec
		template <typename T, typename GETR_FN, typename NAME_FN>
		score::score_classn_value_list make_val_list_of_scores_entries(const T                             &prm_entries,          ///< The ssap_scores_entry_vec data from which to generate the score_classn_value_list
		                                                               const score::label_pair_is_positive &prm_positive_fn,      ///< The label_pair_is_positive to specify which of the ssap_scores_entries corresponds to a positive pair
		                                                               const bool                          &prm_higher_is_better, ///< Whether higher is better for this score
		                                                               const std::string                   &prm_name,             ///< The name to use for this field
		                                                               GETR_FN                              prm_score_getter_fn,  ///< The ssap_scores_entry getter for the required field
		                                                               NAME_FN                              prm_name_getter_fn    ///< The getter for the names
		                                                               ) {
			return make_score_classn_value_list(
				common::transform_build<score::score_classn_value_vec>(
					prm_entries,
					detail::score_classn_value_maker<GETR_FN, NAME_FN>( prm_positive_fn, prm_score_getter_fn, prm_name_getter_fn )
				),
				prm_higher_is_better,
				prm_name
			);
		}

	} // namespace file
} // namespace cath

#endif
