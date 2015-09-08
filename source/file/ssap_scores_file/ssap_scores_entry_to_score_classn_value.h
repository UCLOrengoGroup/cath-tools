/// \file
/// \brief The ssap_scores_entry_to_score_classn_value header

#ifndef SSAP_SCORES_ENTRY_TO_SCORE_CLASSN_VALUE_H_INCLUDED
#define SSAP_SCORES_ENTRY_TO_SCORE_CLASSN_VALUE_H_INCLUDED

#include "common/algorithm/transform_build.h"
#include "file/file_type_aliases.h"
#include "score/score_classification/label_pair_is_positive/label_pair_is_positive.h"
#include "score/score_classification/score_classn_value.h"
#include "score/score_classification/score_classn_value_list.h"

namespace cath {
	namespace file {
		namespace detail {

			/// \brief Helper class for building score_classn_value entries from ssap_scores_entry objects using the
			///        contained ssap_scores_entry getter and label_pair_is_positive to determine which pairs are positive
			class score_classn_value_maker final {
			private:
				/// \brief The ssap_scores_entry getter member function with which to extract the score value
				std::function<double(const ssap_scores_entry&)> f_getter;

				/// \brief A label_pair_is_positive to determine which pairs are positive
				score::label_pair_is_positive the_is_positive;

			public:
				/// \brief Ctor from the ssap_scores_entry getter and label_pair_is_positive
				template <typename GETR_FN>
				score_classn_value_maker(const GETR_FN                       &arg_getter,     ///< The ssap_scores_entry getter member function with which to extract the score value
				                         const score::label_pair_is_positive &arg_is_positive ///< A label_pair_is_positive to determine which pairs are positive
				                         ) : f_getter       ( arg_getter      ),
				                             the_is_positive( arg_is_positive ) {
				}

				/// \brief Build a score_classn_value from the specified ssap_scores_entry using the stored field getter and label_pair_is_positive
				score::score_classn_value operator()(const ssap_scores_entry &arg_entry ///< The ssap_scores_entry containing the data from which the score_classn_value should be built
				                                     ) {
					const std::string &name_1 = arg_entry.get_name_1();
					const std::string &name_2 = arg_entry.get_name_2();
					return {
						f_getter( arg_entry ),
						the_is_positive.is_positive( name_1, name_2 ),
						name_1 + " " + name_2
					};
				}
			};

		}

		/// \brief Build a score_classn_value_list for the specified ssap_scores_entry field from the specified ssap_scores_entry_vec
		template <typename GETR_FN>
		score::score_classn_value_list make_val_list_of_ssap_scores_entries(const ssap_scores_entry_vec         &arg_entries,          ///< The ssap_scores_entry_vec data from which to generate the score_classn_value_list
		                                                                    const score::label_pair_is_positive &arg_positive_fn,      ///< The label_pair_is_positive to specify which of the ssap_scores_entries corresponds to a positive pair
		                                                                    GETR_FN                              arg_getter_fn,        ///< The ssap_scores_entry getter for the required field
		                                                                    const bool                          &arg_higher_is_better, ///< Whether higher is better for this score
		                                                                    const std::string                   &arg_name              ///< The name to use for this field
		                                                                    ) {
			return make_score_classn_value_list(
				common::transform_build<score::score_classn_value_vec>(
					arg_entries,
					detail::score_classn_value_maker( arg_getter_fn, arg_positive_fn )
				),
				arg_higher_is_better,
				arg_name
			);
		}

	}
}

#endif
