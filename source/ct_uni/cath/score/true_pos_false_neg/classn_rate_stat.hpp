/// \file
/// \brief The classn_rate_stat class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_TRUE_POS_FALSE_NEG_CLASSN_RATE_STAT_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_TRUE_POS_FALSE_NEG_CLASSN_RATE_STAT_HPP

#include <boost/fusion/container/generation/map_tie.hpp>
#include <boost/fusion/container/generation/map_tie.hpp>
#include <boost/fusion/include/map_tie.hpp>

#include "cath/common/algorithm/constexpr_find.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/score/true_pos_false_neg/classn_outcome.hpp"
#include "cath/score/true_pos_false_neg/classn_stat.hpp"
#include "cath/score/true_pos_false_neg/true_false_pos_neg.hpp"

using namespace ::cath::common::literals;

namespace cath {
	namespace score {

		/// \brief The standard types of classn_rate_stat that are available
		enum class std_classn_rate_stat : char {
			SENSITIVITY,
			RECALL,
			HIT_RATE,
			TRUE_POSITIVE_RATE,
			SPECIFICITY,
			TRUE_NEGATIVE_RATE,
			PRECISION,
			POSITIVE_PREDICTIVE_VALUE,
			FALL_OUT,
			FALSE_POSITIVE_RATE,
			FALSE_DISCOVERY_RATE
		};

		namespace detail {

			/// \brief Type alias for the tuple that provides the numerator and denominator for each std_classn_rate_stat
			using rate_stat_outcome_outcome_tuple = std::tuple<std_classn_rate_stat, classn_outcome, classn_outcome>;

			/// \brief Store the properties associated with each std_classn_rate_stat for use by classn_rate_stat
			struct properties_of_classn_rate_stat final {

				/// \brief Return a map from std_classn_rate_stat to the corresponding name string
				std::map<std_classn_rate_stat, std::string> operator()() const {
					return {
						{ std_classn_rate_stat::SENSITIVITY,               "Sensitivity"               },
						{ std_classn_rate_stat::RECALL,                    "Recall"                    },
						{ std_classn_rate_stat::HIT_RATE,                  "Hit Rate"                  },
						{ std_classn_rate_stat::TRUE_POSITIVE_RATE,        "True Positive Rate"        },
						{ std_classn_rate_stat::SPECIFICITY,               "Specificity"               },
						{ std_classn_rate_stat::TRUE_NEGATIVE_RATE,        "True Negative Rate"        },
						{ std_classn_rate_stat::PRECISION,                 "Precision"                 },
						{ std_classn_rate_stat::POSITIVE_PREDICTIVE_VALUE, "Positive Predictive Value" },
						{ std_classn_rate_stat::FALL_OUT,                  "Fall Out"                  },
						{ std_classn_rate_stat::FALSE_POSITIVE_RATE,       "False Positive Rate"       },
						{ std_classn_rate_stat::FALSE_DISCOVERY_RATE,      "False Discovery Rate"      }
					};
				}

				/// \brief Store the numerator and denominator for each std_classn_rate_stat
				///
				/// This is implemented as a std::array which allows constexpr lookups in classn_rate_stat,
				/// which also means that keys missing from here will be flagged up at compile time.
				///
				/// This is used in the tests to ensure that the tests are performed on every value key
				static constexpr std::array<rate_stat_outcome_outcome_tuple, 11u> numerator_and_denominator_of_stat{{
					std::make_tuple( std_classn_rate_stat::SENSITIVITY,               classn_outcome::TRUE_POSITIVE,   classn_outcome::FALSE_NEGATIVE ),
					std::make_tuple( std_classn_rate_stat::RECALL,                    classn_outcome::TRUE_POSITIVE,   classn_outcome::FALSE_NEGATIVE ),
					std::make_tuple( std_classn_rate_stat::HIT_RATE,                  classn_outcome::TRUE_POSITIVE,   classn_outcome::FALSE_NEGATIVE ),
					std::make_tuple( std_classn_rate_stat::TRUE_POSITIVE_RATE,        classn_outcome::TRUE_POSITIVE,   classn_outcome::FALSE_NEGATIVE ),

					std::make_tuple( std_classn_rate_stat::SPECIFICITY,               classn_outcome::TRUE_NEGATIVE,   classn_outcome::FALSE_POSITIVE ),
					std::make_tuple( std_classn_rate_stat::TRUE_NEGATIVE_RATE,        classn_outcome::TRUE_NEGATIVE,   classn_outcome::FALSE_POSITIVE ),

					std::make_tuple( std_classn_rate_stat::PRECISION,                 classn_outcome::TRUE_POSITIVE,   classn_outcome::FALSE_POSITIVE ),
					std::make_tuple( std_classn_rate_stat::POSITIVE_PREDICTIVE_VALUE, classn_outcome::TRUE_POSITIVE,   classn_outcome::FALSE_POSITIVE ),

					std::make_tuple( std_classn_rate_stat::FALL_OUT,                  classn_outcome::FALSE_POSITIVE,  classn_outcome::TRUE_NEGATIVE  ),
					std::make_tuple( std_classn_rate_stat::FALSE_POSITIVE_RATE,       classn_outcome::FALSE_POSITIVE,  classn_outcome::TRUE_NEGATIVE  ),

					std::make_tuple( std_classn_rate_stat::FALSE_DISCOVERY_RATE,      classn_outcome::FALSE_POSITIVE,  classn_outcome::TRUE_POSITIVE  )
				}};

				static constexpr size_t num_std_classn_rate_stats = std::tuple_size<
					decltype( properties_of_classn_rate_stat::numerator_and_denominator_of_stat )
				>::value;
			};
		} // namespace detail

		/// \brief Concrete classn_stat for statistics of the form A/(A+B),
		///        (where A and B are one of {TP, TN, FP, FN})
		///
		/// This allows many rate statistics to be simply defined as type aliases.
		/// (eg: sensitivity, recall, hit_rate, true_positive_rate, specificity,
		///      true_negative_rate, precision, positive_predictive_value, fall_out,
		///      false_positive_rate and false_discovery_rate).
		template <std_classn_rate_stat S>
		class classn_rate_stat final : public classn_stat {
		private:
			size_rational do_calculate(const true_false_pos_neg &) const final;

			std::string do_get_name() const final;

			/// \brief The numerator (looked up from numerator_and_denominator_of_stat at compile-time)
			static constexpr classn_outcome numerator   = std::get<1>( common::constexpr_find(
				detail::properties_of_classn_rate_stat::numerator_and_denominator_of_stat,
				S
			) );

			/// \brief The denominator (looked up from numerator_and_denominator_of_stat at compile-time)
			static constexpr classn_outcome denominator = std::get<2>( common::constexpr_find(
				detail::properties_of_classn_rate_stat::numerator_and_denominator_of_stat,
				S
			) );
		};

		/// \brief Calculate the numerator and denominator of the rate and return a rational
		template <std_classn_rate_stat S>
		size_rational classn_rate_stat<S>::do_calculate(const true_false_pos_neg &prm_true_false_pos_neg ///< The true_false_pos_neg from which to calculate the rate
		                                                ) const {
			const size_t &numerator_val = prm_true_false_pos_neg.get_num< numerator   >();
			const size_t &other_val     = prm_true_false_pos_neg.get_num< denominator >();
			return {
				numerator_val,
				std::max( numerator_val + other_val, 1_z )
			};
		}

		/// \brief Calculate the numerator and denominator of the rate and return a rational
		template <std_classn_rate_stat S>
		std::string classn_rate_stat<S>::do_get_name() const {
			return detail::properties_of_classn_rate_stat()().at( S );
		}

		using sensitivity               = classn_rate_stat< std_classn_rate_stat::SENSITIVITY                >;
		using recall                    = classn_rate_stat< std_classn_rate_stat::RECALL                     >;
		using hit_rate                  = classn_rate_stat< std_classn_rate_stat::HIT_RATE                   >;
		using true_positive_rate        = classn_rate_stat< std_classn_rate_stat::TRUE_POSITIVE_RATE         >;
		using specificity               = classn_rate_stat< std_classn_rate_stat::SPECIFICITY                >;
		using true_negative_rate        = classn_rate_stat< std_classn_rate_stat::TRUE_NEGATIVE_RATE         >;
		using precision                 = classn_rate_stat< std_classn_rate_stat::PRECISION                  >;
		using positive_predictive_value = classn_rate_stat< std_classn_rate_stat::POSITIVE_PREDICTIVE_VALUE  >;
		using fall_out                  = classn_rate_stat< std_classn_rate_stat::FALL_OUT                   >;
		using false_positive_rate       = classn_rate_stat< std_classn_rate_stat::FALSE_POSITIVE_RATE        >;
		using false_discovery_rate      = classn_rate_stat< std_classn_rate_stat::FALSE_DISCOVERY_RATE       >;

		/// \brief The rates for the x and y axes respectively of a ROC curve
		using roc_rates              = std::pair<false_positive_rate, true_positive_rate>;

		/// \brief The rates for the x and y axes respectively of a precision-recall graph
		using precision_recall_rates = std::pair<recall,              precision         >;
	} // namespace score
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_TRUE_POS_FALSE_NEG_CLASSN_RATE_STAT_HPP
