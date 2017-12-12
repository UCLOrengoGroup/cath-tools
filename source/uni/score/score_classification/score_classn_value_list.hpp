/// \file
/// \brief The score_classn_value_list class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_SCORE_SCORE_CLASSIFICATION_SCORE_CLASSN_VALUE_LIST_H
#define _CATH_TOOLS_SOURCE_UNI_SCORE_SCORE_CLASSIFICATION_SCORE_CLASSN_VALUE_LIST_H

#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>

#include "common/boost_addenda/string_algorithm/split_build.hpp"
#include "common/exception/runtime_error_exception.hpp"
#include "common/file/open_fstream.hpp"
#include "common/type_aliases.hpp"
#include "score/score_classification/score_classn_value.hpp"
#include "score/score_classification/score_classn_value_better_value.hpp"
#include "score/score_type_aliases.hpp"

#include <fstream>

namespace cath { namespace score { class classn_stat; } }
namespace cath { namespace score { class score_classn_value_list; } }
namespace cath { namespace score { class named_true_false_pos_neg_list; } }

namespace cath {
	namespace score {
		score_classn_value_list make_score_classn_value_list(const score_classn_value_vec &,
		                                                     const bool &,
		                                                     const std::string &);

		/// \brief TODOCUMENT
		///
		/// Invariants:
		///  * score_classn_values is kept sorted at all times, from best to worst
		class score_classn_value_list final {
		private:
			friend score_classn_value_list cath::score::make_score_classn_value_list(const score_classn_value_vec &,
			                                                                         const bool &,
			                                                                         const std::string &);

			/// \brief A vector of score_classn_value objects
			///
			/// (ie scores for instances, each with an label and a flag for whether it's positive)
			score_classn_value_vec score_classn_values;

			/// \brief A functor that defines which of two score_classn_value objects is better
			score_classn_value_better_value better_than;

			/// \brief Name for this series of data
			///
			/// This will often be the algorithm or scoring-scheme that generated the scores
			std::string name;

			void sort_values();

			score_classn_value_list(score_classn_value_vec,
			                        const bool &,
			                        std::string);

		public:
			using const_iterator = score_classn_value_vec_citr;
			using iterator       = score_classn_value_vec_citr;

			bool empty() const;
			size_t size() const;

			const score_classn_value & operator[](const size_t &) const;

			const_iterator begin() const;
			const_iterator end() const;

			const std::string get_name() const;

			const score_classn_value_better_value & get_better_than() const;

			void add_score_classn_value(const score_classn_value &);
		};

//		void fill_missing_instance_labels_with_score(score_classn_value_list &);

		double best_score(const score_classn_value_list &);
		double worst_score(const score_classn_value_list &);

		value_list_scaling get_scaling(const score_classn_value_list &);

		const score_classn_value & best_scoring_actual_positive(const score_classn_value_list &);
		const score_classn_value & best_scoring_actual_negative(const score_classn_value_list &);
		const score_classn_value & worst_scoring_actual_positive(const score_classn_value_list &);
		const score_classn_value & worst_scoring_actual_negative(const score_classn_value_list &);

		std::ostream & summarise_score_classn_value_list(std::ostream &,
		                                                 const score_classn_value_list &);

		str_set get_sorted_instance_labels(const score_classn_value_list &);
		str_vec get_instance_labels(const score_classn_value_list &);
		str_size_pair_vec get_sorted_instance_labels_and_indices(const score_classn_value_list &);

		score_classn_value_vec get_score_classn_values_of_instance_labels(const score_classn_value_list &,
		                                                                  const str_vec &);

		bool instance_labels_match(const score_classn_value_list &,
		                           const score_classn_value_list &);

		const bool & get_higher_is_better(const score_classn_value_list &);

		double worst_possible_score(const score_classn_value_list &);

		score_classn_value_list make_score_classn_value_list(const score_classn_value_vec &,
		                                                     const bool &,
		                                                     const std::string &);

		doub_doub_pair_vec correlated_data(const score_classn_value_list &,
		                                   const score_classn_value_list &);

		/// \brief TODOCUMENT
		///
		/// \relates score_classn_value_list
		template <typename FN>
		score_classn_value_list read_score_classn_value_list(const boost::filesystem::path &arg_path,             ///< TODOCUMENT
		                                                     const bool                    &arg_higher_is_better, ///< TODOCUMENT
		                                                     const std::string             &arg_name,             ///< TODOCUMENT
		                                                     FN                             arg_is_positive_fn    ///< TODOCUMENT
		                                                     ) {
			std::ifstream input_stream;
			common::open_ifstream( input_stream, arg_path );

			score_classn_value_vec score_classn_values;
			std::string line_string;

			while ( std::getline( input_stream, line_string ) ) {
				// Use an istringstream rather than lexical_cast to allow setting boolalpha
				const auto line_parts = common::split_build<str_vec>( line_string, boost::algorithm::is_space(), boost::algorithm::token_compress_on );
				if ( line_parts.size() < 2 ) {
					BOOST_THROW_EXCEPTION(common::runtime_error_exception("Unable to read score_classn_value_list entry from line that doesn't have at least two parts"));
				}
				// This used to
				const std::string &id1         = line_parts[ 0 ];
				const std::string &id2         = line_parts[ 1 ];
				const double       score       = std::stod( line_parts[ 2 ] );
				const bool         is_positive = arg_is_positive_fn( line_parts );

				score_classn_values.emplace_back( score, is_positive, id1 + " " + id2 );

//				std::istringstream line_ss( line_string );
//				line_ss >> std::boolalpha;
//				T t;
//				line_ss >> t;
//				line_entries.push_back( t );
			}
			input_stream.close();

			return make_score_classn_value_list( score_classn_values, arg_higher_is_better, arg_name );
		}

		score_classn_value_list read_svmlight_predictions_file(const boost::filesystem::path &,
		                                                       const std::string &);

		score_classn_value_list_vec read_svmlight_predictions_files(const std::vector<std::pair<boost::filesystem::path, std::string>> &);

		named_true_false_pos_neg_list make_named_true_false_pos_neg_list(const score_classn_value_list &);

		double area_under_curve(const score_classn_value_list &,
		                        const classn_stat &,
		                        const classn_stat &);

		double area_under_roc_curve(const score_classn_value_list &);

	} // namespace score
} // namespace cath

#endif
