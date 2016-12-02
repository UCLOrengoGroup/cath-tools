/// \file
/// \brief The score_classn_value_results_set class header

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

#ifndef _CATH_TOOLS_SOURCE_SCORE_SCORE_CLASSIFICATION_SCORE_CLASSN_VALUE_RESULTS_SET_H
#define _CATH_TOOLS_SOURCE_SCORE_SCORE_CLASSIFICATION_SCORE_CLASSN_VALUE_RESULTS_SET_H

#include <boost/filesystem.hpp>

#include "exception/invalid_argument_exception.hpp"
#include "score/aligned_pair_score_list/aligned_pair_score_value_list.hpp"
#include "score/aligned_pair_score_list/score_value_list_reader/score_value_reader.hpp"
#include "score/score_type_aliases.hpp"
#include "score/score_classification/score_classn_value_list.hpp"
#include "score/true_pos_false_neg/classn_stat_pair_series_list.hpp"

#include <iostream>
#include <random>

namespace cath { namespace score { class aligned_pair_score_value_list; } }
namespace cath { namespace score { class classn_stat; } }
namespace cath { namespace score { class classn_stat_pair_series_list; } }
namespace cath { namespace score { class named_true_false_pos_neg_list_list; } }

namespace cath {
	namespace score {

		/// \brief TODOCUMENT
		///
		/// This allows data to be added from:
		///  * aligned_pair_score_value_lists (ie different types of score for one instance) and
		///  * score_classn_value_lists (ie scores for different instances for one type of algorithm)
		///
		/// Note that if adding both these types of data, one should add all the aligned_pair_score_value_lists
		/// first and then add any score_classn_value_lists (because that way, the data is a complete matrix
		/// at all times whereas once a score_classn_value_list has been added with all the instances,
		/// adding a further aligned_pair_score_value_list will likely be rejected for creating gaps).
		/// (Why the asymmetry? Because the score_classn_value_lists typically cover all instances whereas the
		///  the aligned_pair_score_value_lists typically don't)
		///
		/// This is currently implemented as a sorted score_classn_value_list_vec but this
		/// means that instance_label strings get repeated in each of the equivalent entries in each
		/// of the score_classn_value_lists. If this turns out to use too much memory, it could be
		/// changed to a design that stores:
		///  * one list of score_classn_value_list names
		///  * a list of instance labels
		///  * data indexed by the labels in the two lists above
		///
		/// Invariant:
		///  * score_classn_value_lists is kept sorted and uniqued on the names
		///
		/// It makes sense to keep the data uniqued on name because the names are used
		/// as indices with which to add data.
		///
		/// If adding data to a score_classn_value_results_set from aligned_pair_score_value_lists (ie different
		/// types of score for one instance) *and* from score_classn_value_lists
		class score_classn_value_results_set final {
		private:
			/// \brief TODOCUMENT
			score_classn_value_list_vec score_classn_value_lists;

			bool is_sorted_uniqued() const;
			void check_is_sorted_uniqued() const;
			void sort_score_classn_value_lists();

			const score_classn_value_list & get_score_classn_value_list_of_name_impl(const std::string &) const;
			score_classn_value_list & get_score_classn_value_list_of_name(const std::string &);

		public:
			using iterator       = score_classn_value_list_vec_citr;
			using const_iterator = score_classn_value_list_vec_citr;

			const score_classn_value_list & get_score_classn_value_list_of_name(const std::string &) const;

			bool empty() const;
			size_t size() const;

			const score_classn_value_list & operator[](const size_t &) const;

			const_iterator begin() const;
			const_iterator end() const;

			void add_score_classn_value_list(const score_classn_value_list &);
			void add_aligned_pair_score_value_list(const aligned_pair_score_value_list &,
			                                       const bool &,
			                                       const std::string &);
		};

		score_classn_value_results_set make_score_classn_value_results_set(const score_classn_value_list_vec &);

		score_classn_value_list_vec make_score_classn_value_list_vec(const score_classn_value_results_set &);

		value_list_scaling_vec get_value_list_scalings(const score_classn_value_results_set &);

		doub_doub_pair_vec get_correlated_data(const score_classn_value_results_set &,
		                                       const std::string &,
		                                       const std::string &);

		namespace detail {
			void write_to_svm_light_data_files_impl(const score_classn_value_vec_vec &,
			                                        const value_list_scaling_vec &,
			                                        const boost::filesystem::path &,
			                                        const size_vec &);
		} // namespace detail

		void write_to_svm_light_data_files(const score_classn_value_results_set &,
		                                   const boost::filesystem::path &,
		                                   const size_t &,
		                                   std::mt19937 &,
		                                   const double & = 0.5);

		void add_score_classn_value_list_and_add_missing(score_classn_value_results_set &,
		                                                 score_classn_value_list,
		                                                 const double &,
		                                                 const bool & = false);

		str_set get_sorted_instance_labels(const score_classn_value_results_set &);

		str_vec get_instance_labels(const score_classn_value_results_set &);

		size_t get_num_instances(const score_classn_value_results_set &);

		std::string get_name_of_index(const score_classn_value_results_set &,
		                              const size_t &);

		str_vec get_names(const score_classn_value_results_set &);

		const bool & get_higher_is_better_of_index(const score_classn_value_results_set &,
		                                           const size_t &);

		const score_classn_value_list & find_score_classn_value_list_of_name(const score_classn_value_results_set &,
		                                                                     const std::string &);

		named_true_false_pos_neg_list_list make_named_true_false_pos_neg_list_list(const score_classn_value_results_set &);

		named_true_false_pos_neg_list_list make_named_true_false_pos_neg_list_list(const score_classn_value_list_vec &);

		classn_stat_pair_series_list make_classn_stat_pair_series_list(const score_classn_value_results_set &,
		                                                               const classn_stat &,
		                                                               const classn_stat &);

		classn_stat_pair_series_list make_classn_stat_pair_series_list(const score_classn_value_list_vec &,
		                                                               const classn_stat &,
		                                                               const classn_stat &);

		classn_stat_pair_series_list make_roc_series_list(const score_classn_value_results_set &);
		classn_stat_pair_series_list make_precision_recall_series_list(const score_classn_value_results_set &);

		template <typename T>
		classn_stat_pair_series_list make_standard_classn_stat_pair_series_list(const score_classn_value_results_set &arg_score_classn_value_results_set ///< TODOCUMENT
		                                                                        ) {
			using first_classn_stat  = typename T::first_type;
			using second_classn_stat = typename T::second_type;
			return make_classn_stat_pair_series_list(
				arg_score_classn_value_results_set,
				first_classn_stat(),
				second_classn_stat()
			);
		}


		/// \brief TODOCUMENT
		///
		/// FN should take a single boost::filesystem::path and return optional<pair<bool, string>>,
		/// representing files that should be processed with a bool indicating whether it's a positive instance and a label
		template <typename FN>
		score_classn_value_results_set read_from_dir(const boost::filesystem::path &arg_directory,                     ///< TODOCUMENT
		                                             const FN                       arg_positive_and_label_of_filename ///< TODOCUMENT
		                                             ) {
			if ( ! boost::filesystem::is_directory( arg_directory ) ) {
				BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Cannot read_from_dir() for non-directory path"));
			}

			const auto the_dir_range = boost::make_iterator_range(
				boost::filesystem::directory_iterator( arg_directory ),
				boost::filesystem::directory_iterator(               )
			);

			score_classn_value_results_set the_results;
			size_t read_file_ctr = 0;
			size_t file_ctr      = 0;
			for (const auto &dir_entry : the_dir_range) {
				const auto the_filename       = dir_entry.path();
				const auto positive_and_label = arg_positive_and_label_of_filename( the_filename );
				if ( positive_and_label ) {
					std::cerr << "About to attempt to read file " << std::right << std::setw(6) << read_file_ctr << " / " << file_ctr << " : " << the_filename << std::endl;
					const auto the_scores = score_value_reader::read( the_filename );
					the_results.add_aligned_pair_score_value_list(
						the_scores,
						positive_and_label->first,
						positive_and_label->second
					);
					++read_file_ctr;
				}
				++file_ctr;

//				if ( read_file_ctr > 5000 ) {
//				if ( read_file_ctr > 2000 ) {
//				if ( read_file_ctr > 500 ) {
//					break;
//				}
			}
			return the_results;
		}

	} // namespace score
} // namespace cath

#endif
