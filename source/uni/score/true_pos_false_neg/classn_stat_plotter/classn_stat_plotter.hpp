/// \file
/// \brief The classn_stat_plotter class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_SCORE_TRUE_POS_FALSE_NEG_CLASSN_STAT_PLOTTER_CLASSN_STAT_PLOTTER_H
#define _CATH_TOOLS_SOURCE_UNI_SCORE_TRUE_POS_FALSE_NEG_CLASSN_STAT_PLOTTER_CLASSN_STAT_PLOTTER_H

#include <boost/filesystem/path.hpp>

#include "common/type_aliases.hpp"
#include "score/score_type_aliases.hpp"

namespace cath { namespace score { class classn_stat; } }
namespace cath { namespace score { class classn_stat_pair_series; } }
namespace cath { namespace score { class classn_stat_pair_series_list; } }
namespace cath { namespace score { class classn_stat_plotter_spec; } }
namespace cath { namespace score { class named_true_false_pos_neg_list; } }
namespace cath { namespace score { class score_classn_value_results_set; } }

namespace cath {
	namespace score {

		/// \brief TODOCUMENT
		class classn_stat_plotter final {
		private:
			/// \brief TODOCUMENT
			bool suppress_execution = false;

			std::string process_legend_name_copy(const std::string &,
			                                     const bool &) const;

		public:
			classn_stat_plotter() = default;
			explicit classn_stat_plotter(const bool &);

			void plot(const boost::filesystem::path &,
			          const score_classn_value_list_vec &,
			          const classn_stat &,
			          const classn_stat &,
			          const classn_stat_plotter_spec &) const;

			void plot(const boost::filesystem::path &,
			          const score_classn_value_results_set &,
			          const classn_stat &,
			          const classn_stat &,
			          const classn_stat_plotter_spec &) const;

			void plot(const boost::filesystem::path &,
			          const named_true_false_pos_neg_list &,
			          const classn_stat &,
			          const classn_stat &,
			          const classn_stat_plotter_spec &) const;

			void plot(const boost::filesystem::path &,
			          const classn_stat_pair_series &,
			          const std::string &,
			          const std::string &,
			          const classn_stat_plotter_spec &) const;

			void plot(const boost::filesystem::path &,
			          const classn_stat_pair_series_list &,
			          const std::string &,
			          const std::string &,
			          const classn_stat_plotter_spec &) const;
		};

		template <typename T>
		void plot_classn_stat(const classn_stat_plotter         &arg_classn_stat_plotter,      ///< TODOCUMENT
		                      const boost::filesystem::path     &arg_output_file_stem,         ///< TODOCUMENT
		                      const score_classn_value_list_vec &arg_score_classn_value_lists, ///< TODOCUMENT
		                      const classn_stat_plotter_spec    &arg_plot_spec                 ///< TODOCUMENT
		                      ) {
			using first_classn_stat  = typename T::first_type;
			using second_classn_stat = typename T::second_type;
			arg_classn_stat_plotter.plot(
				arg_output_file_stem,
				arg_score_classn_value_lists,
				first_classn_stat(),
				second_classn_stat(),
				arg_plot_spec
			);
		}

		template <typename T>
		void plot_classn_stat(const classn_stat_plotter            &arg_classn_stat_plotter,      ///< TODOCUMENT
		                      const boost::filesystem::path        &arg_output_file_stem,         ///< TODOCUMENT
		                      const score_classn_value_results_set &arg_scv_results_set,          ///< TODOCUMENT
		                      const classn_stat_plotter_spec       &arg_plot_spec                 ///< TODOCUMENT
		                      ) {
			using first_classn_stat  = typename T::first_type;
			using second_classn_stat = typename T::second_type;
			arg_classn_stat_plotter.plot(
				arg_output_file_stem,
				arg_scv_results_set,
				first_classn_stat(),
				second_classn_stat(),
				arg_plot_spec
			);
		}

		template <typename T>
		void plot_classn_stat(const classn_stat_plotter           &arg_classn_stat_plotter,      ///< TODOCUMENT
		                      const boost::filesystem::path       &arg_output_file_stem,         ///< TODOCUMENT
		                      const named_true_false_pos_neg_list &arg_named_tfpn,               ///< TODOCUMENT
		                      const classn_stat_plotter_spec      &arg_plot_spec                 ///< TODOCUMENT
		                      ) {
			using first_classn_stat  = typename T::first_type;
			using second_classn_stat = typename T::second_type;
			arg_classn_stat_plotter.plot(
				arg_output_file_stem,
				arg_named_tfpn,
				first_classn_stat(),
				second_classn_stat(),
				arg_plot_spec
			);
		}

		void plot_roc(const classn_stat_plotter &,
		              const boost::filesystem::path &,
		              const score_classn_value_list_vec &,
		              const classn_stat_plotter_spec &);

		void plot_roc(const classn_stat_plotter &,
		              const boost::filesystem::path &,
		              const score_classn_value_results_set &,
		              const classn_stat_plotter_spec &);

		void plot_roc(const classn_stat_plotter &,
		              const boost::filesystem::path &,
		              const named_true_false_pos_neg_list &,
		              const classn_stat_plotter_spec &);

		void plot_precision_recall(const classn_stat_plotter &,
		                           const boost::filesystem::path &,
		                           const score_classn_value_list_vec &,
		                           const classn_stat_plotter_spec &);

		void plot_precision_recall(const classn_stat_plotter &,
		                           const boost::filesystem::path &,
		                           const score_classn_value_results_set &,
		                           const classn_stat_plotter_spec &);

		void plot_precision_recall(const classn_stat_plotter &,
		                           const boost::filesystem::path &,
		                           const named_true_false_pos_neg_list &,
		                           const classn_stat_plotter_spec &);

	} // namespace score
} // namespace cath

#endif
