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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_TRUE_POS_FALSE_NEG_CLASSN_STAT_PLOTTER_CLASSN_STAT_PLOTTER_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_TRUE_POS_FALSE_NEG_CLASSN_STAT_PLOTTER_CLASSN_STAT_PLOTTER_HPP

#include <filesystem>

#include "cath/common/type_aliases.hpp"
#include "cath/score/score_type_aliases.hpp"

// clang-format off
namespace cath::score { class classn_stat; }
namespace cath::score { class classn_stat_pair_series; }
namespace cath::score { class classn_stat_pair_series_list; }
namespace cath::score { class classn_stat_plotter_spec; }
namespace cath::score { class named_true_false_pos_neg_list; }
namespace cath::score { class score_classn_value_results_set; }
// clang-format on

namespace cath::score {

	/// \brief TODOCUMENT
	class classn_stat_plotter final {
	private:
		/// \brief TODOCUMENT
		bool suppress_execution = false;

		[[nodiscard]] std::string process_legend_name_copy( const std::string &, const bool & ) const;

	  public:
		classn_stat_plotter() = default;
		explicit classn_stat_plotter(const bool &);

		void plot(const ::std::filesystem::path &,
		          const score_classn_value_list_vec &,
		          const classn_stat &,
		          const classn_stat &,
		          const classn_stat_plotter_spec &) const;

		void plot(const ::std::filesystem::path &,
		          const score_classn_value_results_set &,
		          const classn_stat &,
		          const classn_stat &,
		          const classn_stat_plotter_spec &) const;

		void plot(const ::std::filesystem::path &,
		          const named_true_false_pos_neg_list &,
		          const classn_stat &,
		          const classn_stat &,
		          const classn_stat_plotter_spec &) const;

		void plot(const ::std::filesystem::path &,
		          const classn_stat_pair_series &,
		          const std::string &,
		          const std::string &,
		          const classn_stat_plotter_spec &) const;

		void plot(const ::std::filesystem::path &,
		          const classn_stat_pair_series_list &,
		          const std::string &,
		          const std::string &,
		          const classn_stat_plotter_spec &) const;
	};

	template <typename T>
	void plot_classn_stat(const classn_stat_plotter         &prm_classn_stat_plotter,      ///< TODOCUMENT
	                      const ::std::filesystem::path     &prm_output_file_stem,         ///< TODOCUMENT
	                      const score_classn_value_list_vec &prm_score_classn_value_lists, ///< TODOCUMENT
	                      const classn_stat_plotter_spec    &prm_plot_spec                 ///< TODOCUMENT
	                      ) {
		using first_classn_stat  = typename T::first_type;
		using second_classn_stat = typename T::second_type;
		prm_classn_stat_plotter.plot(
			prm_output_file_stem,
			prm_score_classn_value_lists,
			first_classn_stat(),
			second_classn_stat(),
			prm_plot_spec
		);
	}

	template <typename T>
	void plot_classn_stat(const classn_stat_plotter            &prm_classn_stat_plotter,      ///< TODOCUMENT
	                      const ::std::filesystem::path        &prm_output_file_stem,         ///< TODOCUMENT
	                      const score_classn_value_results_set &prm_scv_results_set,          ///< TODOCUMENT
	                      const classn_stat_plotter_spec       &prm_plot_spec                 ///< TODOCUMENT
	                      ) {
		using first_classn_stat  = typename T::first_type;
		using second_classn_stat = typename T::second_type;
		prm_classn_stat_plotter.plot(
			prm_output_file_stem,
			prm_scv_results_set,
			first_classn_stat(),
			second_classn_stat(),
			prm_plot_spec
		);
	}

	template <typename T>
	void plot_classn_stat(const classn_stat_plotter           &prm_classn_stat_plotter,      ///< TODOCUMENT
	                      const ::std::filesystem::path       &prm_output_file_stem,         ///< TODOCUMENT
	                      const named_true_false_pos_neg_list &prm_named_tfpn,               ///< TODOCUMENT
	                      const classn_stat_plotter_spec      &prm_plot_spec                 ///< TODOCUMENT
	                      ) {
		using first_classn_stat  = typename T::first_type;
		using second_classn_stat = typename T::second_type;
		prm_classn_stat_plotter.plot(
			prm_output_file_stem,
			prm_named_tfpn,
			first_classn_stat(),
			second_classn_stat(),
			prm_plot_spec
		);
	}

	void plot_roc(const classn_stat_plotter &,
	              const ::std::filesystem::path &,
	              const score_classn_value_list_vec &,
	              const classn_stat_plotter_spec &);

	void plot_roc(const classn_stat_plotter &,
	              const ::std::filesystem::path &,
	              const score_classn_value_results_set &,
	              const classn_stat_plotter_spec &);

	void plot_roc(const classn_stat_plotter &,
	              const ::std::filesystem::path &,
	              const named_true_false_pos_neg_list &,
	              const classn_stat_plotter_spec &);

	void plot_precision_recall(const classn_stat_plotter &,
	                           const ::std::filesystem::path &,
	                           const score_classn_value_list_vec &,
	                           const classn_stat_plotter_spec &);

	void plot_precision_recall(const classn_stat_plotter &,
	                           const ::std::filesystem::path &,
	                           const score_classn_value_results_set &,
	                           const classn_stat_plotter_spec &);

	void plot_precision_recall(const classn_stat_plotter &,
	                           const ::std::filesystem::path &,
	                           const named_true_false_pos_neg_list &,
	                           const classn_stat_plotter_spec &);

} // namespace cath::score

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_TRUE_POS_FALSE_NEG_CLASSN_STAT_PLOTTER_CLASSN_STAT_PLOTTER_HPP
