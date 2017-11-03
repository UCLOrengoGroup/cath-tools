/// \file
/// \brief The classn_stat_plotter_spec class definitions

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

#include "classn_stat_plotter_spec.hpp"

#include <boost/optional.hpp>

#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "score/true_pos_false_neg/classn_stat_pair_series.hpp"
#include "score/true_pos_false_neg/classn_stat_pair_series_list.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::score;
using namespace std;

using boost::none;
using boost::optional;

/// \brief TODOCUMENT
classn_stat_plotter_spec::classn_stat_plotter_spec(str_vec                                       arg_pre_plot_strs,              ///< TODOCUMENT
                                                   std::vector<std::pair<std::string, str_opt>>  arg_series_to_plot,             ///< TODOCUMENT
                                                   const bool                                   &arg_tidy_up_score_based_legends ///< TODOCUMENT
                                                   ) : pre_plot_strs              { std::move( arg_pre_plot_strs               ) },
                                                       series_to_plot             { std::move( arg_series_to_plot              ) },
                                                       tidy_up_score_based_legends{ arg_tidy_up_score_based_legends              } {
}

/// \brief TODOCUMENT
const str_vec & classn_stat_plotter_spec::get_pre_plot_strs() const {
	return pre_plot_strs;
}

/// \brief TODOCUMENT
const vector<pair<string, str_opt>> & classn_stat_plotter_spec::get_series_to_plot() const {
	return series_to_plot;
}

/// \brief TODOCUMENT
const bool & classn_stat_plotter_spec::get_tidy_up_score_based_legends() const {
	return tidy_up_score_based_legends;
}

/// \brief TODOCUMENT
vector<pair<string, str_opt>> cath::score::get_series_to_plot_or_make_default(const classn_stat_plotter_spec     &arg_spec, ///< TODOCUMENT
                                                                              const classn_stat_pair_series_list &arg_list  ///< TODOCUMENT
                                                                              ) {
	const auto &the_series_to_plot = arg_spec.get_series_to_plot();
	if ( ! the_series_to_plot.empty() ) {
		return the_series_to_plot;
	}

	return transform_build<vector<pair<string, str_opt>>>(
		indices( arg_list.size() ),
		[&] (const size_t &x) {
			return make_pair( arg_list[ x].get_name(), none );
		}
	);
}

/// \brief TODOCUMENT
classn_stat_plotter_spec cath::score::make_standard_score_roc_plotter_spec(const vector<pair<string, str_opt>> &arg_series_to_plot ///< TODOCUMENT
                                                                           ) {
	return {
		{
			"1-x with lines               lc rgb \"#EEEEEE\" notitle",
			"x   with filledcurves y1 = 0 lc rgb \"#EEEEEE\" notitle"
		},
		arg_series_to_plot,
		true
	};
}

/// \brief TODOCUMENT
classn_stat_plotter_spec cath::score::make_standard_score_precision_recall_plotter_spec(const vector<pair<string, str_opt>> &arg_series_to_plot ///< TODOCUMENT
                                                                                        ) {
	return {
		{
			"0.5 with filledcurves y1 = 0 lc rgb \"#EEEEEE\" notitle"
		},
		arg_series_to_plot,
		true
	};
}
