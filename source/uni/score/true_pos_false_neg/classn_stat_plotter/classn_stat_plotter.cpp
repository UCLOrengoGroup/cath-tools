/// \file
/// \brief The classn_stat_plotter class definitions

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

#include "classn_stat_plotter.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/join.hpp>

#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/filesystem/replace_extension_copy.hpp"
#include "common/third_party_code/gnuplot-iostream.h"
#include "score/score_classification/score_classn_value_results_set.hpp"
#include "score/true_pos_false_neg/classn_rate_stat.hpp"
#include "score/true_pos_false_neg/classn_stat_pair_series.hpp"
#include "score/true_pos_false_neg/classn_stat_pair_series_list.hpp"
#include "score/true_pos_false_neg/classn_stat_plotter/classn_stat_plotter_spec.hpp"
#include "score/true_pos_false_neg/named_true_false_pos_neg_list.hpp"

using namespace cath::common;
using namespace cath::score;
using namespace std;

using boost::algorithm::icontains;
using boost::algorithm::join;
using boost::algorithm::replace_all_copy;
using boost::filesystem::path;
using boost::lexical_cast;
using boost::range::join;
using gnuplotio::Gnuplot;

/// \brief TODOCUMENT
string classn_stat_plotter::process_legend_name_copy(const string &arg_legend_name, ///< TODOCUMENT
                                                     const bool   &arg_tidy         ///< TODOCUMENT
                                                     ) const {
	if ( ! arg_tidy ) {
		return arg_legend_name;
	}
	return replace_all_copy( arg_legend_name, "_", "\\_" );
}

/// \brief TODOCUMENT
classn_stat_plotter::classn_stat_plotter(const bool &arg_suppress_execution ///< TODOCUMENT
                                         ) : suppress_execution ( arg_suppress_execution ) {
}

/// \brief TODOCUMENT
void classn_stat_plotter::plot(const path                        &arg_output_file_stem,         ///< TODOCUMENT
                               const score_classn_value_list_vec &arg_score_classn_value_lists, ///< TODOCUMENT
                               const classn_stat                 &arg_x_classn_stat,            ///< TODOCUMENT
                               const classn_stat                 &arg_y_classn_stat,            ///< TODOCUMENT
                               const classn_stat_plotter_spec    &arg_plot_spec                 ///< TODOCUMENT
                               ) const {
	plot(
		arg_output_file_stem,
		make_classn_stat_pair_series_list( arg_score_classn_value_lists, arg_x_classn_stat, arg_y_classn_stat ),
		arg_x_classn_stat.get_name(),
		arg_y_classn_stat.get_name(),
		arg_plot_spec
	);
}

/// \brief TODOCUMENT
void classn_stat_plotter::plot(const path                           &arg_output_file_stem, ///< TODOCUMENT
                               const score_classn_value_results_set &arg_scv_results_set,  ///< TODOCUMENT
                               const classn_stat                    &arg_x_classn_stat,    ///< TODOCUMENT
                               const classn_stat                    &arg_y_classn_stat,    ///< TODOCUMENT
                               const classn_stat_plotter_spec       &arg_plot_spec         ///< TODOCUMENT
                               ) const {
	plot(
		arg_output_file_stem,
		make_classn_stat_pair_series_list( arg_scv_results_set, arg_x_classn_stat, arg_y_classn_stat ),
		arg_x_classn_stat.get_name(),
		arg_y_classn_stat.get_name(),
		arg_plot_spec
	);
}

/// \brief TODOCUMENT
void classn_stat_plotter::plot(const path                          &arg_output_file_stem, ///< TODOCUMENT
                               const named_true_false_pos_neg_list &arg_named_tfpn,       ///< TODOCUMENT
                               const classn_stat                   &arg_x_classn_stat,    ///< TODOCUMENT
                               const classn_stat                   &arg_y_classn_stat,    ///< TODOCUMENT
                               const classn_stat_plotter_spec      &arg_plot_spec         ///< TODOCUMENT
                               ) const {
	plot(
		arg_output_file_stem,
		get_classn_stat_pair_series( arg_named_tfpn.get_list(), arg_named_tfpn.get_name(), arg_x_classn_stat, arg_y_classn_stat ),
		arg_x_classn_stat.get_name(),
		arg_y_classn_stat.get_name(),
		arg_plot_spec
	);
}


/// \brief TODOCUMENT
void classn_stat_plotter::plot(const path                     &arg_output_file_stem, ///< TODOCUMENT
                               const classn_stat_pair_series  &arg_series,           ///< TODOCUMENT
                               const string                   &arg_x_axis_label,     ///< TODOCUMENT
                               const string                   &arg_y_axis_label,     ///< TODOCUMENT
                               const classn_stat_plotter_spec &arg_plot_spec         ///< TODOCUMENT
                               ) const {
	plot(
		arg_output_file_stem,
		classn_stat_pair_series_list( { arg_series } ),
		arg_x_axis_label,
		arg_y_axis_label,
		arg_plot_spec
	);
}

/// \brief TODOCUMENT
void classn_stat_plotter::plot(const path                         &arg_output_file_stem, ///< TODOCUMENT
                               const classn_stat_pair_series_list &arg_serieses,         ///< TODOCUMENT
                               const string                       &arg_x_axis_label,     ///< TODOCUMENT
                               const string                       &arg_y_axis_label,     ///< TODOCUMENT
                               const classn_stat_plotter_spec     &arg_plot_spec         ///< TODOCUMENT
                               ) const {
	const path gnuplot_file  = replace_extension_copy( arg_output_file_stem, ".gnuplot"  );
	const path eps_file      = replace_extension_copy( arg_output_file_stem, ".eps"      );
	const path the_data_file = replace_extension_copy( arg_output_file_stem, ".data.txt" );

	const auto gnuplot_pipe = suppress_execution ? string( " > /dev/null " )
	                                             : string( " | gnuplot "   );
	Gnuplot gp( "tee " + gnuplot_file.string() + gnuplot_pipe ); // Write to an intermediate gnuplot file

	gp << "set   terminal postscript eps enhanced color\n";

//	{"#023FA5","#7D87B9","#BEC1D4","#D6BCC0","#BB7784","#FFFFFF"},
//	{"#4A6FE3","#8595E1","#B5BBE3","#E6AFB9","#E07B91","#D33F6A"},
//	{"#11C638","#8DD593","#C6DEC7","#EAD3C6","#F0B98D","#EF9708"},
//	{"#0FCFC0","#9CDED6","#D5EAE7","#F3E1EB","#F6C4E1","#F79CD4"}

	gp << "set   output " << eps_file << "\n";

//	gp << "set title  'Some ROC curves'\n";
	gp << "set xlabel '" << arg_x_axis_label << R"(' font "Helvetica,20")" << "\n";
	gp << "set ylabel '" << arg_y_axis_label << R"(' font "Helvetica,20")" << "\n";

	gp << R"(
set size square 2,2
set xtics 0,0.1
set ytics 0,0.1
set xrange [0:1]
set yrange [0:1]
set style line 11 lc rgb "#808080" lt 1
set border 3 back ls 11
set tics nomirror
set style line 12 lc rgb "#808080" lt 0 lw 1
set grid back ls 12
set key bottom right
set key font ",11"
set key spacing 0.7
set xtics font "Helvetica,18"
set ytics font "Helvetica,18"
)";

	const auto series_to_plot = get_series_to_plot_or_make_default( arg_plot_spec, arg_serieses );

	const auto main_plot_strs = transform_build<str_vec>(
		series_to_plot,
		[&] (const pair<string, str_opt> &x) {
			const auto &the_series    = classn_stat_pair_series_list_of_name( arg_serieses, x.first );
			const auto  data_filename = the_data_file.string() + lexical_cast<string>( x.first );
			const auto  legend_name   = process_legend_name_copy( the_series.get_name(), arg_plot_spec.get_tidy_up_score_based_legends() );

			const auto  series_props  = ( x.second.value_or( ""s ) );
			const auto  legend_props  = icontains( series_props, "title" ) ? string() : ( " title '" + legend_name + "'" );

			return gp.file1d( the_series.data, data_filename ) + " with lines " + series_props + legend_props;
		}
	);
	gp << "plot "
	   << boost::algorithm::join( boost::range::join( arg_plot_spec.get_pre_plot_strs(), main_plot_strs ), ", ")
	   << "\n";
	// To check: has removing the final endl from here caused any problems?

}

/// \brief TODOCUMENT
///
/// \relates classn_stat_plotter
void cath::score::plot_roc(const classn_stat_plotter         &arg_classn_stat_plotter,      ///< TODOCUMENT
                           const path                        &arg_output_file_stem,         ///< TODOCUMENT
                           const score_classn_value_list_vec &arg_score_classn_value_lists, ///< TODOCUMENT
                           const classn_stat_plotter_spec    &arg_plot_spec                 ///< TODOCUMENT
                           ) {
	plot_classn_stat<roc_rates>(
		arg_classn_stat_plotter,
		arg_output_file_stem,
		arg_score_classn_value_lists,
		arg_plot_spec
	);
}

/// \brief TODOCUMENT
///
/// \relates classn_stat_plotter
void cath::score::plot_roc(const classn_stat_plotter            &arg_classn_stat_plotter, ///< TODOCUMENT
                           const path                           &arg_output_file_stem,    ///< TODOCUMENT
                           const score_classn_value_results_set &arg_scv_results_set,     ///< TODOCUMENT
                           const classn_stat_plotter_spec       &arg_plot_spec            ///< TODOCUMENT
                           ) {
	plot_classn_stat<roc_rates>(
		arg_classn_stat_plotter,
		arg_output_file_stem,
		arg_scv_results_set,
		arg_plot_spec
	);
}

/// \brief TODOCUMENT
///
/// \relates classn_stat_plotter
void cath::score::plot_roc(const classn_stat_plotter           &arg_classn_stat_plotter, ///< TODOCUMENT
                           const path                          &arg_output_file_stem,    ///< TODOCUMENT
                           const named_true_false_pos_neg_list &arg_named_tfpn,          ///< TODOCUMENT
                           const classn_stat_plotter_spec      &arg_plot_spec            ///< TODOCUMENT
                           ) {
	plot_classn_stat<roc_rates>(
		arg_classn_stat_plotter,
		arg_output_file_stem,
		arg_named_tfpn,
		arg_plot_spec
	);
}

/// \brief TODOCUMENT
///
/// \relates classn_stat_plotter
void cath::score::plot_precision_recall(const classn_stat_plotter         &arg_classn_stat_plotter,      ///< TODOCUMENT
                                        const path                        &arg_output_file_stem,         ///< TODOCUMENT
                                        const score_classn_value_list_vec &arg_score_classn_value_lists, ///< TODOCUMENT
                                        const classn_stat_plotter_spec    &arg_plot_spec                 ///< TODOCUMENT
                                        ) {
	plot_classn_stat<precision_recall_rates>(
		arg_classn_stat_plotter,
		arg_output_file_stem,
		arg_score_classn_value_lists,
		arg_plot_spec
	);
}

/// \brief TODOCUMENT
///
/// \relates classn_stat_plotter
void cath::score::plot_precision_recall(const classn_stat_plotter            &arg_classn_stat_plotter, ///< TODOCUMENT
                                        const path                           &arg_output_file_stem,    ///< TODOCUMENT
                                        const score_classn_value_results_set &arg_scv_results_set,     ///< TODOCUMENT
                                        const classn_stat_plotter_spec       &arg_plot_spec            ///< TODOCUMENT
                                        ) {
	plot_classn_stat<precision_recall_rates>(
		arg_classn_stat_plotter,
		arg_output_file_stem,
		arg_scv_results_set,
		arg_plot_spec
	);
}

/// \brief TODOCUMENT
///
/// \relates classn_stat_plotter
void cath::score::plot_precision_recall(const classn_stat_plotter           &arg_classn_stat_plotter, ///< TODOCUMENT
                                        const path                          &arg_output_file_stem,    ///< TODOCUMENT
                                        const named_true_false_pos_neg_list &arg_named_tfpn,          ///< TODOCUMENT
                                        const classn_stat_plotter_spec      &arg_plot_spec            ///< TODOCUMENT
                                        ) {
	plot_classn_stat<precision_recall_rates>(
		arg_classn_stat_plotter,
		arg_output_file_stem,
		arg_named_tfpn,
		arg_plot_spec
	);
}
