/// \file
/// \brief The gnuplot_matrix_plotter class definitions

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

#include "gnuplot_matrix_plotter.hpp"

#include <filesystem>

#include <boost/numeric/conversion/cast.hpp>

#include <gnuplot-iostream.h>

#include "cath/alignment/dyn_prog_align/detail/return_path_matrix.hpp"
#include "cath/common/boost_addenda/filesystem/replace_extension_copy.hpp"

using namespace ::cath::align::detail;
using namespace ::cath::common;
using namespace ::std;

using ::boost::numeric_cast;
using ::std::filesystem::path;

/// \brief TODOCUMENT
void gnuplot_matrix_plotter::do_plot_scores(doub_vec_vec &prm_scores,        ///< TODOCUMENT
                                            const double &/*prm_min_score*/, ///< TODOCUMENT
                                            const double &/*prm_max_score*/  ///< TODOCUMENT
                                            ) {
	scores = prm_scores;
}

/// \brief TODOCUMENT
void gnuplot_matrix_plotter::do_plot_minor_corner_arrow(const size_t &prm_from_x, ///< TODOCUMENT
                                                        const size_t &prm_from_y, ///< TODOCUMENT
                                                        const size_t &prm_to_x,   ///< TODOCUMENT
                                                        const size_t &prm_to_y    ///< TODOCUMENT
                                                        ) {
	const size_t max_length = max( get_length_a(), get_length_b() );
	preplot_commands << "set   arrow from ";
	preplot_commands << numeric_cast<double>( prm_from_x ) - 0.5;
	preplot_commands << ", ";
	preplot_commands << numeric_cast<double>( prm_from_y ) - 0.5;
	preplot_commands << " to ";
	preplot_commands << numeric_cast<double>( prm_to_x   ) - 0.5;
	preplot_commands << ", ";
	preplot_commands << numeric_cast<double>( prm_to_y   ) - 0.5;
//	preplot_commands << " front head filled\n";
	preplot_commands << " size 0.20, 15 front head filled linewidth " << 10.0 / numeric_cast<double>( max_length + 1) << "\n";
}

/// \brief TODOCUMENT
void gnuplot_matrix_plotter::do_write_centre_score(const size_t &prm_x,    ///< TODOCUMENT
                                                   const size_t &prm_y,    ///< TODOCUMENT
                                                   const double &prm_score ///< TODOCUMENT
                                                   ) {
	const size_t max_length = max( get_length_a(), get_length_b() );
	preplot_commands << "set label ";
	preplot_commands << counter++;
	preplot_commands << " \"";
	preplot_commands << prm_score;
	preplot_commands << "\" at first ";
	preplot_commands << numeric_cast<double>( prm_x );
	preplot_commands << ", ";
	preplot_commands << numeric_cast<double>( prm_y );
	preplot_commands << " centre front textcolor rgbcolor \"black\" nopoint font \", " << 50.0 / numeric_cast<double>( max_length + 1 ) << "\"\n";
}

/// \brief TODOCUMENT
void gnuplot_matrix_plotter::do_write_corner_score(const size_t &prm_x,    ///< TODOCUMENT
                                                   const size_t &prm_y,    ///< TODOCUMENT
                                                   const double &prm_score ///< TODOCUMENT
                                                   ) {
	const size_t max_length = max( get_length_a(), get_length_b() );
	preplot_commands << "set label ";
	preplot_commands << counter++;
	preplot_commands << " \"";
	preplot_commands << prm_score;
	preplot_commands << "\" at first ";
	preplot_commands << numeric_cast<double>( prm_x ) - 0.5;
	preplot_commands << ", ";
	preplot_commands << numeric_cast<double>( prm_y ) - 0.5;
	preplot_commands << " centre front textcolor rgbcolor \"black\" nopoint offset first 0.200, -0.085 font \", " << 50.0 / numeric_cast<double>( max_length + 1 ) << "\"\n";
}

/// \brief TODOCUMENT
void gnuplot_matrix_plotter::do_finish(const path &prm_output_stem ///< TODOCUMENT
                                       ) const {
	const path gnuplot_file  = replace_extension_copy( prm_output_stem, ".gnuplot"  );
	const path eps_file      = replace_extension_copy( prm_output_stem, ".eps"      );
	const path the_data_file = replace_extension_copy( prm_output_stem, ".data.txt" );
	Gnuplot gp("tee " + gnuplot_file.string() + " | gnuplot"); // Write to an intermediate gnuplot file

	const size_t lcl_length_a   = get_length_a();
	const size_t lcl_length_b   = get_length_b();
	const size_t lcl_max_length = max(lcl_length_a, lcl_length_b);

	gp << "set   terminal postscript color\n";
	gp << "set   output " << eps_file << "\n";
	gp << "set   size ratio " << numeric_cast<double>(lcl_length_a) / numeric_cast<double>(lcl_length_b) << "\n";
	gp << "set   palette defined (0 \"white\", 1 0.1 0.1 0.1)\n";
//	gp << "set   palette defined (0 \"white\", 0.2 0.2 0.2 0.2, 1 \"black\")\n";
	gp << "set   autoscale x2fix\n";
	gp << "set   autoscale yfix\n";
	gp << "set   x2tics 1 rotate font \", " << 350.0 / numeric_cast<double>( lcl_max_length + 1 ) << "\"\n";
	gp << "set   ytics  1 rotate font \", " << 350.0 / numeric_cast<double>( lcl_max_length + 1 ) << "\"\n";
	gp << "set   tics scale 0,0.001\n";
	gp << "unset xtics\n";
	gp << "set   mx2tics 2\n";
	gp << "set   mytics  2\n";
	gp << "unset colorbox\n";
	gp << "set   yrange [] reverse\n";
	gp << "set   grid front mx2tics mytics linewidth " << 20.0 / numeric_cast<double>( lcl_max_length + 1) << " lt -1 lc rgb 'white'\n";

	gp << preplot_commands.str();

	cerr << "About to plot with size " << scores.size() << "\n";

	gp << "plot " << gp.file1d(scores, the_data_file.string()) << " matrix with image notitle\n";

	// To check: has removing the final endl from here caused any problems?
}

/// \brief Ctor for gnuplot_matrix_plotter
gnuplot_matrix_plotter::gnuplot_matrix_plotter(const size_t &prm_length_a,    ///< TODOCUMENT
                                               const size_t &prm_length_b,    ///< TODOCUMENT
                                               const size_t &prm_window_width ///< TODOCUMENT
                                               ) : matrix_plotter( prm_length_a, prm_length_b, prm_window_width ) {
}

/// \brief Ctor for gnuplot_matrix_plotter
gnuplot_matrix_plotter::gnuplot_matrix_plotter(const size_t &prm_length_a,    ///< TODOCUMENT
                                               const size_t &prm_length_b     ///< TODOCUMENT
                                               ) : matrix_plotter( prm_length_a, prm_length_b ) {
}
