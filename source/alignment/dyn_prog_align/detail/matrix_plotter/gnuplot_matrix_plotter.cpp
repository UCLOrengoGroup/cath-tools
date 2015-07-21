/// \file
/// \brief The gnuplot_matrix_plotter class definitions

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "gnuplot_matrix_plotter.h"

#include <boost/filesystem/path.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "alignment/dyn_prog_align/detail/return_path_matrix.h"
#include "common/third_party_code/gnuplot-iostream.h"

using namespace boost::filesystem;
using namespace cath::align::detail;
using namespace std;

using boost::numeric_cast;

/// \brief TODOCUMENT
void gnuplot_matrix_plotter::do_plot_scores(doub_vec_vec &arg_scores,        ///< TODOCUMENT
                                            const double &/*arg_min_score*/, ///< TODOCUMENT
                                            const double &/*arg_max_score*/  ///< TODOCUMENT
                                            ) {
	scores = arg_scores;
}

/// \brief TODOCUMENT
void gnuplot_matrix_plotter::do_plot_minor_corner_arrow(const size_t &arg_from_x, ///< TODOCUMENT
                                                        const size_t &arg_from_y, ///< TODOCUMENT
                                                        const size_t &arg_to_x,   ///< TODOCUMENT
                                                        const size_t &arg_to_y    ///< TODOCUMENT
                                                        ) {
	const size_t length_a   = get_length_a();
	const size_t length_b   = get_length_b();
	const size_t max_length = max(length_a, length_b);
	preplot_commands << "set   arrow from ";
	preplot_commands << numeric_cast<double>( arg_from_x ) - 0.5;
	preplot_commands << ", ";
	preplot_commands << numeric_cast<double>( arg_from_y ) - 0.5;
	preplot_commands << " to ";
	preplot_commands << numeric_cast<double>( arg_to_x   ) - 0.5;
	preplot_commands << ", ";
	preplot_commands << numeric_cast<double>( arg_to_y   ) - 0.5;
//	preplot_commands << " front head filled\n";
	preplot_commands << " size 0.20, 15 front head filled linewidth " << 10.0 / numeric_cast<double>( max_length + 1) << "\n";
}

/// \brief TODOCUMENT
void gnuplot_matrix_plotter::do_write_centre_score(const size_t &arg_x,    ///< TODOCUMENT
                                                   const size_t &arg_y,    ///< TODOCUMENT
                                                   const double &arg_score ///< TODOCUMENT
                                                   ) {
	const size_t length_a   = get_length_a();
	const size_t length_b   = get_length_b();
	const size_t max_length = max(length_a, length_b);
	preplot_commands << "set label ";
	preplot_commands << counter++;
	preplot_commands << " \"";
	preplot_commands << arg_score;
	preplot_commands << "\" at first ";
	preplot_commands << numeric_cast<double>( arg_x );
	preplot_commands << ", ";
	preplot_commands << numeric_cast<double>( arg_y );
	preplot_commands << " centre front textcolor rgbcolor \"black\" nopoint font \", " << 50.0 / numeric_cast<double>( max_length + 1 ) << "\"\n";
}

/// \brief TODOCUMENT
void gnuplot_matrix_plotter::do_write_corner_score(const size_t &arg_x,    ///< TODOCUMENT
                                                   const size_t &arg_y,    ///< TODOCUMENT
                                                   const double &arg_score ///< TODOCUMENT
                                                   ) {
	const size_t length_a   = get_length_a();
	const size_t length_b   = get_length_b();
	const size_t max_length = max(length_a, length_b);
	preplot_commands << "set label ";
	preplot_commands << counter++;
	preplot_commands << " \"";
	preplot_commands << arg_score;
	preplot_commands << "\" at first ";
	preplot_commands << numeric_cast<double>( arg_x ) - 0.5;
	preplot_commands << ", ";
	preplot_commands << numeric_cast<double>( arg_y ) - 0.5;
	preplot_commands << " centre front textcolor rgbcolor \"black\" nopoint offset first 0.200, -0.085 font \", " << 50.0 / numeric_cast<double>( max_length + 1 ) << "\"\n";
}

/// \brief TODOCUMENT
void gnuplot_matrix_plotter::do_finish(const path &arg_output_stem ///< TODOCUMENT
                                       ) const {
	const path gnuplot_file = arg_output_stem.parent_path() / ( path(arg_output_stem.filename()).string() + ".gnuplot"  );
	const path eps_file     = arg_output_stem.parent_path() / ( path(arg_output_stem.filename()).string() + ".eps"      );
	const path the_data_file    = arg_output_stem.parent_path() / ( path(arg_output_stem.filename()).string() + ".data.txt" );
	Gnuplot gp("tee " + gnuplot_file.string() + " | gnuplot"); // Write to an intermediate gnuplot file

	const size_t length_a   = get_length_a();
	const size_t length_b   = get_length_b();
	const size_t max_length = max(length_a, length_b);

	gp << "set   terminal postscript color\n";
	gp << "set   output " << eps_file << "\n";
	gp << "set   size ratio " << numeric_cast<double>(length_a) / numeric_cast<double>(length_b) << "\n";
	gp << "set   palette defined (0 \"white\", 1 0.1 0.1 0.1)\n";
//	gp << "set   palette defined (0 \"white\", 0.2 0.2 0.2 0.2, 1 \"black\")\n";
	gp << "set   autoscale x2fix\n";
	gp << "set   autoscale yfix\n";
	gp << "set   x2tics 1 rotate font \", " << 350.0 / numeric_cast<double>( max_length + 1 ) << "\"\n";
	gp << "set   ytics  1 rotate font \", " << 350.0 / numeric_cast<double>( max_length + 1 ) << "\"\n";
	gp << "set   tics scale 0,0.001\n";
	gp << "unset xtics\n";
	gp << "set   mx2tics 2\n";
	gp << "set   mytics  2\n";
	gp << "unset colorbox\n";
	gp << "set   yrange [] reverse\n";
	gp << "set   grid front mx2tics mytics linewidth " << 20.0 / numeric_cast<double>( max_length + 1) << " lt -1 lc rgb 'white'\n";

	gp << preplot_commands.str();

	cerr << "About to plot with size " << scores.size() << endl;

	gp << "plot " << gp.file1d(scores, the_data_file.string()) << " matrix with image notitle" << endl;
}

/// \brief Ctor for gnuplot_matrix_plotter
gnuplot_matrix_plotter::gnuplot_matrix_plotter(const size_t &arg_length_a,    ///< TODOCUMENT
                                               const size_t &arg_length_b,    ///< TODOCUMENT
                                               const size_t &arg_window_width ///< TODOCUMENT
                                               ) : matrix_plotter( arg_length_a, arg_length_b, arg_window_width ) {
}

/// \brief Ctor for gnuplot_matrix_plotter
gnuplot_matrix_plotter::gnuplot_matrix_plotter(const size_t &arg_length_a,    ///< TODOCUMENT
                                               const size_t &arg_length_b     ///< TODOCUMENT
                                               ) : matrix_plotter( arg_length_a, arg_length_b ) {
}
