/// \file
/// \brief The pair_scatter_plotter class definitions

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

#include "pair_scatter_plotter.h"

//#include "common/third_party_code/gnuplot-iostream.h"

//using namespace boost::filesystem;
//using namespace cath::score;
//using namespace std;

//using gnuplotio::Gnuplot;

//void pair_scatter_plotter::plot(const path               &arg_output_file_stem, ///< TODOCUMENT
//                                const doub_doub_pair_vec &arg_data,             ///< TODOCUMENT
//                                const string             &arg_x_axis_label,     ///< TODOCUMENT
//                                const string             &arg_y_axis_label      ///< TODOCUMENT
//                                ) const {
//	const path gnuplot_file  = arg_output_file_stem.parent_path() / ( path(arg_output_file_stem.filename()).string() + ".gnuplot"  );
//	const path eps_file      = arg_output_file_stem.parent_path() / ( path(arg_output_file_stem.filename()).string() + ".eps"      );
//	const path the_data_file = arg_output_file_stem.parent_path() / ( path(arg_output_file_stem.filename()).string() + ".data.txt" );
//
//	const auto gnuplot_pipe = suppress_execution ? string( " > /dev/null " )
//	                                             : string( " | gnuplot "   );
//	Gnuplot gp( "tee " + gnuplot_file.string() + gnuplot_pipe ); // Write to an intermediate gnuplot file
//
//	gp << "set   terminal postscript eps enhanced color\n";
//
////	{"#023FA5","#7D87B9","#BEC1D4","#D6BCC0","#BB7784","#FFFFFF"},
////	{"#4A6FE3","#8595E1","#B5BBE3","#E6AFB9","#E07B91","#D33F6A"},
////	{"#11C638","#8DD593","#C6DEC7","#EAD3C6","#F0B98D","#EF9708"},
////	{"#0FCFC0","#9CDED6","#D5EAE7","#F3E1EB","#F6C4E1","#F79CD4"}
//
//	gp << "set   output " << eps_file << "\n";
//
////	gp << "set title  'Some ROC curves'\n";
//	gp << "set xlabel '" << arg_x_axis_label << R"(' font "Helvetica,20")" << "\n";
//	gp << "set ylabel '" << arg_y_axis_label << R"(' font "Helvetica,20")" << "\n";
//
//	gp << R"(
//set size square 2,2
//# set xtics 0,0.1
//# set ytics 0,0.1
//# set xrange [0:1]
//# set yrange [0:1]
//set style line 11 lc rgb "#808080" lt 1
//set border 3 back ls 11
//set tics nomirror
//set style line 12 lc rgb "#808080" lt 0 lw 1
//set grid back ls 12
//set key bottom right
//set key font ",11"
//set key spacing 0.7
//set xtics font "Helvetica,18"
//set ytics font "Helvetica,18"
//)";
//
//
//	gp << "plot " << gp.file1d( arg_data, the_data_file.string() ) + " with points \n";
//	// To check: has removing the final endl from here caused any problems?
//}
