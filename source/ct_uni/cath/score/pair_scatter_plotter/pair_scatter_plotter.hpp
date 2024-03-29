/// \file
/// \brief The pair_scatter_plotter class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_PAIR_SCATTER_PLOTTER_PAIR_SCATTER_PLOTTER_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_PAIR_SCATTER_PLOTTER_PAIR_SCATTER_PLOTTER_HPP

#include <filesystem>

#include <gnuplot-iostream.h>

#include "cath/common/boost_addenda/filesystem/replace_extension_copy.hpp"
#include "cath/common/type_aliases.hpp"

namespace cath::score {

	/// \brief TODOCUMENT
	class pair_scatter_plotter final {
	private:
		/// \brief TODOCUMENT
		bool suppress_execution = false;

	public:

		template <typename T>
		void plot(const ::std::filesystem::path &prm_output_file_stem, ///< TODOCUMENT
		          const std::vector<T>          &prm_data,             ///< TODOCUMENT
		          const std::string             &prm_x_axis_label,     ///< TODOCUMENT
		          const std::string             &prm_y_axis_label      ///< TODOCUMENT
		          ) const {
			const auto gnuplot_file  = common::replace_extension_copy( prm_output_file_stem, ".gnuplot"  );
			const auto eps_file      = common::replace_extension_copy( prm_output_file_stem, ".eps"      );
			const auto the_data_file = common::replace_extension_copy( prm_output_file_stem, ".data.txt" );

			const auto gnuplot_pipe = suppress_execution ? std::string( " > /dev/null " )
			                                             : std::string( " | gnuplot "   );
			gnuplotio::Gnuplot gp( "tee " + gnuplot_file.string() + gnuplot_pipe ); // Write to an intermediate gnuplot file

			gp << "set   terminal postscript eps enhanced color\n";

		//	{"#023FA5","#7D87B9","#BEC1D4","#D6BCC0","#BB7784","#FFFFFF"},
		//	{"#4A6FE3","#8595E1","#B5BBE3","#E6AFB9","#E07B91","#D33F6A"},
		//	{"#11C638","#8DD593","#C6DEC7","#EAD3C6","#F0B98D","#EF9708"},
		//	{"#0FCFC0","#9CDED6","#D5EAE7","#F3E1EB","#F6C4E1","#F79CD4"}

			gp << "set   output " << eps_file << "\n";

		//	gp << "set title  'Some ROC curves'\n";
			gp << "set xlabel '" << prm_x_axis_label << R"(' font "Helvetica,20")" << "\n";
			gp << "set ylabel '" << prm_y_axis_label << R"(' font "Helvetica,20")" << "\n";

			gp << R"(
set size square 2,2
# set xtics 0,0.1
# set ytics 0,0.1
# set xrange [0:1]
# set yrange [0:1]
set style line 11 lc rgb "#808080" lt 1
set border 3 back ls 11
set tics nomirror
set style line 12 lc rgb "#808080" lt 0 lw 1

rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)

set grid back ls 12
set key bottom right
set key font ",11"
set key spacing 0.7
set xtics font "Helvetica,18"
set ytics font "Helvetica,18"
		)";

			gp << "plot " << gp.file1d( prm_data, the_data_file.string() ) + " with points \n";
		}

	};

} // namespace cath::score

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_PAIR_SCATTER_PLOTTER_PAIR_SCATTER_PLOTTER_HPP
