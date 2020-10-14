/// \file
/// \brief The gnuplot_matrix_plotter class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_MATRIX_PLOTTER_GNUPLOT_MATRIX_PLOTTER_HPP
#define _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_MATRIX_PLOTTER_GNUPLOT_MATRIX_PLOTTER_HPP

#include "alignment/dyn_prog_align/detail/matrix_plotter/matrix_plotter.hpp"

#include <sstream>

namespace cath {
	namespace align {
		namespace detail {

			/// \brief TODOCUMENT
			class gnuplot_matrix_plotter final : public matrix_plotter  {
			private:
				void do_plot_scores(doub_vec_vec &,
				                    const double &,
				                    const double &) final;
				void do_plot_minor_corner_arrow(const size_t &,
				                                const size_t &,
				                                const size_t &,
				                                const size_t &) final;
				void do_write_centre_score(const size_t &,
				                           const size_t &,
				                           const double &) final;
				void do_write_corner_score(const size_t &,
				                           const size_t &,
				                           const double &) final;
				void do_finish(const boost::filesystem::path &) const final;

				/// \brief TODOCUMENT
				std::ostringstream preplot_commands;

				/// \brief TODOCUMENT
				doub_vec_vec       scores;

				/// \brief TODOCUMENT
				size_t             counter = 1;

			public:
				gnuplot_matrix_plotter(const size_t &,
				                       const size_t &,
				                       const size_t &);
				gnuplot_matrix_plotter(const size_t &,
				                       const size_t &);
			};

		} // namespace detail
	} // namespace align
} // namespace cath

#endif
