/// \file
/// \brief The matrix_plot class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_MATRIX_PLOTTER_MATRIX_PLOT_HPP
#define _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_MATRIX_PLOTTER_MATRIX_PLOT_HPP

#include <boost/filesystem/path.hpp>

#include "cath/alignment/dyn_prog_align/detail/return_path_matrix.hpp"
#include "cath/alignment/dyn_prog_align/detail/score_accumulation_matrix.hpp"
#include "cath/alignment/dyn_prog_align/dyn_prog_score_source/dyn_prog_score_source.hpp"
#include "cath/ssap/windowed_matrix.hpp"

namespace cath {
	namespace align {
		namespace detail {

			/// \brief TODOCUMENT
			template <typename P>
			void matrix_plot(const boost::filesystem::path &prm_output_stem, ///< TODOCUMENT
			                 const dyn_prog_score_source   &prm_scorer       ///< TODOCUMENT
			                 ) {
				const size_t length_a = prm_scorer.get_length_a();
				const size_t length_b = prm_scorer.get_length_b();
				P plotter( length_a, length_b, get_window_width_for_full_matrix(length_a, length_b) );
				plotter.plot_scores( prm_scorer );
				plotter.finish( prm_output_stem );
			}

			/// \brief TODOCUMENT
			template <typename P>
			void matrix_plot(const boost::filesystem::path   &prm_output_stem,              ///< TODOCUMENT
			                 const dyn_prog_score_source     &prm_scorer,                   ///< TODOCUMENT
			                 const return_path_matrix        &prm_return_path_matrix,       ///< TODOCUMENT
			                 const score_accumulation_matrix &prm_score_accumulation_matrix ///< TODOCUMENT
			                 ) {
				const size_t length_a     = prm_return_path_matrix.get_length_a();
				const size_t length_b     = prm_return_path_matrix.get_length_b();
				const size_t window_width = prm_return_path_matrix.get_window_width();
				P plotter( length_a, length_b, window_width );
				plotter.plot_scores( prm_scorer );
				plotter.plot_return_path_matrix( prm_return_path_matrix );
				plotter.plot_accumulated_scores( prm_score_accumulation_matrix );
				plotter.finish( prm_output_stem );
			}
		} // namespace detail
	} // namespace align
} // namespace cath

#endif
