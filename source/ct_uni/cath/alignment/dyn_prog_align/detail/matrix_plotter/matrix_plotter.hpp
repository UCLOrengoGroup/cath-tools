/// \file
/// \brief The matrix_plotter class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_MATRIX_PLOTTER_MATRIX_PLOTTER_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_MATRIX_PLOTTER_MATRIX_PLOTTER_HPP

#include <filesystem>
#include <vector>

#include "cath/common/type_aliases.hpp"

namespace cath {
	namespace align {
		class dyn_prog_score_source;

		namespace detail {
			class return_path_matrix;
			class score_accumulation_matrix;

			/// \brief Provide ABC interface for classes that provide concrete details for specific ways of outputting (eg to Gnuplot)
			class matrix_plotter {
			private:
				/// \brief The length of the first  entry (the number of rows    in the matrix)
				size_t length_a;

				/// \brief The length of the second entry (the number of columns in the matrix)
				size_t length_b;

				bool simplified_windowless = false;

				/// \brief The window width to use when plotting
				size_t window_width = 0;

				virtual void do_plot_scores(doub_vec_vec &,
				                            const double &,
				                            const double &) = 0;
				virtual void do_plot_minor_corner_arrow(const size_t &,
				                                        const size_t &,
				                                        const size_t &,
				                                        const size_t &) = 0;
				virtual void do_write_centre_score(const size_t &,
				                                   const size_t &,
				                                   const double &) = 0;
				virtual void do_write_corner_score(const size_t &,
				                                   const size_t &,
				                                   const double &) = 0;
				virtual void do_finish(const ::std::filesystem::path &) const = 0;

				void check_lengths_match(const size_t &,
				                         const size_t &) const;
				void check_window_width_matches(const size_t &);

			protected:
				size_t get_length_a() const;
				size_t get_length_b() const;
				size_t get_window_width() const;

			public:
				matrix_plotter(const size_t &,
				               const size_t &);

				matrix_plotter(const size_t &,
				               const size_t &,
				               const size_t &);
				virtual ~matrix_plotter() noexcept = default;

				matrix_plotter(const matrix_plotter &) = default;
				matrix_plotter(matrix_plotter &&) noexcept = default;
				matrix_plotter & operator=(const matrix_plotter &) = default;
				matrix_plotter & operator=(matrix_plotter &&) noexcept = default;

				void plot_scores(const dyn_prog_score_source &);
				void plot_return_path_matrix(const return_path_matrix &);
				void plot_accumulated_scores(const score_accumulation_matrix &);
				void finish(const ::std::filesystem::path &) const;
			};

		} // namespace detail
	} // namespace align
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_MATRIX_PLOTTER_MATRIX_PLOTTER_HPP
