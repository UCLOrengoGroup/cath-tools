/// \file
/// \brief The return_path_matrix class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_RETURN_PATH_MATRIX_HPP
#define _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_RETURN_PATH_MATRIX_HPP

#include "cath/alignment/align_type_aliases.hpp"
#include "cath/alignment/dyn_prog_align/detail/path_step.hpp"
#include "cath/common/type_aliases.hpp"

#include <vector>

namespace cath { namespace align { class alignment; } }
namespace cath { namespace align { namespace gap { class gap_penalty; } } }

namespace cath {
	namespace align {
		namespace detail {

			/// \brief TODOCUMENT
			///
			/// NOTE: This is NOT symmetric: it always records return paths from each given point
			///       *to the ends* of the two sequences, *not to their starts*.
			///       (because, it uses the fact that there's a unique path leading out of each point towards the end
			///        whereas there may be up to 3 different paths leading out from a given point to the start, as shown below)
			///
			///     #.....#
			///     .\    |
			///     . \   |
			///     .  \  |
			///     .   \ |
			///     .    4V
			///     #---->#
			///
			/// \todo Construct a space-efficient windowed matrix template and then implement
			///       this class in terms of that
			class return_path_matrix final {
			public:
				/// \brief TODOCUMENT
				using size_type = path_step_vec_vec::size_type;

			private:
				/// \brief TODOCUMENT
				path_step_vec_vec return_path;

				/// \brief TODOCUMENT
				size_type         window_width;

				static void check_length(const size_type &);

				void initialise(const size_type &,
				                const size_type &,
				                const size_type &);

			public:
				return_path_matrix(const size_type &,
				                   const size_type &,
				                   const size_type &);

				void reset(const size_type &,
				           const size_type &,
				           const size_type &);

				size_type get_length_a() const;
				size_type get_length_b() const;
				size_type get_window_width() const;

				void set_path_step_towards_end_at_point(const return_path_matrix::size_type &,
				                                        const return_path_matrix::size_type &,
				                                        const path_step &);

				path_step get_path_step_towards_end_at_point(const return_path_matrix::size_type &,
				                                             const return_path_matrix::size_type &) const;
			};

			return_path_matrix make_uninitialised_return_path_matrix();

			alignment make_alignment(const return_path_matrix &);

			path_step get_path_dirn_towards_end_from_point_plus_path_step(const return_path_matrix &,
			                                                              const return_path_matrix::size_type &,
			                                                              const return_path_matrix::size_type &,
			                                                              const path_step &);

			score_type get_gap_penalty_for_path_step_from_point(const return_path_matrix &,
			                                                    const gap::gap_penalty &,
			                                                    const path_step &,
			                                                    const return_path_matrix::size_type &,
			                                                    const return_path_matrix::size_type &);

			size_size_pair get_b_window_start_and_stop_for_a_index(const return_path_matrix &,
			                                                       const return_path_matrix::size_type &);

			std::ostream & operator<<(std::ostream &,
			                          const return_path_matrix &);

			score_type max_path_step_score(const path_step_score_map &);
		} // namespace detail
	} // namespace align
} // namespace cath

#endif
