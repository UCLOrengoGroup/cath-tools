/// \file
/// \brief The score_accumulation_matrix class header

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

#ifndef _CATH_TOOLS_SOURCE_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_SCORE_ACCUMULATION_MATRIX_H
#define _CATH_TOOLS_SOURCE_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_SCORE_ACCUMULATION_MATRIX_H

#include "alignment/align_type_aliases.h"
#include "common/container/vector_of_vector.h"

#include <vector>

namespace cath { namespace align { class dyn_prog_score_source; } }
namespace cath { namespace align { namespace gap { class gap_penalty; } } }

namespace cath {
	namespace align {
		namespace detail {
			class return_path_matrix;

			/// \brief TODOCUMENT
			class score_accumulation_matrix final {
			public:
				/// \brief TODOCUMENT
				using size_type = score_vec_of_vec::size_type;

			private:
				/// \brief TODOCUMENT
				score_vec_of_vec scores;

				/// \brief TODOCUMENT
				size_type        window_width;

				static void check_length(const size_type &);

				static void check_index_against_length(const size_type &,
				                                       const size_type &,
				                                       const bool &);

				void check_indices(const size_type &,
				                   const size_type &,
				                   const bool &arg_permit_one_extra = false) const;

				void initialise(const size_type &,
				                const size_type &,
				                const size_type &);

			public:
				score_accumulation_matrix(const size_type &,
				                          const size_type &,
				                          const size_type &);

				void reset(const size_type &,
				           const size_type &,
				           const size_type &);

				size_type get_length_a() const;
				size_type get_length_b() const;
				size_type get_window_width() const;

				void set_score_towards_end_at_point(const score_accumulation_matrix::size_type &,
				                                    const score_accumulation_matrix::size_type &,
				                                    const score_type &);

				score_type get_score_towards_end_at_point(const score_accumulation_matrix::size_type &,
				                                          const score_accumulation_matrix::size_type &) const;
			};

			score_accumulation_matrix make_uninitialised_score_accumulation_matrix();

			score_type get_score_towards_end_from_point_plus_path_step(const score_accumulation_matrix &,
			                                                           const score_accumulation_matrix::size_type &,
			                                                           const score_accumulation_matrix::size_type &,
			                                                           const path_step &);

			path_step_score_map get_total_scores_of_path_steps_from_point(const score_accumulation_matrix &,
			                                                              const return_path_matrix &,
			                                                              const gap::gap_penalty &,
			                                                              const dyn_prog_score_source &,
			                                                              const size_t &,
			                                                              const size_t &);

		} // namespace detail
	} // namespace align
} // namespace cath

#endif
