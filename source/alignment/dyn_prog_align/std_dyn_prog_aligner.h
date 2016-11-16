/// \file
/// \brief The std_dyn_prog_aligner class header

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

#ifndef _CATH_TOOLS_SOURCE_ALIGNMENT_DYN_PROG_ALIGN_STD_DYN_PROG_ALIGNER_H
#define _CATH_TOOLS_SOURCE_ALIGNMENT_DYN_PROG_ALIGN_STD_DYN_PROG_ALIGNER_H

#include "alignment/dyn_prog_align/detail/return_path_matrix.h"
#include "alignment/dyn_prog_align/detail/score_accumulation_matrix.h"
#include "alignment/dyn_prog_align/dyn_prog_aligner.h"
#include "common/type_aliases.h"

namespace cath {
	namespace align {

		/// \brief TODOCUMENT
		class std_dyn_prog_aligner final : public dyn_prog_aligner {
		private:
			/// \brief TODOCUMENT
			mutable detail::return_path_matrix        the_return_path        = detail::make_uninitialised_return_path_matrix();

			/// \brief TODOCUMENT
			mutable detail::score_accumulation_matrix the_accumulated_scores = detail::make_uninitialised_score_accumulation_matrix();

			virtual std::unique_ptr<dyn_prog_aligner> do_clone() const override final;

			virtual score_alignment_pair do_align(const dyn_prog_score_source &,
			                                      const gap::gap_penalty &,
			                                      const size_type &) const override final;

			detail::path_step choose_path_step(const detail::path_step_score_map &,
			                                   const size_t &,
			                                   const size_t &,
			                                   const size_t &,
			                                   const size_t &) const;
		};
	} // namespace align
} // namespace cath
#endif
