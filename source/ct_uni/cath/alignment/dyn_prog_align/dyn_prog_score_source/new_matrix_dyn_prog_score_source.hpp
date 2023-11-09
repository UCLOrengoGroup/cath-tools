/// \file
/// \brief The new_matrix_dyn_prog_score_source class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_SCORE_SOURCE_NEW_MATRIX_DYN_PROG_SCORE_SOURCE_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_SCORE_SOURCE_NEW_MATRIX_DYN_PROG_SCORE_SOURCE_HPP

#include <vector>

#include "cath/alignment/dyn_prog_align/dyn_prog_score_source/dyn_prog_score_source.hpp"

namespace cath::align {

	/// \brief Concrete class that provides scores for dynamic-programming aligning by
	///        retrieving them from a matrix
	class new_matrix_dyn_prog_score_source final : public dyn_prog_score_source {
	private:
		const float_score_vec_vec &matrix;
		const size_t               length_a;
		const size_t               length_b;

		[[nodiscard]] size_t     do_get_length_a() const final;
		[[nodiscard]] size_t     do_get_length_b() const final;
		[[nodiscard]] score_type do_get_score( const size_t &, const size_t & ) const final;

	  public:
		new_matrix_dyn_prog_score_source(const float_score_vec_vec &,
		                                 const size_t &,
		                                 const size_t &);
	};

} // namespace cath::align

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_SCORE_SOURCE_NEW_MATRIX_DYN_PROG_SCORE_SOURCE_HPP
