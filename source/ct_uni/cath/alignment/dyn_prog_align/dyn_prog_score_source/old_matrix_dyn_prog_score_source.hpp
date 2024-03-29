/// \file
/// \brief The old_matrix_dyn_prog_score_source class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_SCORE_SOURCE_OLD_MATRIX_DYN_PROG_SCORE_SOURCE_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_SCORE_SOURCE_OLD_MATRIX_DYN_PROG_SCORE_SOURCE_HPP

#include "cath/alignment/dyn_prog_align/dyn_prog_score_source/dyn_prog_score_source.hpp"
#include "cath/common/container/vector_of_vector.hpp"

#include <vector>

namespace cath::align {

	/// \brief Concrete class that provides scores for dynamic-programming aligning by
	///        retrieving them from a matrix
	class old_matrix_dyn_prog_score_source final : public dyn_prog_score_source {
	private:
		/// \brief TODOCUMENT
		std::reference_wrapper<const score_vec_of_vec> matrix;

		/// \brief TODOCUMENT
		const size_t         length_a;

		/// \brief TODOCUMENT
		const size_t         length_b;

		/// \brief TODOCUMENT
		const size_t         window;

		[[nodiscard]] size_t     do_get_length_a() const final;
		[[nodiscard]] size_t     do_get_length_b() const final;
		[[nodiscard]] score_type do_get_score( const size_t &, const size_t & ) const final;

	  public:
		old_matrix_dyn_prog_score_source(const score_vec_of_vec &,
		                                 const size_t &,
		                                 const size_t &,
		                                 const size_t &);
		old_matrix_dyn_prog_score_source(const score_vec_of_vec &&,
		                                 const size_t &,
		                                 const size_t &,
		                                 const size_t &) = delete;
	};

} // namespace cath::align

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_SCORE_SOURCE_OLD_MATRIX_DYN_PROG_SCORE_SOURCE_HPP
