/// \file
/// \brief The ssap_code_dyn_prog_aligner class header

/// \copyright
/// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 1989, Orengo Group, University College London
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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_SSAP_CODE_DYN_PROG_ALIGNER_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_SSAP_CODE_DYN_PROG_ALIGNER_HPP

#include <boost/tuple/tuple.hpp>

#include "cath/alignment/dyn_prog_align/dyn_prog_aligner.hpp"
#include "cath/common/type_aliases.hpp"

// clang-format off
namespace cath::align { class alignment; }
namespace cath::align { class dyn_prog_score_source; }
// clang-format on

namespace cath::align {

	/// \brief TODOCUMENT
	///
	/// WARNING: This uses static data so using separate instances is not enough
	///          to make this thread-safe.
	class ssap_code_dyn_prog_aligner final : public dyn_prog_aligner {
	  private:
		[[nodiscard]] std::unique_ptr<dyn_prog_aligner> do_clone() const final;

		using size_size_int_int_score_tuple = std::tuple<size_t, size_t, int, int, score_type>;

		static size_size_int_int_score_tuple score_matrix( const dyn_prog_score_source &,
		                                                   const score_type &,
		                                                   const size_type &,
		                                                   int_vec_vec & );

		static alignment traceback( const size_t &, const size_t &, const int &, const int &, const int &, const int &, const int_vec_vec & );

		static void traceback_recursive( alignment &,
		                                 const int &,
		                                 const int &,
		                                 const int &,
		                                 const int &,
		                                 const int &,
		                                 const int &,
		                                 const int_vec_vec & );

		[[nodiscard]] score_alignment_pair do_align( const dyn_prog_score_source &,
		                                             const gap::gap_penalty &,
		                                             const size_type & ) const final;
	};

} // namespace cath::align

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_SSAP_CODE_DYN_PROG_ALIGNER_HPP
