/// \file
/// \brief Header for the type_aliases in the cath::align namespace

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_ALIGN_TYPE_ALIASES_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_ALIGN_TYPE_ALIASES_HPP

#include <map>
#include <optional>
#include <utility>
#include <vector>

#include <boost/config.hpp> /// \todo Come a resolution for Boost Trac tickets 12142 & 12179, remove this #include

#include "cath/alignment/dyn_prog_align/detail/path_step.hpp"

namespace cath {
	namespace align {
		class alignment;

		/// \brief TODOCUMENT
		using aln_posn_type = size_t;
		/// \brief TODOCUMENT
		using aln_posn_vec = std::vector<aln_posn_type>;
		/// \brief TODOCUMENT
		using aln_posn_vec_vec = std::vector<aln_posn_vec>;
		/// \brief TODOCUMENT
		using aln_posn_opt = ::std::optional<aln_posn_type>;
		/// \brief TODOCUMENT
		using aln_posn_opt_vec = std::vector<aln_posn_opt>;
		/// \brief TODOCUMENT
		using aln_posn_opt_vec_vec = std::vector<aln_posn_opt_vec>;

		namespace detail {
			/// \brief TODOCUMENT
			using path_step_score_map = std::map<  path_step, score_type >;
			/// \brief TODOCUMENT
			using path_step_score_pair = std::pair< path_step, score_type >;

			/// \brief TODOCUMENT
			using bool_aln_pair = std::pair<bool, alignment>;
		} // namespace detail

		/// \brief TODOCUMENT
		using score_alignment_pair = std::pair<score_type, align::alignment>;

		/// \brief TODOCUMENT
		using size_size_alignment_tuple = std::tuple<size_t, size_t, alignment>;
		/// \brief TODOCUMENT
		using size_size_alignment_tuple_vec = std::vector<size_size_alignment_tuple>;

		/// \brief TODOCUMENT
		using alignment_opt = ::std::optional<alignment>;

		/// \brief TODOCUMENT
		using alignment_vec = std::vector<alignment>;
	} // namespace align
} // namespace cath


#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_ALIGN_TYPE_ALIASES_HPP
