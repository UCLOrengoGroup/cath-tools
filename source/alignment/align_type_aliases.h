/// \file
/// \brief Header for the type_aliases in the cath::align namespace

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef ALIGN_TYPE_ALIASES_H_INCLUDED
#define ALIGN_TYPE_ALIASES_H_INCLUDED

#include <boost/optional/optional_fwd.hpp>

#include "alignment/dyn_prog_align/detail/path_step.h"

#include <map>
#include <utility>
#include <vector>

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
		using opt_aln_posn = boost::optional<aln_posn_type>;
		/// \brief TODOCUMENT
		using opt_aln_posn_vec = std::vector<opt_aln_posn>;
		/// \brief TODOCUMENT
		using opt_aln_posn_vec_vec = std::vector<opt_aln_posn_vec>;

		namespace detail {
			/// \brief TODOCUMENT
			using path_step_score_map = std::map<  path_step, score_type >;
			/// \brief TODOCUMENT
			using path_step_score_pair = std::pair< path_step, score_type >;

			/// \brief TODOCUMENT
			using bool_aln_pair = std::pair<bool, alignment>;
		}

		/// \brief TODOCUMENT
		using score_alignment_pair = std::pair<score_type, align::alignment>;

		/// \brief TODOCUMENT
		using size_size_alignment_tuple = std::tuple<size_t, size_t, alignment>;
		/// \brief TODOCUMENT
		using size_size_alignment_tuple_vec = std::vector<size_size_alignment_tuple>;

		/// \brief TODOCUMENT
		using opt_alignment = boost::optional<alignment>;

		/// \brief TODOCUMENT
		using alignment_vec = std::vector<alignment>;
	}
}


#endif
