/// \file
/// \brief The gap_penalty class definitions

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

#include "gap_penalty.hpp"

using namespace ::cath;
//using namespace ::cath::align;
using namespace ::cath::align::gap;
// using namespace ::std;

/// \brief Ctor for gap_penalty
gap_penalty::gap_penalty(const score_type &prm_open_gap_penalty,  ///< TODOCUMENT
                         const score_type &prm_extend_gap_penalty ///< TODOCUMENT
                         ) : open_gap_penalty  ( prm_open_gap_penalty   ),
                             extend_gap_penalty( prm_extend_gap_penalty ) {
}

/// \brief TODOCUMENT
score_type gap_penalty::get_open_gap_penalty() const {
	return open_gap_penalty;
}

/// \brief TODOCUMENT
score_type gap_penalty::get_extend_gap_penalty() const {
	return extend_gap_penalty;
}
