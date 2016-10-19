/// \file
/// \brief The selected_pair class header

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

#ifndef _CATH_TOOLS_SOURCE_SSAP_SELECTED_PAIR_H
#define _CATH_TOOLS_SOURCE_SSAP_SELECTED_PAIR_H

#include <boost/operators.hpp>

#include "common/type_aliases.h"

namespace cath {

	/// \brief Represent a pair of entries (probably residues) and their score
	///
	/// This is a pretty simple class that's used to help simplify the code that
	/// selects the top few pairs (ie select_pairs() and update_best_pair_selections() )
	class selected_pair final : public boost::less_than_comparable<selected_pair> {
		size_t     index_a;
		size_t     index_b;
		score_type score;

	public:
		selected_pair(const size_t &,
		              const size_t &,
		              const score_type &);

		size_t     get_index_a() const;
		size_t     get_index_b() const;
		score_type get_score()   const;
	};

	bool operator<(const selected_pair &,
	               const selected_pair &);

} // namespace cath

#endif
