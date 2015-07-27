/// \file
/// \brief The res_pair_index_dirn_criterion header

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

#ifndef RES_PAIR_INDEX_DIRN_CRITERION_H_INCLUDED
#define RES_PAIR_INDEX_DIRN_CRITERION_H_INCLUDED

#include <iosfwd>

namespace cath {
	namespace scan {

		/// \brief Represent whether two res_pairs' directions of indices (increasing/decreasing) must match
		enum class res_pair_index_dirn_criterion : bool {
			MUST_MATCH,    ///< Two res_pairs' directions of indices (increasing/decreasing) must match
			NEED_NOT_MATCH ///< Two res_pairs' directions of indices (increasing/decreasing) need not match
		};

		std::ostream & operator<<(std::ostream &,
		                          const res_pair_index_dirn_criterion &);

	}
}

#endif
