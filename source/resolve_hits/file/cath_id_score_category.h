/// \file
/// \brief The cath_id_score_category header

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

#ifndef CATH_ID_SCORE_CATEGORY_H_INCLUDED
#define CATH_ID_SCORE_CATEGORY_H_INCLUDED

#include <boost/utility/string_ref_fwd.hpp>

namespace cath {
	namespace rslv {

		/// \brief Represent the type of match ID, wrt the CATH-specific scoring that Jon uses
		enum class cath_id_score_category {
			NORMAL,     ///< A normal match (use this as the default if unsure)
			DC_TYPE,    ///< An ID like dc_c869189e57e572c71376c2f3dfe7dc9c that is handled differently
			LATER_ROUND ///< An ID like 2fcwB01_round_2 or 2ezwA00_round_3 (but not 2ffkB00_round_1 because that's first-round)
		};

		cath_id_score_category cath_score_category_of_id(const boost::string_ref &,
		                                                 const bool &);

		std::string to_string(const cath_id_score_category &);

		std::ostream & operator<<(std::ostream &,
		                          const cath_id_score_category &);

	}
}

#endif
