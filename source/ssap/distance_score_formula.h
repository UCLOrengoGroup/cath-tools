/// \file
/// \brief The distance_score_formula header

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

#ifndef DISTANCE_SCORE_FORMULA_H_INCLUDED
#define DISTANCE_SCORE_FORMULA_H_INCLUDED

#include "common/algorithm/constexpr_is_uniq.h"

#include <array>
#include <iosfwd>
#include <map>

namespace cath {

	/// \brief TODOCUMENT
	enum class distance_score_formula {
		FROM_SSAP_PAPER,       ///< TODOCUMENT
		USED_IN_PREVIOUS_CODE, ///< TODOCUMENT
		SIMPLIFIED             ///< TODOCUMENT
	};

	/// \brief TODOCUMENT
	static constexpr std::array<distance_score_formula, 3> all_distance_score_formulae = {{
		distance_score_formula::FROM_SSAP_PAPER,
		distance_score_formula::USED_IN_PREVIOUS_CODE,
		distance_score_formula::SIMPLIFIED
	}};

	static_assert( common::constexpr_is_uniq( all_distance_score_formulae ), "all_distance_score_formulae shouldn't contain repeated values" );

	/// \brief TODOCUMENT
	static constexpr size_t num_distance_score_formulae = std::tuple_size< decltype( all_distance_score_formulae ) >::value;

	/// \brief TODOCUMENT
	struct name_of_distance_score_formula final {
		static std::map<distance_score_formula, std::string> get();
	};

	std::ostream & operator<<(std::ostream &,
	                          const distance_score_formula &);

}

#endif
