/// \file
/// \brief The ssap_score_accuracy class header

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

#ifndef SSAP_SCORE_ACCURACY_H_INCLUDED
#define SSAP_SCORE_ACCURACY_H_INCLUDED

#include "common/algorithm/constexpr_is_uniq.h"

#include <array>
#include <iosfwd>
#include <map>

namespace cath {
	namespace score {

		/// \brief Whether to use high accuracy (floating point numbers; no explicit rounding) or low accuracy (ints)
		///
		/// The low accuracy option is mainly provided to permit regression testing against the old SSAP score code.
		enum class ssap_score_accuracy {
			LOW, ///< Use high accuracy calculations (floating point numbers; no explicit rounding)
			HIGH ///< Use low accuracy calculations (ints)
		};

		/// \brief TODOCUMENT
		static constexpr std::array<ssap_score_accuracy, 2> all_ssap_score_accuracies = {{
			ssap_score_accuracy::LOW,
			ssap_score_accuracy::HIGH
		}};

		static_assert( common::constexpr_is_uniq( all_ssap_score_accuracies ), "all_ssap_score_accuracies shouldn't contain repeated values" );

		/// \brief TODOCUMENT
		static constexpr size_t num_ssap_score_accuracies = std::tuple_size< decltype( all_ssap_score_accuracies ) >::value;

		/// \brief TODOCUMENT
		struct name_of_ssap_score_accuracy final {
			static std::map<ssap_score_accuracy, std::string> get();
		};

		std::ostream & operator<<(std::ostream &,
		                          const ssap_score_accuracy &);

	}
}

#endif
