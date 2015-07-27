/// \file
/// \brief The ssap_score_accuracy class definitions

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

#include "ssap_score_accuracy.h"

#include <iostream>

using namespace cath::score;
using namespace std;

/// \brief TODOCUMENT
map<ssap_score_accuracy, string> name_of_ssap_score_accuracy::get() {
	return {
		{ ssap_score_accuracy::LOW,  "low_accuracy"  },
		{ ssap_score_accuracy::HIGH, "high_accuracy" }
	};
}

/// \brief TODOCUMENT
///
/// \relates ssap_score_accuracy
std::ostream & cath::score::operator<<(ostream                   &arg_os,                 ///< TODOCUMENT
                                       const ssap_score_accuracy &arg_ssap_score_accuracy ///< TODOCUMENT
                                       ) {
	arg_os << name_of_ssap_score_accuracy::get().at( arg_ssap_score_accuracy );
	return arg_os;
}
