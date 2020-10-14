/// \file
/// \brief The score_classn_value class definitions

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

#include "score_classn_value.hpp"

#include <sstream>
#include <utility>

using namespace ::cath::score;
using namespace ::std;

/// \brief Ctor for score_classn_value
score_classn_value::score_classn_value(const double &prm_score_value,          ///< TODOCUMENT
                                       const bool   &prm_instance_is_positive, ///< TODOCUMENT
                                       string        prm_instance_label        ///< TODOCUMENT
                                       ) : score_value          { prm_score_value                 },
                                           instance_is_positive { prm_instance_is_positive        },
                                           instance_label       { std::move( prm_instance_label ) } {
}

/// \brief TODOCUMENT
const double & score_classn_value::get_score_value() const {
	return score_value;
}

/// \brief TODOCUMENT
const bool & score_classn_value::get_instance_is_positive() const {
	return instance_is_positive;
}

/// \brief TODOCUMENT
const string & score_classn_value::get_instance_label() const {
	return instance_label;
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value
ostream & cath::score::operator<<(ostream                  &prm_ostream,           ///< TODOCUMENT
                                  const score_classn_value &prm_score_classn_value ///< TODOCUMENT
                                  ) {
	ostringstream out_ss;
	out_ss << "score_classn_value[score:";
	out_ss << prm_score_classn_value.get_score_value();
	out_ss << ", label:";
	out_ss << prm_score_classn_value.get_instance_label();
	out_ss << ", is_positive:";
	out_ss << boolalpha << prm_score_classn_value.get_instance_is_positive();
	out_ss << "]";

	prm_ostream << out_ss.str();
	return prm_ostream;
}
