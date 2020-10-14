/// \file
/// \brief The hit_extras class definitions

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

#include "hit_extras.hpp"

#include "common/exception/invalid_argument_exception.hpp"

using namespace cath::common;

using std::string;

/// \brief Generate a string describing the specified hit_extra_cat
///
/// \relates hit_extra_cat
string cath::rslv::to_string(const hit_extra_cat &prm_hit_extra_cat ///< The hit_extra_cat to describe
                             ) {
	switch (prm_hit_extra_cat) {
		case ( hit_extra_cat::ALND_RGNS ) : { return "aligned-regions"; }
		case ( hit_extra_cat::COND_EVAL ) : { return "cond-evalue";     }
		case ( hit_extra_cat::INDP_EVAL ) : { return "indp-evalue";     }
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of hit_extra_cat not recognised whilst converting to_string()"));
}