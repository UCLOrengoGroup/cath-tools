/// \file
/// \brief The hit_score_type class definitions

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

#include "hit_score_type.hpp"

#include "common/exception/invalid_argument_exception.hpp"

using std::ostream;
using std::string;

/// \brief Generate a string describing the specified hit_score_type
///
/// \relates hit_score_type
string cath::rslv::to_string(const hit_score_type &prm_hit_score_type ///< The hit_score_type to describe
                             ) {
	switch ( prm_hit_score_type ) {
		case ( hit_score_type::FULL_EVALUE ) : { return "evalue"    ; }
		case ( hit_score_type::BITSCORE    ) : { return "bitscore"  ; }
		case ( hit_score_type::CRH_SCORE   ) : { return "crh-value" ; }
	}
	BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Value of hit_score_type not recognised whilst converting to_string()"));
}

/// \brief Insert a description of the specified hit_score_type into the specified ostream
///
/// \relates hit_score_type
ostream & cath::rslv::operator<<(ostream              &prm_os,            ///< The ostream into which the description should be inserted
                                 const hit_score_type &prm_hit_score_type ///< The hit_score_type to describe
                                 ) {
	prm_os << to_string( prm_hit_score_type );
	return prm_os;
}