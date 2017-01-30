/// \file
/// \brief The domain class definitions

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

#include <boost/range.hpp>

#include "common/boost_addenda/range/front.hpp"
#include "chopping/domain/domain.hpp"
#include "chopping/region/region.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "exception/invalid_argument_exception.hpp"

using namespace cath;
using namespace cath::chop;
using namespace cath::common;
using namespace std;

using boost::none;

/// \brief TODOCUMENT
void domain::sanity_check() const {
	if ( ! segments.empty() ) {
		const residue_locating seg_res_locating = get_residue_locating( segments.front() );
		for (const region &segment : segments) {
			if ( get_residue_locating( segment ) != seg_res_locating ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot construct a domain from segments with different methods of locating residues (names and/or indices)"));
			}
		}
	}
}

/// \brief Ctor for domain
domain::domain(const region_vec &arg_segments, ///< TODOCUMENT
               const string     &arg_domain_id///< TODOCUMENT
               ) : segments ( arg_segments  ),
                   domain_id( arg_domain_id ) {
	sanity_check();
}

/// \brief Ctor for domain
domain::domain(const region_vec &arg_segments ///< TODOCUMENT
               ) : segments( arg_segments ) {
	sanity_check();
}

/// \brief TODOCUMENT
size_t domain::num_segments() const {
	return segments.size();
}

///// \brief TODOCUMENT
//region domain::operator[](const size_t &arg_index ///< TODOCUMENT
//                          ) {
//	return segments[ arg_index ];
//}

/// \brief TODOCUMENT
const region & domain::operator[](const size_t &arg_index ///< TODOCUMENT
                                  ) const {
	return segments[ arg_index ];
}

/// \brief TODOCUMENT
void domain::set_opt_domain_id(const str_opt &arg_opt_domain_id ///< TODOCUMENT
                               ) {
	domain_id = arg_opt_domain_id;
}

/// \brief TODOCUMENT
const str_opt & domain::get_opt_domain_id() const {
	return domain_id;
}

///// \brief TODOCUMENT
//domain::iterator domain::begin() {
//	return std::begin( segments );
//}

///// \brief TODOCUMENT
//domain::iterator domain::end() {
//	return std::end( segments );
//}

/// \brief TODOCUMENT
domain::const_iterator domain::begin() const {
	return common::cbegin( segments );
}

/// \brief TODOCUMENT
domain::const_iterator domain::end() const {
	return common::cend( segments );
}

/// \brief TODOCUMENT
///
/// \relates domain
bool cath::chop::has_domain_id(const domain &arg_domain ///< TODOCUMENT
                               ) {
	return static_cast<bool>( arg_domain.get_opt_domain_id() );
}

/// \brief TODOCUMENT
///
/// \relates domain
string cath::chop::get_domain_id(const domain &arg_domain ///< TODOCUMENT
                                 ) {
	return *arg_domain.get_opt_domain_id();
}

/// \brief TODOCUMENT
///
/// \relates domain
residue_locating_opt cath::chop::get_residue_locating(const domain &arg_domain ///< TODOCUMENT
                                                      ) {
	// If there are no segments then this domain doesn't locate any residues so return none
	if ( arg_domain.num_segments() == 0 ) {
		return none;
	}
	// Otherwise return the residue_locating method of the first segment
	return get_residue_locating( front( arg_domain ) );
}

