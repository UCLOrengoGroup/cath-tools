/// \file
/// \brief The chopping class definitions

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

#include "chopping.hpp"

#include <boost/range.hpp>

#include "chopping/domain/domain.hpp"
#include "chopping/region/region.hpp"
#include "common/cpp14/cbegin_cend.hpp"

#include <utility>

using namespace cath;
using namespace cath::chop;

/// \brief TODOCUMENT
void chopping::sanity_check() const {

}

/// \brief Ctor for chopping
chopping::chopping(domain_vec prm_domains,  ///< TODOCUMENT
                   region_vec prm_fragments ///< TODOCUMENT
                   ) : domains  { std::move( prm_domains   ) },
                       fragments{ std::move( prm_fragments ) } {
}

/// \brief TODOCUMENT
size_t chopping::num_domains() const {
	return domains.size();
}

/// \brief TODOCUMENT
size_t chopping::num_fragments() const {
	return fragments.size();
}

/// \brief TODOCUMENT
const region & chopping::get_fragment_of_index(const size_t &prm_index ///< TODOCUMENT
                                               ) const {
	return fragments[ prm_index ];
}

/// \brief (A const-reference to) the domain at the specified index within this chopping
///
/// \param prm_index The index of the domain of interest
const domain &chopping::operator[]( const size_t &prm_index ) const {
	return domains[ prm_index ];
}

/// \brief TODOCUMENT
auto chopping::begin() const -> const_iterator {
	return common::cbegin( domains );
}

/// \brief TODOCUMENT
auto chopping::end() const -> const_iterator {
	return common::cend( domains );
}
