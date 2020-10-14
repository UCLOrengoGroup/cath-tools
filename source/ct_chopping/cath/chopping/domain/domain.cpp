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

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/equal.hpp>

#include "cath/chopping/domain/domain.hpp"
#include "cath/chopping/region/region.hpp"
#include "cath/common/algorithm/are_same.hpp"
#include "cath/common/boost_addenda/range/front.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/optional/make_optional_if.hpp"

#include <string>

using namespace cath::chop;
using namespace cath::common;
using namespace cath;
using namespace std::literals::string_literals;

using boost::adaptors::transformed;
using boost::algorithm::join;
using boost::none;
using boost::range::equal;
using std::equal_to;
using std::ostream;
using std::string;

/// \brief TODOCUMENT
void domain::sanity_check() const {
	if ( ! are_same( segments, equal_to<>{}, [] (const region &x) { return get_residue_locating( x ); } ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Cannot construct a domain from segments with different methods of locating residues (names and/or indices)"
		));
	}
}

/// \brief Ctor for domain
domain::domain(region_vec prm_segments, ///< TODOCUMENT
               string     prm_domain_id ///< TODOCUMENT
               ) : segments { std::move( prm_segments  ) },
                   domain_id{ std::move( prm_domain_id ) } {
	sanity_check();
}

/// \brief Ctor for domain
domain::domain(region_vec prm_segments ///< TODOCUMENT
               ) : segments{ std::move( prm_segments ) } {
	sanity_check();
}

/// \brief TODOCUMENT
size_t domain::num_segments() const {
	return segments.size();
}

///// \brief TODOCUMENT
//region domain::operator[](const size_t &prm_index ///< TODOCUMENT
//                          ) {
//	return segments[ prm_index ];
//}

/// \brief TODOCUMENT
const region & domain::operator[](const size_t &prm_index ///< TODOCUMENT
                                  ) const {
	return segments[ prm_index ];
}

/// \brief TODOCUMENT
void domain::set_opt_domain_id(const str_opt &prm_opt_domain_id ///< TODOCUMENT
                               ) {
	domain_id = prm_opt_domain_id;
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

/// \brief Get any regions from the specified optional domain or none if none
///
/// \relates domain
///
/// \relatesalso region
region_vec_opt cath::chop::get_regions_opt(const domain_opt &prm_domain_opt ///< The optional domain from which any regions should be extracted
                                           ) {
	return make_optional_if_fn(
		static_cast<bool>( prm_domain_opt ),
		[&] {
			return region_vec{
				common::cbegin( *prm_domain_opt ),
				common::cend  ( *prm_domain_opt ),
			};
		}
	);
}

/// \brief Return whether the two specified domains are identical
///
/// \relates domain
bool cath::chop::operator==(const domain &prm_lhs, ///, The first  domain to compare
                            const domain &prm_rhs  ///, The second domain to compare
                            ) {
	return (
		prm_lhs.get_opt_domain_id() == prm_rhs.get_opt_domain_id()
		&&
		equal( prm_lhs, prm_rhs )
	);
}

/// \brief TODOCUMENT
///
/// \relates domain
bool cath::chop::has_domain_id(const domain &prm_domain ///< TODOCUMENT
                               ) {
	return static_cast<bool>( prm_domain.get_opt_domain_id() );
}

/// \brief TODOCUMENT
///
/// \relates domain
string cath::chop::get_domain_id(const domain &prm_domain ///< TODOCUMENT
                                 ) {
	return *prm_domain.get_opt_domain_id();
}

/// \brief TODOCUMENT
///
/// \relates domain
residue_locating_opt cath::chop::get_residue_locating(const domain &prm_domain ///< TODOCUMENT
                                                      ) {
	// If there are no segments then this domain doesn't locate any residues so return none
	if ( prm_domain.num_segments() == 0 ) {
		return none;
	}
	// Otherwise return the residue_locating method of the first segment
	return get_residue_locating( front( prm_domain ) );
}


/// \brief Generate a string describing the specified domain
///
/// \relates domain
string cath::chop::to_string(const domain &prm_domain ///< The domain to describe
                             ) {
	return "domain["
		+ (
			has_domain_id( prm_domain )
			? ( "name:" + get_domain_id( prm_domain ) + ", " )
			: ""s
		)
		+ join(
			prm_domain
				| transformed( [] (const region &x) { return to_string( x ); } ),
			","
		)
		+ "]";
}

/// \brief Insert a description of the specified domain into the specified ostream
///
/// \relates domain
ostream & cath::chop::operator<<(ostream      &prm_os,    ///< The ostream into which the description should be inserted
                                 const domain &prm_domain ///< The domain to describe
                                 ) {
	prm_os << to_string( prm_domain );
	return prm_os;
}

