/// \file
/// \brief The bifur_hbond_list class definitions

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

#include "bifur_hbond_list.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <string>

using namespace cath::sec;

using boost::adaptors::transformed;
using boost::algorithm::join;
using std::ostream;
using std::string;

/// \brief Generate a string describing the specified hbond_half_opt
///
/// \relates hbond_half_opt
string cath::sec::to_string(const hbond_half_opt &arg_hbond_half_opt ///< The hbond_half_opt to describe
                            ) {
	return arg_hbond_half_opt
		? (
			  "("
			+ ::std::to_string( arg_hbond_half_opt->index )
			+ ", "
			+ ::std::to_string( arg_hbond_half_opt->energy )
			+ ")"
		)
		: "-";
}

/// \brief Generate a string describing the specified bifur_hbond
///
/// \relates bifur_hbond
string cath::sec::to_string(const bifur_hbond &arg_bifur_hbond ///< The bifur_hbond to describe
                            ) {
	return
		"bifur_hbond[nh_1st:"
		+ to_string( arg_bifur_hbond.get_bound_pair_for_this_nh().first  )
		+ ", nh_2nd"
		+ to_string( arg_bifur_hbond.get_bound_pair_for_this_nh().second )
		+ ", co_1st:"
		+ to_string( arg_bifur_hbond.get_bound_pair_for_this_co().first  )
		+ ", co_2nd:"
		+ to_string( arg_bifur_hbond.get_bound_pair_for_this_co().second )
		+ "]";
}

/// \brief Insert a description of the specified bifur_hbond into the specified ostream
///
/// \relates bifur_hbond
ostream & cath::sec::operator<<(ostream           &arg_os,         ///< The ostream into which the description should be inserted
                                const bifur_hbond &arg_bifur_hbond ///< The bifur_hbond to describe
                                ) {
	arg_os << to_string( arg_bifur_hbond );
	return arg_os;
}

/// \brief Generate a string describing the specified bifur_hbond_list
///
/// \relates bifur_hbond_list
string cath::sec::to_string(const bifur_hbond_list &arg_bifur_hbond_list ///< The bifur_hbond_list to describe
                            ) {
	return
		"bifur_hbond_list[\n\t"
		+ join(
			arg_bifur_hbond_list
				| transformed( [] (const bifur_hbond &x) { return to_string( x ); } ),
			"\n\t"
		)
		+ "]";
}

/// \brief Insert a description of the specified bifur_hbond_list into the specified ostream
///
/// \relates bifur_hbond_list
std::ostream & cath::sec::operator<<(ostream                &arg_os,              ///< The ostream into which the description should be inserted
                                     const bifur_hbond_list &arg_bifur_hbond_list ///< The bifur_hbond_list to describe
                                     ) {
	arg_os << to_string( arg_bifur_hbond_list );
	return arg_os;
}
