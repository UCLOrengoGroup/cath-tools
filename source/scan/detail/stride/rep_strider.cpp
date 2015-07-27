/// \file
/// \brief The re_strider class definitions

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

#include "rep_strider.h"

#include "common/algorithm/constexpr_modulo_fns.h"
#include "exception/invalid_argument_exception.h"

using namespace cath::common;
using namespace cath::scan;
using namespace cath::scan::detail;
using namespace std;

using boost::none;
using boost::optional;

/// \brief TODOCUMENT
///
/// \pre arg_next_centre_index >= arg_entry_index else throws an invalid_argument_exception
optional<index_type> cath::scan::detail::centre_index_of_index_and_next_centre_index(const index_type &arg_entry_index,      ///< TODOCUMENT
                                                                                     const index_type &arg_co_stride,        ///< TODOCUMENT
                                                                                     const index_type &arg_next_centre_index ///< TODOCUMENT
                                                                                     ) {
	if ( arg_next_centre_index < arg_entry_index ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("centre_index_of_index_and_next_centre_index() requires that arg_next_centre_index >= arg_entry_index"));
	}
//	const auto       num_in_stride_range          = num_in_stride_neighbour_range           ( arg_co_stride );
	const auto       centre_stride_neighour_index = detail::stride_neighbour_index_of_centre( arg_co_stride );
	const index_type steps_to_next_centre         = arg_next_centre_index - arg_entry_index;
	return ( steps_to_next_centre <= centre_stride_neighour_index ) ? optional<index_type>( arg_next_centre_index           ) :
	       ( arg_entry_index      >= arg_co_stride                ) ? optional<index_type>( arg_entry_index - arg_co_stride ) :
	                                                                  none;
	}

/// \brief TODOCUMENT
///
/// \relates rep_strider
rep_rep_pair_opt cath::scan::detail::get_rep_of_indices(const rep_strider &arg_rep_strider_a, ///< TODOCUMENT
                                                        const index_type  &arg_index_a,       ///< TODOCUMENT
                                                        const rep_strider &arg_rep_strider_b, ///< TODOCUMENT
                                                        const index_type  &arg_index_b        ///< TODOCUMENT
                                                        ) {
	const auto coprime_pair = chinese_remainder_coprime_pair(
		arg_index_a,
		arg_index_b,
		arg_rep_strider_a.get_stride() + 1,
		arg_rep_strider_b.get_stride() + 1
	);
	assert( coprime_pair.first + arg_index_b == coprime_pair.second + arg_index_a);
	const auto the_co_stride = co_stride( arg_rep_strider_a.get_stride(), arg_rep_strider_b.get_stride() );
	const auto centre_index_a = centre_index_of_index_and_next_centre_index(
		arg_index_a,
		the_co_stride,
		coprime_pair.first
	);
	const auto centre_index_b = centre_index_of_index_and_next_centre_index(
		arg_index_b,
		the_co_stride,
		coprime_pair.second
	);
	if ( ! centre_index_a || ! centre_index_b ) {
		return none;
	}
	return { make_pair(
		static_cast<res_rep_index_type>( *centre_index_a / ( arg_rep_strider_a.get_stride() + 1 ) ),
		static_cast<res_rep_index_type>( *centre_index_b / ( arg_rep_strider_b.get_stride() + 1 ) )
	) };
}
