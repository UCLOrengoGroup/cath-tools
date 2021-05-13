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

#include "rep_strider.hpp"

#include <optional>

#include "cath/common/algorithm/constexpr_modulo_fns.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"

using namespace ::cath::common;
using namespace ::cath::scan;
using namespace ::cath::scan::detail;

using ::std::make_pair;
using ::std::nullopt;
using ::std::optional;

/// \brief TODOCUMENT
///
/// \pre prm_next_centre_index >= prm_entry_index else throws an invalid_argument_exception
optional<index_type> cath::scan::detail::centre_index_of_index_and_next_centre_index(const index_type &prm_entry_index,      ///< TODOCUMENT
                                                                                     const index_type &prm_co_stride,        ///< TODOCUMENT
                                                                                     const index_type &prm_next_centre_index ///< TODOCUMENT
                                                                                     ) {
	if ( prm_next_centre_index < prm_entry_index ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("centre_index_of_index_and_next_centre_index() requires that prm_next_centre_index >= prm_entry_index"));
	}
//	const auto       num_in_stride_range          = num_in_stride_neighbour_range           ( prm_co_stride );
	const auto       centre_stride_neighour_index = detail::stride_neighbour_index_of_centre( prm_co_stride );
	const index_type steps_to_next_centre         = prm_next_centre_index - prm_entry_index;
	return ( steps_to_next_centre <= centre_stride_neighour_index ) ? optional<index_type>( prm_next_centre_index           ) :
	       ( prm_entry_index      >= prm_co_stride                ) ? optional<index_type>( prm_entry_index - prm_co_stride ) :
	                                                                  nullopt;
	}

/// \brief TODOCUMENT
///
/// \relates rep_strider
rep_rep_pair_opt cath::scan::detail::get_rep_of_indices(const rep_strider &prm_rep_strider_a, ///< TODOCUMENT
                                                        const index_type  &prm_index_a,       ///< TODOCUMENT
                                                        const rep_strider &prm_rep_strider_b, ///< TODOCUMENT
                                                        const index_type  &prm_index_b        ///< TODOCUMENT
                                                        ) {
	const auto coprime_pair = chinese_remainder_coprime_pair(
		prm_index_a,
		prm_index_b,
		prm_rep_strider_a.get_stride() + 1,
		prm_rep_strider_b.get_stride() + 1
	);
	assert( coprime_pair.first + prm_index_b == coprime_pair.second + prm_index_a);
	const auto the_co_stride = co_stride( prm_rep_strider_a.get_stride(), prm_rep_strider_b.get_stride() );
	const auto centre_index_a = centre_index_of_index_and_next_centre_index(
		prm_index_a,
		the_co_stride,
		coprime_pair.first
	);
	const auto centre_index_b = centre_index_of_index_and_next_centre_index(
		prm_index_b,
		the_co_stride,
		coprime_pair.second
	);
	if ( ! centre_index_a || ! centre_index_b ) {
		return nullopt;
	}
	return { make_pair(
		static_cast<res_rep_index_type>( *centre_index_a / ( prm_rep_strider_a.get_stride() + 1 ) ),
		static_cast<res_rep_index_type>( *centre_index_b / ( prm_rep_strider_b.get_stride() + 1 ) )
	) };
}
