/// \file
/// \brief The pick_random_paiR header

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Tools project and then tweaked, eg namespaced in cath)
/// Copyright (C) 2007, Tony Lewis
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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_RANDOM_PICK_RANDOM_PAIR_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_RANDOM_PICK_RANDOM_PAIR_HPP

#include "cath/common/exception/invalid_argument_exception.hpp"

#include <random>

namespace cath::common {

	/// \brief Randomly pick a pair of distinct integers in the specified range
	///
	/// To ensure distinctness, this picks the second value from a one-smaller range and
	/// then increments if the result's >= the first value.
	///
	/// If there's reason, this could be generalised in a few ways:
	///  * add a version that selects n-values without replacement
	///  * add a version that takes a range and returns a vector of
	///    (references to?) values in the range
	template <typename T>
	inline std::pair<T, T> pick_random_pair(const T      &prm_min, ///< TODOCUMENT
	                                        const T      &prm_max, ///< TODOCUMENT
	                                        std::mt19937 &prm_rng  ///< TODOCUMENT
	                                        ) {
		if ( ! ( prm_max > prm_min ) ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot pick random pair if the max isn't greater than the min"));
		}

		const T value_a = std::uniform_int_distribution<T>{ prm_min, prm_max     }( prm_rng );
		const T raw_b   = std::uniform_int_distribution<T>{ prm_min, prm_max - 1 }( prm_rng );
		const T value_b = ( raw_b >= value_a ) ? raw_b + 1 : raw_b;

		return {
			value_a,
			value_b
		};
	}

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_RANDOM_PICK_RANDOM_PAIR_HPP
