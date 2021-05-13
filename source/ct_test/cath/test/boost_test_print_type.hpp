/// \file
/// \brief Header containing useful boost_test_print_type()s

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_TEST_BOOST_TEST_PRINT_TYPE_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_TEST_BOOST_TEST_PRINT_TYPE_HPP

#include <optional>
#include <tuple>
#include <utility>

#include "cath/common/string/cath_to_string.hpp"

namespace std {

	/// It's rather unpleasant to be putting these into the std namespace, which is bad practice.
	/// Unfortunately, this is the way to get things to work in the Boost Test library ATM.
	///
	/// \TODO Move this into a module for test code - this shouldn't be included with any non-test code

	/// \brief A boost_test_print_type for a nullopt_t
	///
	/// \param prm_os    The ostream to which the value should be written
	/// \param prm_value The nullopt_t to write
	inline ::std::ostream &boost_test_print_type( ::std::ostream &prm_os, [[maybe_unused]] const ::std::nullopt_t &prm_value ) {
		prm_os << "nullopt";
		return prm_os;
	}

	/// \brief A boost_test_print_type for an optional
	///
	/// \param prm_os    The ostream to which the value should be written
	/// \param prm_value The optional to write
	template <typename T>
	inline ::std::ostream &boost_test_print_type( ::std::ostream &prm_os, const ::std::optional<T> &prm_value ) {
		prm_os << ::cath::common::cath_to_string( prm_value );
		return prm_os;
	}

	/// \brief A boost_test_print_type for a pair
	///
	/// \param prm_os   The ostream to which the value should be written
	/// \param prm_pair The pair to write
	template <typename T, typename U>
	inline ::std::ostream &boost_test_print_type( ::std::ostream &prm_os, const ::std::pair<T, U> &prm_pair ) {
		prm_os << ::cath::common::cath_to_string( prm_pair );
		return prm_os;
	}

	/// \brief A boost_test_print_type for a tuple
	///
	/// \param prm_os    The ostream to which the value should be written
	/// \param prm_tuple The tuple to write
	template <typename... Ts>
	inline ::std::ostream &boost_test_print_type( ::std::ostream &prm_os, const ::std::tuple<Ts...> &prm_tuple ) {
		prm_os << ::cath::common::cath_to_string( prm_tuple );
		return prm_os;
	}

} // namespace std

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_TEST_BOOST_TEST_PRINT_TYPE_HPP
