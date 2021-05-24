/// \file
/// \brief The equal_grouped_holder class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_DETAIL_EQUAL_GROUPED_HOLDER_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_DETAIL_EQUAL_GROUPED_HOLDER_HPP

#include "cath/common/boost_addenda/range/range_concept_type_aliases.hpp"

#include <functional>

namespace cath::common::detail {

	/// \brief TODOCUMENT
	template <typename FN>
	class equal_grouped_holder final {
	private:
		/// \brief TODOCUMENT
		FN unequal_function;

	public:
		explicit equal_grouped_holder(FN);
		equal_grouped_holder() = default;

		template <typename RNG>
		std::function<bool(const range_value_t<RNG> &, const range_value_t<RNG> &)> get_function() const;
	};
	/// \brief TODOCUMENT
	template <typename FN>
	equal_grouped_holder<FN>::equal_grouped_holder(FN prm_unequal_function ///< TODOCUMENT
	                                               ) : unequal_function( prm_unequal_function ) {
	}

	/// \brief TODOCUMENT
	///
	/// Will probably need to make this a compile-time switch based on whether FN is of type nullptr_t
	template <typename FN>
	template <typename RNG>
	std::function<bool(const range_value_t<RNG> &, const range_value_t<RNG> &)> equal_grouped_holder<FN>::get_function() const {
		return { unequal_function };
	}

	/// \brief TODOCUMENT
	///
	/// Will probably need to make this a compile-time switch based on whether FN is of type nullptr_t
	template <>
	template <typename RNG>
	std::function<bool(const range_value_t<RNG> &, const range_value_t<RNG> &)> equal_grouped_holder<std::nullptr_t>::get_function() const {
		return { std::not_equal_to<range_value_t<RNG>>() };
	}

} // namespace cath::common::detail

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_DETAIL_EQUAL_GROUPED_HOLDER_HPP
