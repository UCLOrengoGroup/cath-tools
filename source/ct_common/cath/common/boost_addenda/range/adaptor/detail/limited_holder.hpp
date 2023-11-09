/// \file
/// \brief The limited_holder class header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_DETAIL_LIMITED_HOLDER_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_DETAIL_LIMITED_HOLDER_HPP

#include "cath/common/boost_addenda/range/range_concept_type_aliases.hpp"

#include <functional>

namespace cath::common::detail {

	/// \brief TODOCUMENT
	class limited_holder final {
	private:
		/// \brief TODOCUMENT
		size_t max_num_elements;

	public:
		explicit limited_holder(const size_t &);

		[[nodiscard]] const size_t &get_max_num_elements() const;
	};

	/// \brief TODOCUMENT
	inline limited_holder::limited_holder(const size_t &prm_max_num_elements ///< TODOCUMENT
	                                      ) : max_num_elements( prm_max_num_elements ) {
	}

	/// \brief TODOCUMENT
	inline const size_t & limited_holder::get_max_num_elements() const {
		return max_num_elements;
	}

} // namespace cath::common::detail

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_DETAIL_LIMITED_HOLDER_HPP
