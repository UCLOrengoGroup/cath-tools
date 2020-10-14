/// \file
/// \brief The exception_is_equivalent class header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_EXCEPTION_EXCEPTION_IS_EQUIVALENT_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_EXCEPTION_EXCEPTION_IS_EQUIVALENT_HPP

#include <boost/exception/all.hpp>

namespace cath {
	namespace test {

		/// \brief A simple predicate to establish whether a boost::exception is equivalent to this one
		///       (can be cast to this one's type and same what() string).
		///
		/// The motivation for this class is to be used as the last argument in BOOST_CHECK_EXCEPTION()
		///
		/// Note: don't try to create a exception_is_equivalent with a temporary exceptionToCompareTo because it will have been destroyed
		///       by the time the operator() is called.
		template <typename T>
		class exception_is_equivalent final {
			const T &exception_to_compare_to;

		public:
			/// \brief Ctor for exception_is_equivalent
			explicit exception_is_equivalent(const T &prm_exception_to_compare_to
			                                 ) : exception_to_compare_to(prm_exception_to_compare_to) {
			}

			/// Function call operator (operator()) to compare the two boost::exception objects
			bool operator()(const boost::exception &prm_exception_to_compare ///< The boost::exception that is to be compared
			                ) const {
				try {
					const T &exception_to_compare = dynamic_cast<const T &>(prm_exception_to_compare);
					const std::string exception_to_compare_to_what_string(exception_to_compare_to.what());
					const std::string exception_to_compare_what_string(exception_to_compare.what());
					return ( exception_to_compare_to_what_string == exception_to_compare_what_string );
				}
				catch(const std::bad_cast &) {
					return false;
				}
			}
		};

	} // namespace test
} // namespace cath

#endif
