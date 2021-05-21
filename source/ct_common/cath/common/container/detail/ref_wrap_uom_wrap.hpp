/// \file
/// \brief The ref_wrap_uom_wrap header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CONTAINER_DETAIL_REF_WRAP_UOM_WRAP_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CONTAINER_DETAIL_REF_WRAP_UOM_WRAP_HPP

#include <functional>

namespace cath {
	namespace common {
		namespace detail {

			/// \brief Wrap a reference_wrapper<T> for use in an unordered_map with hashing on the T value
			///
			/// Equality is also performed on the T value.
			///
			/// \todo If this continues to be used, move it to its own header
			template <typename T>
			class ref_wrap_uom_wrap final {
			private:
				/// \brief The reference_wrapper of the value
				std::reference_wrapper<T> value;

			public:
				/// \brief Ctor from a reference
				explicit ref_wrap_uom_wrap(T &prm_value ///< The value from which to construct
				                  ) : value{ prm_value } {
				}

				/// \brief Prevent construction from an rvalue
				ref_wrap_uom_wrap(T&& x ) = delete;

				/// \brief Ctor from a reference_wrapper<T>
				explicit ref_wrap_uom_wrap(const std::reference_wrapper<T> &prm_value ///< The reference_wrapper<T> from which this should be constructed
				                           ) : value{ prm_value } {
				}
				/// \brief Copy ctor
				ref_wrap_uom_wrap(const ref_wrap_uom_wrap<T> &prm_value ///< The ref_wrap_uom_wrap from which this should be constructed
				                  ) : value{ prm_value.get_ref_wrap() } {
				}

				/// \brief Assignment from a reference_wrapper<T>
				ref_wrap_uom_wrap & operator=(const std::reference_wrapper<T> &prm_rhs ///< The reference_wrapper<T> to assign
				                              ) {
					value = prm_rhs;
					return *this;
				}
				/// \brief Copy assignment operator
				ref_wrap_uom_wrap & operator=(const ref_wrap_uom_wrap<T> &prm_rhs ///< The ref_wrap_uom_wrap to assign
				                              ) {
					value = prm_rhs.get_ref_wrap();
					return *this;
				}

				/// \brief Get the reference_wrapper
				const std::reference_wrapper<T> & get_ref_wrap() {
					return value;
				}

				/// \brief Conversion to a reference to T
				explicit operator T & () const {
					return value.get();
				}

				/// \brief Get a reference to T
				T & get() const {
					return value.get();
				}

				/// \brief Invoke the wrapped value on the specified arguments
				template  <typename... Ts>
				std::result_of_t< T & (Ts &&...) > operator()(Ts &&... args ///< The values to pass to the value
				                                              ) const {
					return value( std::forward<Ts>( args )... );
				}
			};

			/// \brief Equality operator for ref_wrap_uom_wrap
			template <typename T>
			inline bool operator==(const ref_wrap_uom_wrap<T> &prm_lhs, ///< The first  ref_wrap_uom_wrap to compare
			                       const ref_wrap_uom_wrap<T> &prm_rhs  ///< The second ref_wrap_uom_wrap to compare
			                       ) {
				return ( prm_lhs.get() == prm_rhs.get() );
			}
		} // namespace detail
	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CONTAINER_DETAIL_REF_WRAP_UOM_WRAP_HPP
