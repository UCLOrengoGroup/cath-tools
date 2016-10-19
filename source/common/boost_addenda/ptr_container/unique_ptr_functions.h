/// \file
/// \brief The unique_ptr functions header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_BOOST_ADDENDA_PTR_CONTAINER_UNIQUE_PTR_FUNCTIONS_H
#define _CATH_TOOLS_SOURCE_COMMON_BOOST_ADDENDA_PTR_CONTAINER_UNIQUE_PTR_FUNCTIONS_H

#include <boost/ptr_container/ptr_map.hpp>
#include <boost/ptr_container/ptr_set.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include <memory>

namespace cath {
namespace common {

		/// \brief Provide a common way to insert a unique_ptr into a ptr_map
		///
		/// This is required because ptr_container doesn't yet support unique_ptr.
		///
		/// This only takes the key parameter by non-const value because ptr_map's
		/// pointer-based insert() inexplicably takes the first argument by non-const
		/// reference.
		///
		/// \todo If ptr_map gets fixed then change the key argument to a const reference
		template <typename K, typename V>
		void insert(boost::ptr_map<K, V>  &arg_ptr_map,  ///< The ptr_map into which the key/value should be inserted
		            K                      arg_key,      ///< The key to insert
		            std::unique_ptr<V>   &&arg_ptr_value ///< The associated value to insert, held by unique_ptr
		            ) {
			arg_ptr_map.insert( arg_key, arg_ptr_value.release() );
		}

		/// \brief Provide a common way to insert a unique_ptr into a ptr_set
		///
		/// This is required because ptr_container doesn't yet support unique_ptr.
		template <typename T>
		void insert(boost::ptr_set<T>   &arg_ptr_set,  ///< The ptr_set into which the value should be inserted
		            std::unique_ptr<T> &&arg_ptr_value ///< The associated value to insert, held by unique_ptr
		            ) {
			arg_ptr_set.insert( arg_ptr_value.release() );
		}

		/// \brief Provide a common way to push a unique_ptr onto the back of a ptr_vector
		///
		/// This is required because ptr_container doesn't yet support unique_ptr.
		template <typename T>
		void push_back(boost::ptr_vector<T>  &arg_ptr_vector, ///< The ptr_vector onto which the value should be pushed
		               std::unique_ptr<T>   &&arg_ptr_value   ///< The value to push_back, held by unique_ptr
		               ) {
			arg_ptr_vector.push_back( arg_ptr_value.release() );
		}

		/// \brief Provide a common way to replace a value in a ptr_vector from a unique_ptr
		///
		/// This is required because ptr_container doesn't yet support unique_ptr.
		template <typename T>
		void replace(boost::ptr_vector<T>  &arg_ptr_vector, ///< The ptr_vector onto which the value should be pushed
		             const size_t          &idx,            ///< The index specifying the position to be replaced
					 std::unique_ptr<T>   &&arg_ptr_value   ///< The value to push_back, held by unique_ptr
					 ) {
			arg_ptr_vector.replace( idx, arg_ptr_value.release() );
		}

	} // namespace common
} // namespace cath

#endif
