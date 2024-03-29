/// \file
/// \brief The clone_ptr header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CLONE_CLONE_PTR_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CLONE_CLONE_PTR_HPP

#include <boost/serialization/nvp.hpp> // not needed for compilation but kept as is used in this header's templates
#include <boost/serialization/unique_ptr.hpp> // not needed for compilation but kept as is used in this header's templates

#include "cath/common/clone/detail/is_cloneable.hpp"

#include <cassert>
#include <memory>
#include <string>
#include <typeinfo>

namespace cath::common {

	/// \brief A light wrapper for std::unique_ptr that provides a copy ctor and assignment operator
	///        by using the fact that the template class provides a sensible clone method.
	///
	/// is_cloneable is used to check (at compile time) that the template class does indeed
	/// have a sensible clone method.
	template <typename cloneable_pointee_type>
	class clone_ptr final {
	private:
		friend class boost::serialization::access;
		template<typename archive> void serialize(archive & ar,
		                                          const size_t /*version*/
		                                          ) {
			ar & BOOST_SERIALIZATION_NVP( ptr );
		}

		BOOST_CONCEPT_ASSERT(( detail::is_cloneable<cloneable_pointee_type> ));

		std::unique_ptr<cloneable_pointee_type> ptr;

	public:
		/// \brief A type alias for the type of the object managed by this clone_ptr
		using element_type = cloneable_pointee_type;

		/// \brief Pass through to unique_ptr's ctor from a pointer
		explicit clone_ptr(cloneable_pointee_type * prm_pointer = nullptr ///< TODOCUMENT
		                   ) : ptr( prm_pointer ) {
		}

		/// \brief A copy constructor which uses the clone method of cloneable_pointee_type
		clone_ptr(const clone_ptr &prm_clone_ptr ///< TODOCUMENT
		          ) : ptr( detail::make_clone( *prm_clone_ptr ) ) {
			assert( typeid( *ptr.get() ) == typeid( *prm_clone_ptr.get() ) ); // Check the dynamic type returned by clone()
		}

		/// \brief A move-ctor that uses the default behaviour (ie the move-ctor of unique_ptr)
		clone_ptr(clone_ptr &&) noexcept = default;

		/// \brief Ctor from a unique_ptr const lvalue reference, which takes a clone from it
		explicit clone_ptr(const std::unique_ptr<cloneable_pointee_type> &prm_unique_ptr ///< TODOCUMENT
		                   ) : ptr( detail::make_clone( *prm_unique_ptr ) ) {
		}

		/// \brief Ctor from a unique_ptr rvalue reference, which moves it
		explicit clone_ptr(std::unique_ptr<cloneable_pointee_type> &&prm_unique_ptr ///< TODOCUMENT
		                   ) noexcept : ptr( std::move( prm_unique_ptr ) ) {
		}

		/// \brief An assignment operator which uses the clone method of cloneable_pointee_type
		clone_ptr & operator=(const clone_ptr &prm_clone_ptr ///< TODOCUMENT
		                      ) {
			// Check for self-assignment
			if ( this != &prm_clone_ptr ) {
				std::unique_ptr<cloneable_pointee_type> temp_ptr( prm_clone_ptr->clone() );
				*this = temp_ptr;
				assert( typeid( *ptr.get() ) == typeid( *prm_clone_ptr.get() ) ); // Check the dynamic type returned by clone()
			}
			return *this;
		}

		/// \brief A move assignment operator that uses the default behaviour (ie the move assignment operator of unique_ptr)
		clone_ptr & operator=(clone_ptr &&) noexcept = default;

		/// \brief A copy assignment operator from a unique_ptr lvalue reference
		clone_ptr & operator=(const std::unique_ptr<cloneable_pointee_type> &prm_unique_ptr ///< TODOCUMENT
		                      ) {
			// Check for self-assignment
			if ( &ptr != &prm_unique_ptr ) {
				ptr = prm_unique_ptr->clone();
				assert( typeid( *ptr ) == typeid( *prm_unique_ptr ) ); // Check the dynamic type returned by clone()
			}
			return *this;
		}

		/// \brief A move assignment operator from a unique_ptr rvalue reference
		clone_ptr & operator=(std::unique_ptr<cloneable_pointee_type> &&prm_unique_ptr ///< TODOCUMENT
		                      ) noexcept {
			ptr = std::move( prm_unique_ptr );
			return *this;
		}

		/// \brief Pass through to std::unique_ptr's reset method
		void reset(cloneable_pointee_type *prm_pointer = nullptr ///< TODOCUMENT
		           ) {
			ptr.reset(prm_pointer);
		}
		/// \brief Pass through to std::unique_ptr's operator* method
		cloneable_pointee_type & operator*() const {
			return ptr.operator*();
		}
		/// \brief Pass through to std::unique_ptr's operator-> method
		cloneable_pointee_type * operator->() const {
			return ptr.operator->();
		}
		/// \brief Pass through to std::unique_ptr's get method
		cloneable_pointee_type * get() const {
			return ptr.get();
		}
		/// \brief TODOCUMENT
		explicit operator bool() const {
			return static_cast<bool>( ptr );
		}
		/// \brief Pass through to std::unique_ptr's swap method
		void swap(clone_ptr &prm_clone_ptr ///< TODOCUMENT
		          ) {
			ptr.swap(prm_clone_ptr.ptr);
		}
	};

	template <typename T, typename U>
	inline bool operator==(const clone_ptr<T> &prm_clone_ptr_a, ///< TODOCUMENT
	                       const clone_ptr<T> &prm_clone_ptr_b  ///< TODOCUMENT
	                       ) {
		return prm_clone_ptr_a.get() == prm_clone_ptr_b.get();
	}
} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CLONE_CLONE_PTR_HPP
