/// \file
/// \brief The sort_uniq_copy header

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
///
/// Note that sort() isn't required because that's already provided in boost range.
///
/// This file contains convenient wrapper functions for the sort() and unique() functions
/// in Boost range. These add two things:
///  * the [...]_copy functions provide convenient factories (and do so in such a way
///    as to hopefully avoid needless copy ctns in C++11 where move ctors can be used instead).
///  * the [...]uniq[...] functions perform not only the unique() operation to rearrange the elements
///    but also the erase to remove any duplicate leftover elements

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_SORT_UNIQ_COPY_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_SORT_UNIQ_COPY_HPP

#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext/erase.hpp>

namespace cath::common {

	// prototypes to help anyone skimming this file
	template <typename R            > R           sort_copy      ( R      );
	template <typename R, typename P> R           sort_copy      ( R,   P );

	template <typename R            > R    stable_sort_copy      ( R      );
	template <typename R, typename P> R    stable_sort_copy      ( R,   P );

	template <typename R            > void        uniq           ( R &    );

	template <typename R            > R           uniq_copy      ( R      );

	template <typename R            > void        sort_uniq      ( R &    );

	template <typename R, typename P> void        sort_uniq      ( R &, P );

	template <typename R            > void stable_sort_uniq      ( R &    );

	template <typename R            > R           sort_uniq_copy ( R      );

	template <typename R, typename P> R           sort_uniq_copy ( R,   P );

	template <typename R            > R    stable_sort_uniq_copy ( R      );

	/// \brief Convenience function for making a sorted copy of a range
	///
	/// \sa copy_build() generate_n_build() random_sample_n_build() sort_copy() sort_uniq_copy() transform_build() uniq_copy()
	///
	/// Note that the _copy suffix is used to mean something quite different in names like
	/// unique_copy() in Boost and std (ie writing output to a specified output iterator).
	///
	/// Come C++11, this should be able to avoid making a copy when its passed an rvalue argument.
	/// This is achieved by the argument being passed by non-const value
	/// (rather than by const reference) so that the compiler handles making
	/// the local, modifiable copy and can use a move ctor where appropriate.
	///
	/// \todo Write a test to check that these functions don't call needlessly
	///       call copy ctors. Do this with a generalised test class that throws on attempts
	///       to call its copy ctor. (Perhaps: `template <typename T> class copy_ctor_thrower final : public T { ... };` )
	///
	/// ATM, this doesn't perform any concept checks and leaves that to Boost Range's sort().
	template <typename R>
	R sort_copy(R prm_range ///< The range on which the sorted copy should be based
	            ) {
		boost::range::sort( prm_range );
		return prm_range;
	}

	/// \overload
	template <typename R, typename P>
	R sort_copy(R prm_range,   ///< The range on which the sorted copy should be based
	            P prm_bin_pred ///< The binary predicate to use as a less-than operator for sorting
	            ) {
		boost::range::sort( prm_range, prm_bin_pred );
		return prm_range;
	}

	/// \brief Convenience function for making a stable_sorted copy of a range
	///
	/// \sa copy_build() generate_n_build() random_sample_n_build() sort_copy() stable_sort_uniq_copy() transform_build() uniq_copy()
	///
	/// Note that the _copy suffix is used to mean something quite different in names like
	/// unique_copy() in Boost and std (ie writing output to a specified output iterator).
	///
	/// Come C++11, this should be able to avoid making a copy when its passed an rvalue argument.
	/// This is achieved by the argument being passed by non-const value
	/// (rather than by const reference) so that the compiler handles making
	/// the local, modifiable copy and can use a move ctor where appropriate.
	///
	/// \todo Write a test to check that these functions don't call needlessly
	///       call copy ctors. Do this with a generalised test class that throws on attempts
	///       to call its copy ctor. (Perhaps: `template <typename T> class copy_ctor_thrower final : public T { ... };` )
	///
	/// ATM, this doesn't perform any concept checks and leaves that to Boost Range's stable_sort().
	template <typename R>
	R stable_sort_copy(R prm_range ///< The range on which the stable_sorted copy should be based
	                   ) {
		boost::range::stable_sort( prm_range );
		return prm_range;
	}

	/// \overload
	template <typename R, typename P>
	R stable_sort_copy(R prm_range,   ///< The range on which the stable_sorted copy should be based
	                   P prm_bin_pred ///< The binary predicate to use as a less-than operator for stable_sorting
	                   ) {
		boost::range::stable_sort( prm_range, prm_bin_pred );
		return prm_range;
	}

	/// \brief Convenience function for uniquing a range and
	///        (unlike the std/boost unique() functions) erasing leftover elements
	template <typename R>
	void uniq(R &prm_range ///< The range to unique(-erase)
	          ) {
		boost::range::erase( prm_range, boost::range::unique<boost::return_found_end>( prm_range ) );
	}

	/// \brief Convenience function for making a uniqued copy of a range
	///        (which, unlike the std/boost unique() functions, actually erases leftover elements)
	///
	/// \sa copy_build() generate_n_build() random_sample_n_build() sort_copy() sort_uniq_copy() transform_build() uniq_copy()
	///
	/// This is very similar to sort_copy() so its details are copied here...
	///
	/// \copydetails sort_copy()
	template <typename R>
	R uniq_copy(R prm_range ///< The range on which the uniqued copy should be based
	            ) {
		uniq( prm_range );
		return prm_range;
	}

	/// \brief Convenience function for sorting+uniquing a range and
	///        (unlike the std/boost unique() functions) erasing leftover elements
	template <typename R>
	void sort_uniq(R &prm_range ///< The range to sort and unique(-erase)
	               ) {
		boost::range::sort( prm_range );
		uniq( prm_range );
	}

	/// \brief Convenience function for sorting+uniquing a range and
	///        (unlike the std/boost unique() functions) erasing leftover elements
	template <typename R, typename P>
	void sort_uniq(R &prm_range,  ///< The range to sort and unique(-erase)
	               P prm_bin_pred ///< The binary predicate to use as a less-than operator for sorting
	               ) {
		boost::range::sort( prm_range, prm_bin_pred );
		uniq( prm_range );
	}

	/// \brief Convenience function for stable_sorting+uniquing a range and
	///        (unlike the std/boost unique() functions) erasing leftover elements
	template <typename R>
	void stable_sort_uniq(R &prm_range ///< The range to stable_sort and unique(-erase)
	                      ) {
		boost::range::stable_sort( prm_range );
		uniq( prm_range );
	}

	/// \brief Convenience function for making a sorted, uniqued copy of a range
	///
	/// \sa copy_build() generate_n_build() random_sample_n_build() sort_copy() sort_uniq_copy() transform_build() uniq_copy()
	///
	/// \copydetails sort_copy()
	template <typename R>
	R sort_uniq_copy(R prm_range ///< The range on which the sorted, uniqued copy should be based
	                 ) {
		sort_uniq( prm_range );
		return prm_range;
	}

	/// \overload
	template <typename R, typename P>
	R sort_uniq_copy(R prm_range,   ///< The range on which the sorted, uniqued copy should be based
	                 P prm_bin_pred ///< The binary predicate to use as a less-than operator for sorting
	                 ) {
		sort_uniq( prm_range, prm_bin_pred );
		return prm_range;
	}

	/// \brief Convenience function for making a stable_sorted, uniqued copy of a range
	///
	/// \sa copy_build() generate_n_build() random_sample_n_build() stable_sort_copy() stable_sort_uniq_copy() transform_build() uniq_copy()
	///
	/// \copydetails stable_sort_copy()
	template <typename R>
	R stable_sort_uniq_copy(R prm_range ///< The range on which the stable_sorted, uniqued copy should be based
	                        ) {
		stable_sort_uniq( prm_range );
		return prm_range;
	}

} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_SORT_UNIQ_COPY_HPP
