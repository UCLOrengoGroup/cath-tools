/// \file
/// \brief The sort_uniq_build header

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
///  * the [...]_build functions provide convenient factories (and do so in such a way
///    as to hopefully avoid needless copy ctns in C++11 where move ctors can be used instead).
///  * the [...]uniq[...] functions perform not only the unique() operation to rearrange the elements
///    but also the erase to remove any duplicate leftover elements

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_SORT_UNIQ_BUILD_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_SORT_UNIQ_BUILD_HPP

#include "cath/common/algorithm/sort_uniq_copy.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"

namespace cath {
	namespace common {

		// prototypes to help anyone skimming this file
		template <typename Container, typename R            > Container        sort_build     ( R        );
		template <typename Container, typename R, typename P> Container        sort_build     ( R ,    P );
		template <typename Container, typename I            > Container        sort_build     ( I , I    );
		template <typename Container, typename I, typename P> Container        sort_build     ( I , I, P );
		template <typename Container, typename R            > Container stable_sort_build     ( R        );
		template <typename Container, typename R, typename P> Container stable_sort_build     ( R ,    P );
		template <typename Container, typename R            > Container        uniq_build     ( R        );
		template <typename Container, typename R            > Container        sort_uniq_build( R        );
		template <typename Container, typename R            > Container stable_sort_uniq_build( R        );

		/// \brief Convenience function for making a sorted copy of a range
		///
		/// \sa copy_build() generate_n_build() random_sample_n_build() sort_build() sort_uniq_build() transform_build() uniq_build()
		///
		/// Note that the _build suffix is used to mean something quite different in names like
		/// unique_build() in Boost and std (ie writing output to a specified output iterator).
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
		template <typename Container, typename R>
		Container sort_build(R prm_range ///< The range on which the sorted copy should be based
		                     ) {
			return sort_build<Container>(
				common::cbegin( prm_range ),
				common::cend  ( prm_range )
			);
		}

		/// \overload
		template <typename Container, typename R, typename P>
		Container sort_build(R prm_range,   ///< The range on which the sorted copy should be based
		                     P prm_bin_pred ///< The binary predicate to use as a less-than operator for sorting
		                     ) {
			return sort_build<Container>(
				common::cbegin( prm_range ),
				common::cend  ( prm_range ),
				prm_bin_pred
			);
		}

		template <typename Container, typename I>
		Container sort_build(I prm_begin, ///< The begin iterator for the range on which the sorted copy should be based
		                     I prm_end    ///< The begin iterator for the range on which the sorted copy should be based
		                     ) {
			return sort_copy(
				Container{ prm_begin, prm_end }
			);
		}

		/// \overload
		template <typename Container, typename I, typename P>
		Container sort_build(I prm_begin, ///< The begin iterator for the range on which the sorted copy should be based
		                     I prm_end,   ///< The begin iterator for the range on which the sorted copy should be based
		                     P prm_bin_pred ///< The binary predicate to use as a less-than operator for sorting
		                     ) {
			return sort_copy(
				Container{ prm_begin, prm_end },
				prm_bin_pred
			);
		}

		/// \brief Convenience function for making a sorted copy of a range
		///
		/// \sa copy_build() generate_n_build() random_sample_n_build() sort_build() sort_uniq_build() transform_build() uniq_build()
		///
		/// Note that the _build suffix is used to mean something quite different in names like
		/// unique_build() in Boost and std (ie writing output to a specified output iterator).
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
		template <typename Container, typename R>
		Container stable_sort_build(R prm_range ///< The range on which the sorted copy should be based
		                            ) {
			return stable_sort_copy(
				Container(
					common::cbegin( prm_range ),
					common::cend  ( prm_range )
				)
			);
		}

		/// \overload
		template <typename Container, typename R, typename P>
		Container stable_sort_build(R prm_range,   ///< The range on which the sorted copy should be based
		                            P prm_bin_pred ///< The binary predicate to use as a less-than operator for sorting
		                            ) {
			return stable_sort_copy(
				Container(
					common::cbegin( prm_range ),
					common::cend  ( prm_range )
				),
				prm_bin_pred
			);
		}

		/// \brief Convenience function for making a uniqued copy of a range
		///        (which, unlike the std/boost unique() functions, actually erases leftover elements)
		///
		/// \sa copy_build() generate_n_build() random_sample_n_build() sort_build() sort_uniq_build() transform_build() uniq_build()
		///
		/// This is very similar to sort_build() so its details are copied here...
		///
		/// \copydetails sort_build()
		template <typename Container, typename R>
		Container uniq_build(R prm_range ///< The range on which the uniqued copy should be based
		                     ) {
			return uniq_copy(
				Container(
					common::cbegin( prm_range ),
					common::cend  ( prm_range )
				)
			);
		}

		/// \brief Convenience function for making a sorted, uniqued copy of a range
		///
		/// \sa copy_build() generate_n_build() random_sample_n_build() sort_build() sort_uniq_build() transform_build() uniq_build()
		///
		/// \copydetails sort_build()
		template <typename Container, typename R>
		Container sort_uniq_build(R prm_range ///< The range on which the sorted, uniqued copy should be based
		                          ) {
			return sort_uniq_copy(
				Container(
					common::cbegin( prm_range ),
					common::cend  ( prm_range )
				)
			);
		}

		/// \brief Convenience function for making a sorted, uniqued copy of a range
		///
		/// \sa copy_build() generate_n_build() random_sample_n_build() sort_build() sort_uniq_build() transform_build() uniq_build()
		///
		/// \copydetails sort_build()
		template <typename Container, typename R, typename P>
		Container sort_uniq_build(R prm_range,   ///< The range on which the sorted, uniqued copy should be based
		                          P prm_bin_pred ///< The binary predicate to use as a less-than operator for sorting
		                          ) {
			return sort_uniq_copy(
				Container(
					common::cbegin( prm_range ),
					common::cend  ( prm_range )
				),
				prm_bin_pred
			);
		}


		/// \brief Convenience function for making a stable_sorted, uniqued copy of a range
		///
		/// \sa copy_build() generate_n_build() random_sample_n_build() sort_build() sort_uniq_build() transform_build() uniq_build()
		///
		/// \copydetails sort_build()
		template <typename Container, typename R>
		Container stable_sort_uniq_build(R prm_range ///< The range on which the stable_sorted, uniqued copy should be based
		                                 ) {
			return stable_sort_uniq_copy(
				Container(
					common::cbegin( prm_range ),
					common::cend  ( prm_range )
				)
			);
		}

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_SORT_UNIQ_BUILD_HPP
