/// \file
/// \brief The generate_n_build() header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_GENERATE_N_BUILD_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_GENERATE_N_BUILD_HPP

#include <boost/concept_check.hpp>
#include <boost/range/algorithm.hpp>

namespace cath {
	namespace common {

		/// \brief Factory wrapper for std::generate_n() that constructs, populates and returns the output container
		///
		/// \sa copy_build() generate_n_build() gcc_random_sample_n_build() sort_copy() sort_uniq_copy() transform_build() uniq_copy()
		///
		/// Motivation: to help simplify client code by removing one more obstacle to it being a list of
		/// assignments to const variables.
		///
		/// If you're not already familiar with std::generate_n(), you should
		/// probably read about the basic idea before reading this.
		///
		/// Example usage (using Boost Lambda):
		///
		/// \code
		/// const doub_vec ten_rands = generate_n_build<doub_vec>(
		///     10,
		///     [] { return rand(); }
		/// );
		///
		/// Pre-C++11, this will invoke at least one inefficient call to the container's copy-ctor to return it
		/// but any such calls should be replaced with efficient calls to its move-ctor in C++11.
		///
		/// This populates the container by using the output iterator returned by `inserter( container, std::end( container ) )`.
		/// This should work for associative containers whilst also being efficient in sequences, such as vector (because the
		/// insertion is always happening at the end).
		///
		/// At the moment, there isn't any attempt to call the output container's reserve() method (if it has one) so this code will
		/// potentially invoke multiple allocations in an output vector, even if the input's also as simple as a vector.
		///
		/// \todo Get this to call container.reserve( prm_size ) if:
		///  * Container has a reserve() method (create a new has_reserve_method concept) and
		///
		/// \pre Container        is a model of the Mutable_Container concept
		/// \pre Generator        is a model of the Generator that returns Container::value_type
		template <typename Container,
		          typename Size,
		          typename Gen>
		inline Container generate_n_build(const Size &prm_size, ///< The number of elements to generate
		                                  Gen         prm_fun   ///< A nullary function to create each of the elements
		                                  ) {
			using value_type = typename Container::value_type;
			// Static-check that Container is a Mutable_Container
			BOOST_CONCEPT_ASSERT(( boost::Mutable_Container< Container > ));

			// Check that SinglePassRange is a SinglePassRangeConcept
			BOOST_CONCEPT_ASSERT(( boost::Generator< Gen, value_type > ));

			// Construct an instance of the container
			Container container;

			// Call std::generate_n()
			std::generate_n(
				inserter( container, std::end( container ) ),
				prm_size,
				prm_fun
			);

			// Return the populated container
			return container;
		}

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_GENERATE_N_BUILD_HPP
