/// \file
/// \brief The transform_build() header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_ALGORITHM_TRANSFORM_BUILD_H
#define _CATH_TOOLS_SOURCE_COMMON_ALGORITHM_TRANSFORM_BUILD_H

#include <boost/range/algorithm.hpp>

#include "common/cpp14/cbegin_cend.hpp"

namespace cath {
	namespace common {

		/// \brief Factory wrapper for unary transform() that constructs, populates and returns the output container
		///
		/// \sa copy_build() generate_n_build() random_sample_n_build() sort_copy() sort_uniq_copy() transform_build() uniq_copy()
		///
		/// \todo Consider adding a version of each that also sorts the created container before returning it?
		///       Is there a way to avoid the combinatorial explosion of appending sort (with and without
		///       less-than-predicate) and other mutating algorithms to all three versions of transform_build
		///       (range (see above todo, range+unary, range+range+binary)?
		///       Perhaps introduce mutating algorithm as a function to be applied to the container
		///       so that sort (with or without less-than-predicate) can be conveniently applied via
		///       a C++11 lambda. A default value of an identity function could probably be very efficient.
		///       On second thoughts, maybe this shouldn't be complicated with that. Alternatively,
		///       the caller can just pass the output of the transform_build() through a lambda function
		///       (which shouldn't need any copies, come C++11). Eg:
		///           const size_vec a = [](auto &&x) { sort(x); return x; }( transform_build( ... ) );
		///       (not sure about rvalue ref there).
		///       It might be clearer to create a modify_rvalue() template:
		///           const size_vec a = modify_rvalue(
		///             transform_build( ... ),
		///             [](auto &x) { sort(x); }
		///           );
		///       ...coming back to this comment after a while, I can't see what's wrong with just
		///       passing the result through sort_copy(), which shouldn't require any copies. Eg:
		///       ~~~~~.cpp
		///       const size_vec a = sort_copy( transform_build( ... ) );
		///       ~~~~~
		///
		/// Motivation: to help simplify client code by removing one more obstacle to it being a list of
		/// assignments to const variables.
		///
		/// These functions wrap their Boost Range counterparts, which in turn wrap their std counterparts.
		/// If you're not already familiar with std::transform() (and boost::range::transform()), you should
		/// probably read about the basic idea of transform() before reading this.
		///
		/// Example usage (using Boost Lambda):
		///
		/// \code
		/// const size_vec three_four_five  = { 3, 4, 5 };
		///
		/// // Use unary form of transform_build() to make vector<size_t> by
		/// // adding four to each member of three_four_five
		/// const size_vec seven_eight_nine = transform_build<size_vec>(
		///    three_four_five,
		///    _1 + 4
		/// );
		///
		/// // Use binary form of transform_build() to make set<size_t> by
		/// // pairwise subtracting each member of three_four_five from seven_eight_nine
		/// const size_set four_four_four   = transform_build<size_set>(
		///    seven_eight_nine,
		///    three_four_five,
		///    _1 - _2
		/// );
		/// \endcode
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
		/// \todo Get this to call container.reserve( size( rng1 ) ) if:
		///  * Container has a reserve() method (create a new has_reserve_method concept) and
		///  * RandomAccessRangeConcept<SinglePassRange1> (so that size can be calculated in constant time)
		///
		/// Since these are implemented using the Boost Range transform(), the binary version stops on hitting the end of rng1 or rng2,
		/// and so cannot overrun the end of rng2, unlike the std::transform().
		///
		/// These functions don't concept-check the UnaryOperation/BinaryOperation but then nor do the Boost Range transform() implementations.
		///
		/// \pre Container        is a model of the Mutable_Container concept
		/// \pre SinglePassRange1 is a model of the SinglePassRangeConcept
		/// \pre SinglePassRange2 is a model of the SinglePassRangeConcept ( binary version only )
		/// \pre UnaryOperation   is a model of the UnaryFunctionConcept   ( unary  version only )
		/// \pre BinaryOperation  is a model of the BinaryFunctionConcept  ( binary version only )
		template <typename Container,
		          typename SinglePassRange1,
		          typename UnaryOperation>
		inline Container transform_build(const SinglePassRange1 &rng1, ///< A single-pass input range
		                                 UnaryOperation          fun  ///< A unary function to execute on each of the elements of rng
		                                 ) {
			// Static-check that Container is a Container
			BOOST_CONCEPT_ASSERT(( boost::Container< Container > ));
			/// \todo It seems as though this should test for Mutable_Container but that fails
			///       for maps (eg size_size_map). Does that indicate that Mutable_Container
			///       isn't the suitable concept here or is something else going on?

			// Check that SinglePassRange is a SinglePassRangeConcept
//			BOOST_RANGE_CONCEPT_ASSERT(( boost::SinglePassRangeConcept< SinglePassRange1 > ));

			// Construct an instance of the container
			Container container;

			// Call the normal Boost Range transform()
			boost::range::transform(
				rng1,
				inserter( container, std::end( container ) ),
				fun
			);

			// Return the populated container
			return container;
		}

		/// \brief Factory wrapper for binary transform() that constructs, populates and returns the output container
		///
		/// \copydetails transform_build()
		template <typename Container,
		          typename SinglePassRange1,
		          typename SinglePassRange2,
		          typename BinaryOperation>
		inline Container transform_build(const SinglePassRange1 &rng1, ///< A single-pass input range
		                                 const SinglePassRange2 &rng2, ///< A second single-pass input range
		                                 BinaryOperation         fun   ///< A binary function to execute on the pairwise list of elements of rng1 and rng2
		                                 ) {
			// Static-check that Container is a Mutable_Container
			BOOST_CONCEPT_ASSERT(( boost::Mutable_Container< Container > ));

			// Static-check that SinglePassRange is a SinglePassRangeConcept
//			BOOST_RANGE_CONCEPT_ASSERT(( boost::SinglePassRangeConcept< const SinglePassRange1 > ));
//			BOOST_RANGE_CONCEPT_ASSERT(( boost::SinglePassRangeConcept< const SinglePassRange2 > ));

			// Construct an instance of the container
			Container container;

			// Call the normal Boost Range transform()
			boost::range::transform(
				rng1,
				rng2,
				inserter( container, std::end( container ) ),
				fun
			);

			// Return the populated container
			return container;
		}
	} // namespace common

} // namespace cath
#endif
