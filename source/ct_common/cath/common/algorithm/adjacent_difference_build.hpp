/// \file
/// \brief The adjacent_difference_build() header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_ADJACENT_DIFFERENCE_BUILD_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_ADJACENT_DIFFERENCE_BUILD_HPP

#include <boost/range/numeric.hpp>

namespace cath::common {

	/// \brief Constructs, populates and return an output container from a range
	///
	/// WARNING: adjacent_difference() is a bit rubbish - it copies the first element of the input range into
	///          the output range before doing calculating the adjacent_differences. In practice, this means
	///          the output range must have the same value type and the input range.
	///
	/// \sa adjacent_difference_build() generate_n_build() random_sample_n_build() sort_adjacent_difference() sort_uniq_adjacent_difference() transform_build() uniq_adjacent_difference()
	///
	/// \pre Container        is a model of the Mutable_Container concept
	/// \pre SinglePassRange1 is a model of the SinglePassRangeConcept
	template <typename Container,
	          typename SinglePassRange1>
	inline Container adjacent_difference_build(const SinglePassRange1 &rng1 ///< A single-pass input range
	                                           ) {
		// Static-check that Container is a Mutable_Container
		BOOST_CONCEPT_ASSERT(( boost::Mutable_Container< Container > ));

		// Check that SinglePassRange is a SinglePassRangeConcept
		BOOST_RANGE_CONCEPT_ASSERT(( boost::SinglePassRangeConcept< const SinglePassRange1 > ));

		// Construct an instance of the container
		Container container;

		// Call the normal Boost Range adjacent_difference()
		boost::adjacent_difference(
			rng1,
			inserter( container, std::end( container ) )
		);

		// Return the populated container
		return container;
	}

	/// \brief Constructs, populates and return an output container from a range
	////
	/// WARNING: adjacent_difference() is a bit rubbish - it copies the first element of the input range into
	///          the output range before doing calculating the adjacent_differences. In practice, this means
	///          the output range must have the same value type and the input range.
	///
	/// \sa adjacent_difference_build() generate_n_build() random_sample_n_build() sort_adjacent_difference() sort_uniq_adjacent_difference() transform_build() uniq_adjacent_difference()
	///
	/// \pre Container        is a model of the Mutable_Container concept
	/// \pre SinglePassRange1 is a model of the SinglePassRangeConcept
	template <typename Container,
	          typename SinglePassRange1,
	          typename BinaryOperation>
	inline Container adjacent_difference_build(const SinglePassRange1 &rng1, ///< A single-pass input range
	                                           BinaryOperation         fun   ///< A binary function to execute on the pairwise list of elements of rng1 and rng2
	                                           ) {
		// Static-check that Container is a Mutable_Container
		BOOST_CONCEPT_ASSERT(( boost::Mutable_Container< Container > ));

		// Check that SinglePassRange is a SinglePassRangeConcept
		BOOST_RANGE_CONCEPT_ASSERT(( boost::SinglePassRangeConcept< const SinglePassRange1 > ));

		// Construct an instance of the container
		Container container;

		// Call the normal Boost Range adjacent_difference()
		boost::adjacent_difference(
			rng1,
			inserter( container, std::end( container ) ),
			fun
		);

		// Return the populated container
		return container;
	}

} // namespace cath::common
#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_ADJACENT_DIFFERENCE_BUILD_HPP
