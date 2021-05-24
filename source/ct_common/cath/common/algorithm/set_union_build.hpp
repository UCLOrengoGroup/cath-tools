/// \file
/// \brief The set_union_build() header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_SET_UNION_BUILD_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_SET_UNION_BUILD_HPP

#include <boost/range/algorithm/set_algorithm.hpp>

namespace cath::common {

	/// \brief Constructs, populates and return an output container from a range
	///
	/// WARNING: set_union() is a bit rubbish - it copies the first element of the input range into
	///          the output range before doing calculating the set_unions. In practice, this means
	///          the output range must have the same value type and the input range.
	///
	/// \sa set_union_build() generate_n_build() random_sample_n_build() sort_set_union() sort_uniq_set_union() transform_build() uniq_set_union()
	///
	/// \pre Container        is a model of the Mutable_Container concept
	/// \pre SinglePassRange1 is a model of the SinglePassRangeConcept
	template <typename Container,
	          typename SinglePassRange1,
	          typename SinglePassRange2>
	inline Container set_union_build(const SinglePassRange1 &rng1, ///< A single-pass input range
	                                 const SinglePassRange2 &rng2  ///< A single-pass input range
	                                 ) {
		// Static-check that Container is a Mutable_Container
		BOOST_CONCEPT_ASSERT(( boost::Mutable_Container< Container > ));

		// Check that SinglePassRange1 and SinglePassRange2 meet SinglePassRangeConcept
		BOOST_RANGE_CONCEPT_ASSERT(( boost::SinglePassRangeConcept< const SinglePassRange1 > ));
		BOOST_RANGE_CONCEPT_ASSERT(( boost::SinglePassRangeConcept< const SinglePassRange2 > ));

		// Construct an instance of the container
		Container container;

		// Call the normal Boost Range set_union()
		boost::range::set_union(
			rng1,
			rng2,
			inserter( container, std::end( container ) )
		);

		// Return the populated container
		return container;
	}

	/// \brief Constructs, populates and return an output container from a range
	////
	/// WARNING: set_union() is a bit rubbish - it copies the first element of the input range into
	///          the output range before doing calculating the set_unions. In practice, this means
	///          the output range must have the same value type and the input range.
	///
	/// \sa set_union_build() generate_n_build() random_sample_n_build() sort_set_union() sort_uniq_set_union() transform_build() uniq_set_union()
	///
	/// \pre Container        is a model of the Mutable_Container concept
	/// \pre SinglePassRange1 is a model of the SinglePassRangeConcept
	template <typename Container,
	          typename SinglePassRange1,
	          typename SinglePassRange2,
	          typename BinaryOperation>
	inline Container set_union_build(const SinglePassRange1 &rng1, ///< A single-pass input range
	                                 const SinglePassRange2 &rng2, ///< A single-pass input range
	                                 BinaryOperation         fun   ///< A binary function to execute on the pairwise list of elements of rng1 and rng2
	                                 ) {
		// Static-check that Container is a Mutable_Container
		BOOST_CONCEPT_ASSERT(( boost::Mutable_Container< Container > ));

		// Check that SinglePassRange1 and SinglePassRange2 meet SinglePassRangeConcept
		BOOST_RANGE_CONCEPT_ASSERT(( boost::SinglePassRangeConcept< const SinglePassRange1 > ));
		BOOST_RANGE_CONCEPT_ASSERT(( boost::SinglePassRangeConcept< const SinglePassRange2 > ));

		// Construct an instance of the container
		Container container;

		// Call the normal Boost Range set_union()
		boost::range::set_union(
			rng1,
			rng2,
			inserter( container, std::end( container ) ),
			fun
		);

		// Return the populated container
		return container;
	}

} // namespace cath::common
#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_SET_UNION_BUILD_HPP
