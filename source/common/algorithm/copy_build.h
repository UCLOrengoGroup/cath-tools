/// \file
/// \brief The copy_build() header

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

#ifndef COPY_BUILD_H_INCLUDED
#define COPY_BUILD_H_INCLUDED

#include <boost/range/algorithm/copy.hpp>

namespace cath {
	namespace common {

		/// \brief Constructs, populates and return an output container from a range
		///
		/// \sa copy_build() generate_n_build() random_sample_n_build() sort_copy() sort_uniq_copy() transform_build() uniq_copy()
		///
		/// \tparam Container is a model of the Mutable_Container concept
		/// \tparam ITER      is TODOCUMENT
		template <typename Container,
		          typename ITER>
		inline Container copy_build(const ITER &arg_begin_itr, ///< TODOCUMENT
		                            const ITER &arg_end_itr    ///< TODOCUMENT
		                            ) {
			// Static-check that Container is a Mutable_Container
			BOOST_CONCEPT_ASSERT(( boost::Mutable_Container< Container > ));

			// Construct an instance of the container
			Container container;

			// Call the normal std::copy
			std::copy(
				arg_begin_itr,
				arg_end_itr,
				inserter( container, std::end( container ) )
			);

			// Return the populated container
			return container;
		}

		/// \brief Constructs, populates and return an output container from a range
		///
		/// \sa copy_build() generate_n_build() random_sample_n_build() sort_copy() sort_uniq_copy() transform_build() uniq_copy()
		///
		/// \tparam Container        is a model of the Mutable_Container concept
		/// \tparam SinglePassRange1 is a model of the SinglePassRangeConcept
		template <typename Container,
		          typename SinglePassRange1>
		inline Container copy_build(const SinglePassRange1 &rng1 ///< A single-pass input range
		                            ) {
			// Static-check that Container is a Mutable_Container
			BOOST_CONCEPT_ASSERT(( boost::Mutable_Container< Container > ));

			// Check that SinglePassRange is a SinglePassRangeConcept
			BOOST_RANGE_CONCEPT_ASSERT(( boost::SinglePassRangeConcept< SinglePassRange1 > ));

			// Construct an instance of the container
			Container container;

			// Call the normal Boost Range copy()
			boost::range::copy(
				rng1,
				inserter( container, std::end( container ) )
			);

			// Return the populated container
			return container;
		}

	}
}
#endif
