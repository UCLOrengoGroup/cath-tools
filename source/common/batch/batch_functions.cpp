/// \file
/// \brief The batch definitions

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Binaries project and then tweaked, eg namespaced in cath)
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

#include "batch_functions.h"

#include <boost/lexical_cast.hpp>

#include "exception/invalid_argument_exception.h"

#include <algorithm>

using namespace cath;
using namespace cath::common;
using namespace std;

using boost::lexical_cast;

/// \brief TODOCUMENT
size_t cath::common::batch_size(const size_t           &arg_num_items,       ///< TODOCUMENT
                                const size_t           &arg_num_batches,     ///< TODOCUMENT
                                const broken_batch_tol &arg_broken_batch_tol ///< TODOCUMENT
                                ) {
	if ( arg_num_batches == 0 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot form zero batches of size zero"));
	}

	const size_t batch_size_if_not_broken = arg_num_items / arg_num_batches;
	if ( arg_num_items % arg_num_batches == 0 ) {
		return batch_size_if_not_broken;
	}
	if ( arg_broken_batch_tol == broken_batch_tol::FORBID ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Batching "      + lexical_cast<string>( arg_num_items   )
			+ " items into " + lexical_cast<string>( arg_num_batches )
			+ " batches leaves a broken batch, which has been forbidden"
		));
	}
	return batch_size_if_not_broken + 1;
}

/// \brief TODOCUMENT
size_t cath::common::num_batches(const size_t           &arg_num_items,       ///< TODOCUMENT
                                 const size_t           &arg_batch_size,      ///< TODOCUMENT
                                 const broken_batch_tol &arg_broken_batch_tol ///< TODOCUMENT
                                 ) {
	if ( arg_batch_size == 0 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot form batches of size zero"));
	}

	const size_t num_full_batches = arg_num_items / arg_batch_size;
	if ( arg_num_items % arg_batch_size == 0 ) {
		return num_full_batches;
	}
	if ( arg_broken_batch_tol == broken_batch_tol::FORBID ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Batching "                      + lexical_cast<string>( arg_num_items  )
			+ " items into batches of size " + lexical_cast<string>( arg_batch_size )
			+ " leaves a broken batch, which has been forbidden"
		));
	}
	return num_full_batches + 1;
}


/// \brief TODOCUMENT
void cath::common::check_batch_index(const size_t           &arg_num_items,       ///< TODOCUMENT
                                     const size_t           &arg_batch_size,      ///< TODOCUMENT
                                     const size_t           &arg_batch_index,     ///< TODOCUMENT
                                     const broken_batch_tol &arg_broken_batch_tol ///< TODOCUMENT
                                     ) {
	const size_t number_of_batches = num_batches( arg_num_items, arg_batch_size, arg_broken_batch_tol );
	if ( arg_batch_index >= number_of_batches ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Batch index "
			+ lexical_cast<string>( arg_batch_index   )
			+ " should be less than the number of batches ("
			+ lexical_cast<string>( number_of_batches )
			+ ")"
		));
	}
}

/// \brief TODOCUMENT
size_t cath::common::batch_begin(const size_t           &arg_num_items,       ///< TODOCUMENT
                                 const size_t           &arg_batch_size,      ///< TODOCUMENT
                                 const size_t           &arg_batch_index,     ///< TODOCUMENT
                                 const broken_batch_tol &arg_broken_batch_tol ///< TODOCUMENT
                                 ) {
	check_batch_index( arg_num_items, arg_batch_size, arg_batch_index, arg_broken_batch_tol );
	return arg_batch_index * arg_batch_size;
}

/// \brief TODOCUMENT
///
/// Note: The end value is one-past-the-last-value (like C++ iterators)
size_t cath::common::batch_end(const size_t           &arg_num_items,       ///< TODOCUMENT
                               const size_t           &arg_batch_size,      ///< TODOCUMENT
                               const size_t           &arg_batch_index,     ///< TODOCUMENT
                               const broken_batch_tol &arg_broken_batch_tol ///< TODOCUMENT
                               ) {
	const size_t begin = batch_begin( arg_num_items, arg_batch_size, arg_batch_index, arg_broken_batch_tol );
	return min( arg_num_items, begin + arg_batch_size );
}

/// \brief TODOCUMENT
///
/// Note: The end value is one-past-the-last-value (like C++ iterators)
size_size_pair cath::common::batch_begin_and_end(const size_t           &arg_num_items,       ///< TODOCUMENT
                                                 const size_t           &arg_batch_size,      ///< TODOCUMENT
                                                 const size_t           &arg_batch_index,     ///< TODOCUMENT
                                                 const broken_batch_tol &arg_broken_batch_tol ///< TODOCUMENT
                                                 ) {
	const size_t begin = batch_begin( arg_num_items, arg_batch_size, arg_batch_index, arg_broken_batch_tol );
	const size_t end   = batch_end  ( arg_num_items, arg_batch_size, arg_batch_index, arg_broken_batch_tol );
	assert( begin < end );
	return make_pair( begin, end );
}

/// \brief TODOCUMENT
size_size_pair cath::common::batch_start_and_stop(const size_t           &arg_num_items,       ///< TODOCUMENT
                                                  const size_t           &arg_batch_size,      ///< TODOCUMENT
                                                  const size_t           &arg_batch_index,     ///< TODOCUMENT
                                                  const broken_batch_tol &arg_broken_batch_tol ///< TODOCUMENT
                                                  ) {
	const size_size_pair begin_and_end = batch_begin_and_end( arg_num_items, arg_batch_size, arg_batch_index, arg_broken_batch_tol );
	return make_pair( begin_and_end.first, begin_and_end.second - 1 );
}
