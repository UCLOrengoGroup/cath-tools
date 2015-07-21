/// \file
/// \brief The batch_functions header

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

#ifndef BATCH_FUNCTIONS_H_INCLUDED
#define BATCH_FUNCTIONS_H_INCLUDED

#include "common/batch/broken_batch_tol.h"
#include "common/type_aliases.h"

namespace cath {
	namespace common {

		size_t batch_size(const size_t &,
		                  const size_t &,
		                  const broken_batch_tol &);

		size_t num_batches(const size_t &,
		                   const size_t &,
		                   const broken_batch_tol &);

		void check_batch_index(const size_t &,
		                       const size_t &,
		                       const size_t &,
		                       const broken_batch_tol &);

		size_t batch_begin(const size_t &,
		                   const size_t &,
		                   const size_t &,
		                   const broken_batch_tol &);

		size_t batch_end(const size_t &,
		                 const size_t &,
		                 const size_t &,
		                 const broken_batch_tol &);

		size_size_pair batch_begin_and_end(const size_t &,
		                                   const size_t &,
		                                   const size_t &,
		                                   const broken_batch_tol &);

		size_size_pair batch_start_and_stop(const size_t &,
		                                    const size_t &,
		                                    const size_t &,
		                                    const broken_batch_tol &);

	}
}

#endif
