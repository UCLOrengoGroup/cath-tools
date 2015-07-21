/// \file
/// \brief The superposition_acquirer class header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 2011, Orengo Group, University College London
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

#ifndef SUPERPOSITION_ACQUIRER_H_INCLUDED
#define SUPERPOSITION_ACQUIRER_H_INCLUDED

#include <boost/tuple/tuple.hpp>

#include "common/type_aliases.h"

#include <iostream>

namespace cath { namespace opts { class cath_superpose_options; } }
namespace cath { namespace sup { class superposition_context; } }

namespace cath {
	namespace opts {

		/// \brief TODOCUMENT
		class superposition_acquirer {
		private:
			virtual sup::superposition_context do_get_superposition(std::ostream &) const = 0;

		public:
			virtual ~superposition_acquirer() noexcept = default;

			sup::superposition_context get_superposition(std::ostream &) const;

			static const double PERCENT_TOLERANCE_FOR_EQUAL_RMSDS;
		};

	}
}

#endif
