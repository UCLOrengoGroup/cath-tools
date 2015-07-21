/// \file
/// \brief The spanning_tree_greater class header

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

#ifndef SPANNING_TREE_GREATER_H_INCLUDED
#define SPANNING_TREE_GREATER_H_INCLUDED

#include "common/type_aliases.h"

namespace cath {
	namespace sup {
		class superpose_orderer;

		namespace detail {

			/// \brief Provide greater-than function operator for sorting pairs of indices by their score in a superpose_orderer
			class spanning_tree_greater final {
			private:
				const superpose_orderer &the_superpose_orderer;

			public:
				spanning_tree_greater(const superpose_orderer &);

				bool operator()(const size_size_pair &,
				                const size_size_pair &) const;
			};

		}
	}
}

#endif
