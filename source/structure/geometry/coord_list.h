/// \file
/// \brief The coord_list class header

/// \copyright
/// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef COORD_LIST_H_INCLUDED
#define COORD_LIST_H_INCLUDED

#include <boost/operators.hpp>
#include <boost/range.hpp>

#include "structure/structure_type_aliases.h"

#include <vector>

namespace cath {
	namespace geom {

		/// \brief TODOCUMENT
		class coord_list final : boost::additive<coord_list, coord> {
		private:
			/// \brief TODOCUMENT
			coord_vec coords;

		public:
			coord_list() = default;
			explicit coord_list(const coord_vec &);

			void reserve(const size_t &);
			bool empty() const noexcept;
			size_t size() const;
			void push_back(const coord &);
			coord & operator[](const size_t &);
			const coord & operator[](const size_t &) const;

			void operator+=(const coord &);
			void operator-=(const coord &);

			// Provide iterators to make this into a range
			using iterator = coord_vec::iterator;
			using const_iterator = coord_vec::const_iterator;
			iterator begin();
			iterator end();
			const_iterator begin() const;
			const_iterator end() const;
		};

		coord_list flatten_coord_lists(const coord_list_vec &);

		coord sum(const coord_list &);

		size_t check_non_empty_and_equal_size(const coord_list &,
		                                      const coord_list &);

		coord centre_of_gravity(const coord_list &);

		double calc_mean_deviation(const coord_list &,
		                           const coord_list &);

		double calc_rmsd(const coord_list &,
		                 const coord_list &);

		std::ostream & operator<<(std::ostream &,
		                          const coord_list &);
	}
}

#endif
