/// \file
/// \brief The hierarchy class header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_CLUSTAGGLOM_CLUSTAGGLOM_HIERARCHY_HIERARCHY_GROUP_HPP
#define _CATH_TOOLS_SOURCE_SRC_CLUSTAGGLOM_CLUSTAGGLOM_HIERARCHY_HIERARCHY_GROUP_HPP

#include "cath/clustagglom/hierarchy/hierarchy_value.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"

namespace cath {
	namespace clust {

		/// \brief A group of hierarchy_values within a hierarchy_layer, within a hierarchy
		class hierarchy_group final {
		private:
			/// \brief The hierarchy_values
			hierarchy_value_vec values;

		public:
			/// \brief An iterator type alias as part of making this a range over hierarchy_values
			using iterator       = hierarchy_value_vec::iterator;

			/// \brief A const_iterator type alias as part of making this a range over hierarchy_values
			using const_iterator = hierarchy_value_vec::const_iterator;

			bool empty() const;
			size_t size() const;

			const hierarchy_value & operator[](const size_t &) const;

			template <typename... Ts>
			hierarchy_group & emplace_back(Ts &&...);

			iterator begin();
			iterator end();
			const_iterator begin() const;
			const_iterator end() const;
		};

		/// \brief Return whether this is empty
		inline bool hierarchy_group::empty() const {
			return values.empty();
		}

		/// \brief Return the number of hierarchy_values
		inline size_t hierarchy_group::size() const {
			return values.size();
		}

		/// \brief Get the hierarchy_value associated with the sequence with the specified index
		inline const hierarchy_value & hierarchy_group::operator[](const size_t &prm_index ///< The index of the hierarchy_value to retrieve
		                                                           ) const {
			return values[ prm_index ];
		}

		/// \brief Emplace a hierarchy_value to be constructed with the specified arguments
		template <typename... Ts>
		inline hierarchy_group & hierarchy_group::emplace_back(Ts &&...args ///< The arguments to perfect-forward to hierarchy_value's ctor
		                                                       ) {
			values.emplace_back( std::forward<Ts>( args )... );
			return *this;
		}

		/// \brief Standard non-const begin() method, as part of making this into a range over the hierarchy_value entries
		///
		/// \TODO See if it makes more sense to drop the non-const begin()/end() and support whatever
		///       currently uses that some other way
		inline auto hierarchy_group::begin() -> iterator {
			return std::begin( values );
		}

		/// \brief Standard non-const end() method, as part of making this into a range over the hierarchy_value entries
		///
		/// \TODO See if it makes more sense to drop the non-const begin()/end() and support whatever
		///       currently uses that some other way
		inline auto hierarchy_group::end() -> iterator {
			return std::end( values );
		}

		/// \brief Standard const begin() method, as part of making this into a range over the hierarchy_value entries
		inline auto hierarchy_group::begin() const -> const_iterator {
			return common::cbegin( values );
		}

		/// \brief Standard const end() method, as part of making this into a range over the hierarchy_value entries
		inline auto hierarchy_group::end() const -> const_iterator {
			return common::cend( values );
		}

		std::string to_string(const hierarchy_group &);
		std::ostream & operator<<(std::ostream &,
		                          const hierarchy_group &);

	} // namespace clust
} // namespace cath

#endif
