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

#ifndef _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_HIERARCHY_HIERARCHY_LAYER_HPP
#define _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_HIERARCHY_HIERARCHY_LAYER_HPP

#include "cath/clustagglom/hierarchy/hierarchy_group.hpp"
#include "cath/common/type_aliases.hpp"

namespace cath {
	namespace clust {

		/// \brief One layer of groups of values within a hierarchy
		class hierarchy_layer final {
		private:
			/// \brief The groups
			hierarchy_group_vec groups;

		public:
			/// \brief An iterator type alias as part of making this a range over hierarchy_groups
			using iterator = hierarchy_group_vec::iterator;

			/// \brief A const_iterator type alias as part of making this a range over hierarchy_groups
			using const_iterator = hierarchy_group_vec::const_iterator;

			hierarchy_layer() = default;
			explicit hierarchy_layer(const size_t &);
			explicit hierarchy_layer(hierarchy_group_vec);

			bool empty() const;
			size_t size() const;

			hierarchy_group & operator[](const size_t &);
			const hierarchy_group & operator[](const size_t &) const;

			template <typename... Ts>
			hierarchy_layer & emplace_back(Ts &&...);

			iterator begin();
			iterator end();
			const_iterator begin() const;
			const_iterator end() const;
		};

		/// \brief Construct from a count of hierarchy_groups to default-insert
		inline hierarchy_layer::hierarchy_layer(const size_t &prm_count ///< The number of hierarchy_groups to default insert
		                                        ) : groups( prm_count ) {
		}

		/// \brief Construct from a vector of hierarchy_group values
		inline hierarchy_layer::hierarchy_layer(hierarchy_group_vec prm_hierarchy_group_vec ///< The vector of hierarchy_group values from which this should be constructed
		                                        ) : groups{ std::move( prm_hierarchy_group_vec ) } {
		}

		/// \brief Return whether this is empty
		inline bool hierarchy_layer::empty() const {
			return groups.empty();
		}

		/// \brief Return the number of hierarchy_groups
		inline size_t hierarchy_layer::size() const {
			return groups.size();
		}

		/// \brief Get the hierarchy_group associated with the sequence with the specified index
		inline hierarchy_group & hierarchy_layer::operator[](const size_t &prm_index ///< The index of the hierarchy_group to retrieve
		                                                     ) {
			return groups[ prm_index ];
		}

		/// \brief Get the hierarchy_group associated with the sequence with the specified index
		inline const hierarchy_group & hierarchy_layer::operator[](const size_t &prm_index ///< The index of the hierarchy_group to retrieve
		                                                           ) const {
			return groups[ prm_index ];
		}

		/// \brief Emplace a hierarchy_group to be constructed with the specified arguments
		template <typename... Ts>
		inline hierarchy_layer & hierarchy_layer::emplace_back(Ts &&...args ///< The arguments to perfect-forward to hierarchy_group's ctor
		                                                       ) {
			groups.emplace_back( std::forward<Ts>( args )... );
			return *this;
		}

		/// \brief Standard non-const begin() method, as part of making this into a range over the hierarchy_group entries
		///
		/// \TODO See if it makes more sense to drop the non-const begin()/end() and support whatever
		///       currently uses that some other way
		inline auto hierarchy_layer::begin() -> iterator {
			return std::begin( groups );
		}

		/// \brief Standard non-const end() method, as part of making this into a range over the hierarchy_group entries
		///
		/// \TODO See if it makes more sense to drop the non-const begin()/end() and support whatever
		///       currently uses that some other way
		inline auto hierarchy_layer::end() -> iterator {
			return std::end( groups );
		}

		/// \brief Standard const begin() method, as part of making this into a range over the hierarchy_group entries
		inline auto hierarchy_layer::begin() const -> const_iterator {
			return ::std::cbegin( groups );
		}

		/// \brief Standard const end() method, as part of making this into a range over the hierarchy_group entries
		inline auto hierarchy_layer::end() const -> const_iterator {
			return ::std::cend( groups );
		}

		str_vec to_strings(const hierarchy_layer &);
		std::string to_string(const hierarchy_layer &,
		                      const std::string & = ", ");
		std::ostream & operator<<(std::ostream &,
		                          const hierarchy_layer &);

	} // namespace clust
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_HIERARCHY_HIERARCHY_LAYER_HPP
