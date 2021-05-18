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

#ifndef _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_HIERARCHY_HIERARCHY_VALUE_HPP
#define _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_HIERARCHY_HIERARCHY_VALUE_HPP

#include "cath/clustagglom/clustagglom_type_aliases.hpp"
#include "cath/clustagglom/hierarchy/hierarchy_ref.hpp"

#include <iosfwd>

namespace cath {
	namespace clust {

		/// \brief A value in a hierarchy that contains one entry (ie leaf node) or
		///        one cluster from the next deepest layer
		class hierarchy_value final {
		private:
			/// \brief Whether this is an entry (ie leaf node) or a cluster
			///        from the next deepest layer of the hierarchy
			///
			/// This is in std::variant-ish territory and the implementation could switch
			/// that way depending in the future
			hierarchy_ref type;

			/// \brief The index of the entry (leaf node) or the index of the group
			///        in the next deepest layer
			item_idx      index;

		public:
			constexpr hierarchy_value(const hierarchy_ref &,
			                          const item_idx &) noexcept;

			[[nodiscard]] constexpr const hierarchy_ref &get_type() const;
			[[nodiscard]] constexpr const item_idx &     get_index() const;

			hierarchy_value & set_type(const hierarchy_ref &);
			hierarchy_value & set_index(const item_idx &);
		};

		/// \brief Ctor from the type and index
		inline constexpr hierarchy_value::hierarchy_value(const hierarchy_ref &prm_type, ///< Whether this is an entry (ie leaf node) or a cluster from the next deepest layer of the hierarchy
		                                                  const item_idx      &prm_index ///< The index of the entry (leaf node) or the index of the group in the next deepest layer
		                                                  ) noexcept : type { prm_type  },
		                                                               index{ prm_index } {
		}

		/// \brief Getter for whether this is an entry (ie leaf node) or a cluster from the next deepest layer of the hierarchy
		inline constexpr const hierarchy_ref & hierarchy_value::get_type() const {
			return type;
		}

		/// \brief Getter for the index of the entry (leaf node) or the index of the group in the next deepest layer
		inline constexpr const item_idx & hierarchy_value::get_index() const {
			return index;
		}

		/// \brief Setter for whether this is an entry (ie leaf node) or a cluster from the next deepest layer of the hierarchy
		inline hierarchy_value & hierarchy_value::set_type(const hierarchy_ref &prm_type ///< Whether this is an entry (ie leaf node) or a cluster from the next deepest layer of the hierarchy
		                                                   ) {
			type = prm_type;
			return *this;
		}

		/// \brief Setter for the index of the entry (leaf node) or the index of the group in the next deepest layer
		inline hierarchy_value & hierarchy_value::set_index(const item_idx &prm_index ///< The index of the entry (leaf node) or the index of the group in the next deepest layer
		                                                    ) {
			index = prm_index;
			return *this;
		}

		std::string to_string(const hierarchy_value &);
		std::ostream & operator<<(std::ostream &,
		                          const hierarchy_value &);

	} // namespace clust
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_HIERARCHY_HIERARCHY_VALUE_HPP
