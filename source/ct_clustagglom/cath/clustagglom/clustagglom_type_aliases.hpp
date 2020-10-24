/// \file
/// \brief The cluster type_aliases header

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

#ifndef _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_CLUSTAGGLOM_TYPE_ALIASES_HPP
#define _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_CLUSTAGGLOM_TYPE_ALIASES_HPP

#include <boost/optional/optional_fwd.hpp>

#include "cath/common/boost_addenda/range/range_concept_type_aliases.hpp"

#include <vector>

namespace cath { namespace clust { class hierarchy_group; } }
namespace cath { namespace clust { class hierarchy_layer; } }
namespace cath { namespace clust { class hierarchy_value; } }
namespace cath { namespace clust { class link_list; } }
namespace cath { namespace clust { struct link; } }
namespace cath { namespace clust { struct merge; } }

namespace cath {
	namespace clust {

		/// \brief Type alias for a vector of hierarchy_group objects
		using hierarchy_group_vec        = std::vector<hierarchy_group>;

		/// \brief Type alias for a vector of hierarchy_layer objects
		using hierarchy_layer_vec        = std::vector<hierarchy_layer>;

		/// \brief Type alias for        a vector of hierarchy_value objects
		using hierarchy_value_vec        = std::vector<hierarchy_value>;

		/// \brief Type alias for the type to using for items in clustering code
		using item_idx                   = uint32_t;

		/// \brief Type alias for a vector of item_idx values
		using item_vec                   = std::vector<item_idx>;

		/// \brief Type alias for an optional item_idx
		using item_opt                   = boost::optional<item_idx>;

		/// \brief Type alias for a vector of item_opt values
		using item_opt_vec               = std::vector<item_opt>;

		/// \brief Type alias for similarities/dissimilarities in clustering code
		using strength                   = float;

		/// \brief Type alias for an optional strength
		using strength_opt               = boost::optional<strength>;

		/// \brief Type alias for a vector of strengths
		using strength_vec               = std::vector<strength>;

		/// \brief Type alias for a tuple of item_idx, item_idx and strength
		using item_item_strength_tpl     = std::tuple<item_idx, item_idx, strength>;

		/// \brief Type alias for a vector of item_item_strength_tpl values
		using item_item_strength_tpl_vec = std::vector<item_item_strength_tpl>;

		/// \brief Type alias for a vector of link values
		using link_vec                   = std::vector<link>;

		/// \brief Type alias for a vector of link_list values
		using link_list_vec              = std::vector<link_list>;

		/// \brief Type alias for link_list_vec's const_iterator type
		using link_list_vec_citr         = common::range_const_iterator_t<link_list_vec>;

		/// \brief Type alias for a vector of merge values
		using merge_vec                  = std::vector<merge>;

		/// \brief Type alias for merge_vec's const_iterator type
		using merge_vec_citr             = common::range_const_iterator_t<merge_vec>;

		/// \brief Type alias for a vector of merge_vec_citr values
		using merge_vec_citr_vec         = std::vector<merge_vec_citr>;
	} // namespace clust
} // namespace cath


#endif // _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_CLUSTAGGLOM_TYPE_ALIASES_HPP
