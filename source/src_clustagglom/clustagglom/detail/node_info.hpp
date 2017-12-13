/// \file
/// \brief The node_info class header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_CLUSTAGGLOM_CLUSTAGGLOM_DETAIL_NODE_INFO_HPP
#define _CATH_TOOLS_SOURCE_SRC_CLUSTAGGLOM_CLUSTAGGLOM_DETAIL_NODE_INFO_HPP

#include <boost/optional/optional_fwd.hpp>

#include <vector>

namespace cath {
	namespace clust {
		namespace detail {

			/// \brief Represent the information on a cluster (ie a root node within a cutoff level)
			struct node_info {
				/// \brief The location of the cluster (ie the index under which it's being accessed)
				size_t locn;

				/// \brief The layer of the cluster
				size_t layer;

				constexpr node_info(const size_t &,
				                    const size_t &) noexcept;
			};

			/// \brief Standard ctor like aggregate initialisation
			inline constexpr node_info::node_info(const size_t &arg_locn, ///< The location of the cluster (ie the index under which it's being accessed)
			                                      const size_t &arg_layer ///< The layer of the cluster
			                                      ) noexcept : locn  { arg_locn  },
			                                                   layer { arg_layer } {
			}

			/// \brief Type alias for an optional node_info
			using node_info_opt     = boost::optional<node_info>;

			/// \brief Type alias for a vector of node_info_opt values
			using node_info_opt_vec = std::vector<node_info_opt>;

		} // namespace detail
	} // namespace clust
} // namespace cath

#endif
