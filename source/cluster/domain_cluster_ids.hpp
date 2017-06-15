/// \file
/// \brief The domain_cluster_ids class header

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

#ifndef _CATH_TOOLS_SOURCE_CLUSTER_DOMAIN_CLUSTER_IDS_H
#define _CATH_TOOLS_SOURCE_CLUSTER_DOMAIN_CLUSTER_IDS_H

#include <boost/optional.hpp>

#include "cluster/cluster_type_aliases.hpp"
#include "seq/seq_seg_run.hpp"

namespace cath {
	namespace clust {

		/// \brief The segments of a domain (or none if the domain covers the full entry (sequence)) and the domain's cluster ID
		struct domain_cluster_id {

			/// \brief The segments of a domain (or none if the domain covers the full entry (sequence))
			boost::optional<seq::seq_seg_run> segments;

			/// \brief The cluster ID of the domain
			cluster_id_t cluster_id;

			/// \brief Ctor from the segments and cluster ID of a domain
			domain_cluster_id(boost::optional<seq::seq_seg_run>  arg_segments,  ///< The segments of a domain (or none if the domain covers the full entry (sequence))
			                  const cluster_id_t                &arg_cluster_id ///< The cluster ID of the domain
			                  ) : segments   { std::move( arg_segments ) },
			                      cluster_id { arg_cluster_id            } {
			}
		};

		/// \brief A list of domain_cluster_ids
		///
		/// \TODO Consider making this enforce:
		///   * No overlapping segments
		///   * Keep in some sort of order?
		class domain_cluster_ids {
		private:
			/// \brief The vector of domain_cluster_id objects
			std::vector<domain_cluster_id> dom_clust_ids;

		public:
			/// \brief A const_iterator type alias as part of making this a range over domain_cluster_ids
			using const_iterator = domain_cluster_id_vec::const_iterator;

			/// \brief An iterator type alias as part of making this a range over domain_cluster_ids
			///
			/// This is just an alias to the const_iterator because this isn't a genuinely mutable range
			using iterator = const_iterator;

			/// \brief Emplace_back a domain_cluster_id
			template <typename... Ts>
			inline void emplace_back(Ts &&...args ///< The arguments to pass to the domain_cluster_id ctor
			                         ) {
				dom_clust_ids.emplace_back( std::forward<Ts>( args )... );
			}

			/// \brief Whether this is empty
			inline bool empty() const {
				return dom_clust_ids.empty();
			}

			/// \brief The number of domain_cluster_ids
			inline size_t size() const {
				return dom_clust_ids.size();
			}

			/// \brief Const-overload of standard subscript operator
			inline const domain_cluster_id & operator[](const size_t &arg_index ///< The index of the domain_cluster_id to access
			                                            ) {
				return dom_clust_ids[ arg_index ];
			}

			/// \brief Standard const begin() method, as part of making this a range over domain_cluster_ids
			inline auto begin() const -> const_iterator {
				return common::cbegin( dom_clust_ids );
			}

			/// \brief Standard const end() method, as part of making this a range over domain_cluster_ids
			inline auto end() const -> const_iterator {
				return common::cend( dom_clust_ids );
			}
		};

		std::string to_string(const domain_cluster_id &);
		std::string to_string(const domain_cluster_ids &);

	} // namespace clust
} // namespace cath

#endif
