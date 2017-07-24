/// \file
/// \brief The cluster_domains class header

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

#ifndef _CATH_TOOLS_SOURCE_CLUSTER_CLUSTER_DOMAINS_H
#define _CATH_TOOLS_SOURCE_CLUSTER_CLUSTER_DOMAINS_H

#include "cluster/domain_cluster_ids.hpp"
#include "common/container/id_of_string.hpp"

namespace cath {
	namespace clust {



		namespace detail {

			/// \brief Store a seq_id and the corresponding domain_cluster_ids
			struct seq_id_and_domain_cluster_ids_pair final {
				/// \brief The seq_id of a sequence
				common::id_of_string::id_type seq_id;

				/// \brief The domain_cluster_ids of the domains on the sequence
				domain_cluster_ids dom_cluster_ids;
			};

		}




		/// \brief Store a (sparse, sorted) list of seq_id_and_domain_cluster_ids_pair entries, sorted
		class cluster_domains final {
		private:
			/// \brief The sparse, sorted list of seq_id_and_domain_cluster_ids_pair entries
			///
			/// \todo Does this really benefit from using domain_cluster_ids
			///       over a vector of seq::seq_seg_run_opt?
			///
			/// \todo Is it OK that this is sparse and maintains order with lower_bound->insertion?
			detail::seq_id_and_domain_cluster_ids_pair_vec seq_domains;

		public:
			/// \brief A const_iterator type alias as part of making this a range over seq_id_and_domain_cluster_ids_pairs
			using const_iterator = detail::seq_id_and_domain_cluster_ids_pair_vec_citr;

			cluster_domains & add_domain(const common::id_of_string::id_type &,
			                             seq::seq_seg_run_opt,
			                             const cluster_id_t &);

			bool empty() const;
			size_t size() const;

			const_iterator begin() const;
			const_iterator end() const;
		};

		/// \brief Add a domain with the specified (optional) segments in the specified cluster
		///        under the specified sequence ID
		inline cluster_domains & cluster_domains::add_domain(const common::id_of_string::id_type &arg_seq_id,    ///< The ID of the sequence on which the domain to add appears
		                                                     seq::seq_seg_run_opt                 arg_segments,  ///< The (optional) segments of the domain to add
		                                                     const cluster_id_t                  &arg_cluster_id ///< The cluster ID of the domain to add
		                                                     ) {
			// Find an iterator to the part of seq_domains where this entry this entry's
			// seq_id_and_domain_cluster_ids_pair is or would be
			const auto find_itr = boost::range::lower_bound(
				seq_domains,
				arg_seq_id,
				[] (const auto &x, const auto &y) {
					return x.seq_id < y;
				}
			);
			// Insert a seq_id_and_domain_cluster_ids_pair if necessary and, either way, grab an iterator to
			// the one to which the domain should be added
			const auto emplace_itr = ( find_itr == common::cend( seq_domains ) || find_itr->seq_id != arg_seq_id )
				? seq_domains.insert( find_itr, detail::seq_id_and_domain_cluster_ids_pair{ arg_seq_id, domain_cluster_ids{} } )
				: find_itr;

			// Add the domain and return *this
			emplace_itr->dom_cluster_ids.emplace_back( arg_segments, arg_cluster_id );
			return *this;
		}

		/// \brief Return whether this is empty
		inline bool cluster_domains::empty() const {
			return seq_domains.empty();
		}

		/// \brief Return the number of seq_id_and_domain_cluster_ids_pairs
		inline size_t cluster_domains::size() const {
			return seq_domains.size();
		}

		/// \brief Standard const begin() method, as part of making this into a range over the seq_id_and_domain_cluster_ids_pairs
		inline auto cluster_domains::begin() const -> const_iterator {
			return common::cbegin( seq_domains );
		}

		/// \brief Standard const end() method, as part of making this into a range over the seq_id_and_domain_cluster_ids_pairs
		inline auto cluster_domains::end() const -> const_iterator {
			return common::cend  ( seq_domains );
		}

		std::string to_string(const cluster_domains &,
		                      const common::id_of_string &);

	} // namespace clust
} // namespace cath

#endif
