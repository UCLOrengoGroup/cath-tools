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

#ifndef _CATH_TOOLS_SOURCE_CLUSTER_CLUSTER_DOMAINS_HPP
#define _CATH_TOOLS_SOURCE_CLUSTER_CLUSTER_DOMAINS_HPP

#include "cluster/domain_cluster_ids.hpp"

#include <iostream>

namespace cath { namespace common { class id_of_str_bidirnl; } }

namespace cath {
	namespace clust {



		namespace detail {

			/// \brief Store a seq_id and the corresponding domain_cluster_ids
			struct seq_id_and_domain_cluster_ids_pair final {
				/// \brief The seq_id of a sequence
				cluster_id_t seq_id;

				/// \brief The domain_cluster_ids of the domains on the sequence
				domain_cluster_ids dom_cluster_ids;

				/// \brief Standard ctor so this can be used in emplace_back
				seq_id_and_domain_cluster_ids_pair(const cluster_id_t &arg_seq_id,
				                                   domain_cluster_ids  arg_dom_cluster_ids
				                                   ) : seq_id         { arg_seq_id          },
				                                       dom_cluster_ids{ arg_dom_cluster_ids } {
				}


			};

		} // namespace detail




		/// \brief Store a (sparse, sorted) list of seq_id_and_domain_cluster_ids_pair entries, sorted
		class cluster_domains final {
		private:
			/// \brief Lookup from seq_id to the index in the seq_domains vector
			std::unordered_map<cluster_id_t, size_t> index_of_seq_id;

			/// \brief The seq_id_and_domain_cluster_ids_pair entries
			///
			/// This is unsorted; the above unordered_map is used to find the elements
			detail::seq_id_and_domain_cluster_ids_pair_vec seq_domains;

		public:
			/// \brief A const_iterator type alias as part of making this a range over seq_id_and_domain_cluster_ids_pairs
			using const_iterator = detail::seq_id_and_domain_cluster_ids_pair_vec_citr;

			clust_entry_problem add_domain(const cluster_id_t &,
			                               seq::seq_seg_run_opt,
			                               const cluster_id_t &);

			std::vector<cluster_id_t> sorted_seq_ids() const;

			const domain_cluster_ids & domain_cluster_ids_of_seq_id(const cluster_id_t &) const;

			bool empty() const;
			size_t num_seqs() const;

			const_iterator begin() const;
			const_iterator end() const;
		};

		/// \brief Add a domain with the specified (optional) segments in the specified cluster
		///        under the specified sequence ID
		inline clust_entry_problem cluster_domains::add_domain(const cluster_id_t   &arg_seq_id,    ///< The ID of the sequence on which the domain to add appears
		                                                       seq::seq_seg_run_opt  arg_segments,  ///< The (optional) segments of the domain to add
		                                                       const cluster_id_t   &arg_cluster_id ///< The cluster ID of the domain to add
		                                                       ) {
			// Use the unordered map to find the correct seq_id_and_domain_cluster_ids_pair_vec for this seq_id
			// std::cerr << "Looking for seq_id " << arg_seq_id << "\n";
			const auto index_find_itr = index_of_seq_id.find( arg_seq_id );
			auto &the_pair = ( index_find_itr != common::cend( index_of_seq_id ) )
				? seq_domains[ index_find_itr->second ]
				: [&] () -> decltype( auto ) {
					const size_t new_index = seq_domains.size();
					seq_domains.emplace_back( arg_seq_id, domain_cluster_ids{} );
					index_of_seq_id.emplace( arg_seq_id, new_index );
					return seq_domains.back();
				} ();

			// Add the domain
			auto &the_domain_cluster_ids = the_pair.dom_cluster_ids;

			for (const auto &the_domain_cluster_id : the_domain_cluster_ids) {
				const auto intrcn = interaction( arg_segments, the_domain_cluster_id.segments );
				if ( intrcn != clust_entry_problem::NONE ) {
					return intrcn;
				}
			}
			the_domain_cluster_ids.emplace_back( arg_segments, arg_cluster_id );
			return clust_entry_problem::NONE;
		}

		/// \brief Get the sorted list of seq_ids
		inline std::vector<cluster_id_t> cluster_domains::sorted_seq_ids() const {
			std::vector<cluster_id_t> seq_ids;
			for (const auto &x : index_of_seq_id) {
				seq_ids.push_back( x.first );
			}
			std::sort( seq_ids.begin(), seq_ids.end() );
			return seq_ids;
		}

		/// \brief Get the domain_cluster_ids associated with the specified seq ID
		inline const domain_cluster_ids & cluster_domains::domain_cluster_ids_of_seq_id(const cluster_id_t &arg_seq_id ///< The ID of the sequence to query
		                                                                                ) const {
			return seq_domains[ index_of_seq_id.find( arg_seq_id )->second ].dom_cluster_ids;
		}

		/// \brief Return whether this is empty
		inline bool cluster_domains::empty() const {
			return seq_domains.empty();
		}

		/// \brief Return the number of seq_id_and_domain_cluster_ids_pairs
		///
		/// This is called num_seqs() rather than size() to avoid it mistakenly being used for num_entries())
		inline size_t cluster_domains::num_seqs() const {
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

		size_t num_entries(const cluster_domains &);

		std::string to_string(const cluster_domains &,
		                      const common::id_of_str_bidirnl &);

	} // namespace clust
} // namespace cath

#endif
