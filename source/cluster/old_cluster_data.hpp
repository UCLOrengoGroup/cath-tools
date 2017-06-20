/// \file
/// \brief The old_cluster_data class header

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

#ifndef _CATH_TOOLS_SOURCE_CLUSTER_OLD_CLUSTER_DATA_H
#define _CATH_TOOLS_SOURCE_CLUSTER_OLD_CLUSTER_DATA_H

#include <boost/range/algorithm/lower_bound.hpp>

#include "cluster/domain_cluster_ids_by_seq.hpp"
#include "cluster/cluster_info.hpp"

namespace cath {
	namespace clust {




		namespace detail {

			/// \brief Store a seq_id and the corresponding domain_cluster_ids
			struct seq_id_and_domain_cluster_ids_pair {
				/// \brief The seq_id of a sequence
				common::id_of_string::id_type seq_id;

				/// \brief The domain_cluster_ids of the domains on the sequence
				domain_cluster_ids dom_cluster_ids;
			};

			/// \brief Type alias for a vector of seq_id_and_domain_cluster_ids_pair objects
			using seq_id_and_domain_cluster_ids_pair_vec = std::vector<seq_id_and_domain_cluster_ids_pair>;

			/// \brief Type alias for seq_id_and_domain_cluster_ids_pair_vec's const_iterator type
			using seq_id_and_domain_cluster_ids_pair_vec_citr = seq_id_and_domain_cluster_ids_pair_vec::const_iterator;

		}




		/// \brief Store a (sparse, sorted) list of seq_id_and_domain_cluster_ids_pair entries, sorted
		class cluster_domains {
		private:
			/// \brief The sparse, sorted list of seq_id_and_domain_cluster_ids_pair entries
			///
			/// \todo Does this really benefit from using domain_cluster_ids
			///       over a vector of boost::optional<seq::seq_seg_run>?
			///
			/// \todo Is it OK that this is sparse and maintains order with lower_bound->insertion?
			detail::seq_id_and_domain_cluster_ids_pair_vec seq_domains;

		public:
			/// \brief A const_iterator type alias as part of making this a range over seq_id_and_domain_cluster_ids_pairs
			using const_iterator = detail::seq_id_and_domain_cluster_ids_pair_vec_citr;

			cluster_domains & add_domain(const common::id_of_string::id_type &,
			                             boost::optional<seq::seq_seg_run>,
			                             const cluster_id_t &);

			const_iterator begin() const;

			const_iterator end() const;
		};

		/// \brief Add a domain with the specified (optional) segments in the specified cluster
		///        under the specified sequence ID
		inline cluster_domains & cluster_domains::add_domain(const common::id_of_string::id_type &arg_seq_id,    ///< The ID of the sequence on which the domain to add appears
		                                                     boost::optional<seq::seq_seg_run>    arg_segments,  ///< The (optional) segments of the domain to add
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

		/// \brief Type alias for a vector of cluster_domains entries
		using cluster_domains_vec = std::vector<cluster_domains>;
		
		/// \brief Type alias for cluster_domains_vec's const_iterator type
		using cluster_domains_vec_citr = cluster_domains_vec::const_iterator;




		/// \brief Store a list of cluster_domains entries, one for each cluster
		class cluster_list {
		private:
			/// \brief The vector of cluster_domains entries, one for each cluster
			cluster_domains_vec cluster_seq_domains;

			/// \brief Ensure that there is a cluster_domains for the specified index and return it
			cluster_domains & ensure_and_get_cluster_domains_of_cluster_id(const size_t &arg_index ///< The index of the cluster_domains to ensure and return
			                                                               ) {
				if ( arg_index >= cluster_seq_domains.size() ) {
					cluster_seq_domains.resize( arg_index + 1 );
				}
				return cluster_seq_domains[ arg_index ];
			}

		public:
			/// \brief A const_iterator type alias as part of making this a range over cluster_domains entries
			using const_iterator = cluster_domains_vec_citr;

			/// \brief Add a domain with the specified (optional) segments in the specified cluster
			///        under the specified sequence ID
			cluster_list & add_domain_to_cluster(const cluster_id_t                  &arg_cluster_id, ///< The cluster ID of the domain to add
			                                     const common::id_of_string::id_type &arg_seq_id,     ///< The ID of the sequence on which the domain to add appears
			                                     boost::optional<seq::seq_seg_run>    arg_segments    ///< The (optional) segments of the domain to add
			                                     ) {
				ensure_and_get_cluster_domains_of_cluster_id( arg_cluster_id ).add_domain(
					arg_seq_id,
					arg_segments,
					arg_cluster_id
				);
				return *this;
			}

			/// \brief Return whether this is empty (ie stores info for no clusters)
			bool empty() const {
				return cluster_seq_domains.empty();
			}

			/// \brief Return the number of cluster_domains entries (ie the number of clusters)
			size_t size() const {
				return cluster_seq_domains.size();
			}

			/// \brief Get the cluster_domains associated with the sequence with the specified index
			const cluster_domains & operator[](const size_t &arg_index ///< The index of the cluster_domains to retrieve
			                                   ) const {
				return cluster_seq_domains[ arg_index ];
			}

			/// \brief Standard const begin() method, as part of making this into a range over the cluster_domains entries
			const_iterator begin() const {
				return common::cbegin( cluster_seq_domains );
			}

			/// \brief Standard const end() method, as part of making this into a range over the cluster_domains entries
			const_iterator end() const {
				return common::cend( cluster_seq_domains );
			}
		};

		/// \brief Type alias for cluster_list's const_iterator type
		using cluster_list_ctir = cluster_list::const_iterator;



		/// \brief Store the data associated with old "to" clusters
		///
		/// This uses a reference to an external id_of_string for mapping
		/// from sequences names to IDs. This allows the same IDs to be shared
		/// with other data structures.
		///
		/// \todo Currently unable to check for clashing re-use of sequences?
		class old_cluster_data final {
		private:
			/// \brief A reference_wrapper to an id_of_string object
			std::reference_wrapper<common::id_of_string> id_of_seq_name;

			/// \brief Basic information on the cluster
			///
			/// \todo Is it worthwhile to have this rather than just cluster_name_ider?
			///
			/// Probably not because the clusters store the sizes. Consider switching.
			cluster_info clust_info;

			/// \brief The detailed info on the contents of the clusters
			cluster_list clusters;

		public:
			/// \brief A const_iterator type alias as part of making this a range over cluster_domains entries
			using const_iterator = cluster_list_ctir;

			/// \brief An iterator type alias that just duplicates const_iterator to appease some Boost code (Range?)
			using iterator = const_iterator;

			explicit old_cluster_data(common::id_of_string &) noexcept;

			/// Prevent construction from an id_of_string rvalue
			old_cluster_data(const common::id_of_string &&) = delete;

			old_cluster_data & add_entry(const boost::string_ref &,
			                             const boost::string_ref &,
			                             boost::optional<seq::seq_seg_run>);

			const cluster_info & get_clust_info() const;
			const common::id_of_string & get_id_of_seq_name() const;

			bool empty() const;
			size_t size() const;
			const cluster_domains & operator[](const size_t &x) const;

			const_iterator begin() const;
			const_iterator end() const;
		};

		/// \brief Ctor from an id_of_string reference
		inline old_cluster_data::old_cluster_data(common::id_of_string &arg_id_of_str ///< The id_of_string to use to map from sequences names to IDs
		                                          ) noexcept : id_of_seq_name{ arg_id_of_str } {
		}

		/// \brief Add an entry with the specified sequence name and (optional) segments to the cluster with
		///        the specified name
		inline old_cluster_data & old_cluster_data::add_entry(const boost::string_ref           &arg_clust_name, ///< The name of the cluster of the entry
		                                                      const boost::string_ref           &arg_seq_id,     ///< The name of the sequence within which this entry appears
		                                                      boost::optional<seq::seq_seg_run>  arg_segments    ///< The (optional) segments of the entry within the sequence
		                                                      ) {
			const auto cluster_id = increment_and_get_id_for_cluster_of_name(
				clust_info,
				arg_clust_name
			);
			const common::id_of_string::id_type &seq_id = id_of_seq_name.get().emplace( arg_seq_id.to_string() ).second;
			clusters.add_domain_to_cluster(
				cluster_id,
				seq_id,
				std::move( arg_segments )
			);
			return *this;
		}

		/// \brief Get the cluster info
		inline const cluster_info & old_cluster_data::get_clust_info() const {
			return clust_info;
		}

		/// \brief Get the map from sequence name to sequence ID from the specified old_cluster_data
		inline const common::id_of_string & old_cluster_data::get_id_of_seq_name() const {
			return id_of_seq_name.get();
		}

		/// \brief Return whether this is empty (ie stores info for no clusters)
		inline bool old_cluster_data::empty() const {
			return clusters.empty();
		}

		/// \brief Return the number of cluster_domains entries (ie the number of clusters)
		inline size_t old_cluster_data::size() const {
			return clusters.size();
		}

		/// \brief Get the cluster_domains associated with the sequence with the specified index
		inline const cluster_domains & old_cluster_data::operator[](const size_t &arg_index ///< The index of the cluster_domains to get
		                                                            ) const {
			return clusters[ arg_index ];
		}

		/// \brief Standard const begin() method, as part of making this into a range over the cluster_domains entries
		inline auto old_cluster_data::begin() const -> const_iterator {
			return common::cbegin( clusters );
		}

		/// \brief Standard const end() method, as part of making this into a range over the cluster_domains entries
		inline auto old_cluster_data::end() const -> const_iterator {
			return common::cend( clusters );
		}

		/// \brief Get the number of clusters from the specified old_cluster_data
		///
		/// \relates old_cluster_data
		inline size_t get_num_clusters(const old_cluster_data &arg_old_cluster_data ///< The old_cluster_data to query
		                               ) {
			return arg_old_cluster_data.get_clust_info().get_num_clusters();
		}

		/// \brief Get the size of the cluster with the specified ID from the specified old_cluster_data
		///
		/// \relates old_cluster_data
		inline size_t get_size_of_cluster_of_id(const old_cluster_data &arg_old_cluster_data, ///< The old_cluster_data to query
		                                        const cluster_id_t     &arg_cluster_id        ///< The ID of the cluster of interest
		                                        ) {
			return arg_old_cluster_data.get_clust_info().get_size_of_cluster_of_id( arg_cluster_id );
		}

		/// \brief Get the name of the cluster with the specified ID in the specified old_cluster_data
		///
		/// \relates old_cluster_data
		inline const std::string & get_name_of_cluster_of_id(const old_cluster_data &arg_old_cluster_data, ///< The old_cluster_data to query
		                                                     const cluster_id_t     &arg_cluster_id        ///< The ID of the cluster of interest
		                                                     ) {
			return arg_old_cluster_data.get_clust_info().get_name_of_cluster_of_id( arg_cluster_id );
		}

		std::string to_string(const old_cluster_data &);

	} // namespace clust
} // namespace cath

#endif
