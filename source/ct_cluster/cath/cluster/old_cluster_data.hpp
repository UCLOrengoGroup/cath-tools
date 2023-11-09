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

#ifndef CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_OLD_CLUSTER_DATA_HPP
#define CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_OLD_CLUSTER_DATA_HPP

#include <functional>
#include <optional>

#include "cath/cluster/cluster_list.hpp"
#include "cath/cluster/clusters_info.hpp"
#include "cath/common/boost_addenda/range/accumulate_proj.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/container/id_of_str_bidirnl.hpp"

namespace cath::clust {

	/// \brief Store the data associated with old "to" clusters
	///
	/// This uses a reference to an external id_of_str_bidirnl for mapping
	/// from sequences names to IDs. This allows the same IDs to be shared
	/// with other data structures.
	///
	/// \todo Currently unable to check for clashing re-use of sequences?
	class old_cluster_data final {
	private:
		/// \brief A reference_wrapper to an id_of_str_bidirnl object
		std::reference_wrapper<common::id_of_str_bidirnl> id_of_seq_name;

		/// \brief Basic information on the cluster
		///
		/// \todo Is it worthwhile to have this rather than just id_of_str_bidirnl?
		///
		/// Probably not because the clusters store the sizes. Consider switching.
		clusters_info clust_info;

		/// \brief The detailed info on the contents of the clusters
		cluster_list clusters;

	public:
		/// \brief A const_iterator type alias as part of making this a range over cluster_domains entries
		using const_iterator = cluster_list::const_iterator;

		/// \brief An iterator type alias that just duplicates const_iterator to appease some Boost code (Range?)
		using iterator = const_iterator;

		explicit old_cluster_data(common::id_of_str_bidirnl &) noexcept;

		/// Prevent construction from an id_of_str_bidirnl rvalue
		old_cluster_data(const common::id_of_str_bidirnl &&) = delete;

		clust_entry_problem add_entry(const ::std::string_view &,
		                              const ::std::string_view &,
		                              const ::std::string_view &,
		                              seq::seq_seg_run_opt);

		const clusters_info & get_clust_info() const;
		const common::id_of_str_bidirnl & get_id_of_seq_name() const;

		bool empty() const;
		size_t size() const;
		const cluster_domains & operator[](const size_t &x) const;

		const_iterator begin() const;
		const_iterator end() const;
	};

	/// \brief Ctor from an id_of_str_bidirnl reference
	inline old_cluster_data::old_cluster_data(common::id_of_str_bidirnl &prm_id_of_str ///< The id_of_str_bidirnl to use to map from sequences names to IDs
	                                          ) noexcept : id_of_seq_name{ prm_id_of_str } {
	}

	/// \brief Add an entry with the specified sequence name and (optional) segments to the cluster with
	///        the specified name
	inline clust_entry_problem old_cluster_data::add_entry(const ::std::string_view &prm_clust_name, ///< The name of the cluster of the entry
	                                                       const ::std::string_view &prm_seq_id,     ///< The name of the sequence within which this entry appears
	                                                       const ::std::string_view &prm_domain_id,  ///< The name of the entry
	                                                       seq::seq_seg_run_opt      prm_segments    ///< The (optional) segments of the entry within the sequence
	                                                       ) {
		try {
			const auto cluster_id = update_info_and_get_id_for_cluster_of_name(
				clust_info,
				prm_clust_name,
				prm_domain_id,
				prm_segments
			);
			const size_t &seq_id = id_of_seq_name.get().add_name( ::std::string( prm_seq_id ) );
			return clusters.add_domain_to_cluster(
				cluster_id,
				seq_id,
				std::move( prm_segments )
			);
		}
		catch (...) {
			return clust_entry_problem::CLASH;
		}
	}

	/// \brief Get the cluster info
	inline const clusters_info & old_cluster_data::get_clust_info() const {
		return clust_info;
	}

	/// \brief Get the map from sequence name to sequence ID from the specified old_cluster_data
	inline const common::id_of_str_bidirnl & old_cluster_data::get_id_of_seq_name() const {
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
	inline const cluster_domains & old_cluster_data::operator[](const size_t &prm_index ///< The index of the cluster_domains to get
	                                                            ) const {
		return clusters[ prm_index ];
	}

	/// \brief Standard const begin() method, as part of making this into a range over the cluster_domains entries
	inline auto old_cluster_data::begin() const -> const_iterator {
		return ::std::cbegin( clusters );
	}

	/// \brief Standard const end() method, as part of making this into a range over the cluster_domains entries
	inline auto old_cluster_data::end() const -> const_iterator {
		return ::std::cend( clusters );
	}

	/// \brief Get the number of clusters from the specified old_cluster_data
	///
	/// \relates old_cluster_data
	inline size_t get_num_clusters(const old_cluster_data &prm_old_cluster_data ///< The old_cluster_data to query
	                               ) {
		return prm_old_cluster_data.get_clust_info().get_num_clusters();
	}

	/// \brief Get the cluster_info of the cluster with the specified ID from the specified old_cluster_data
	///
	/// \relates old_cluster_data
	inline const cluster_info & get_info_of_cluster_of_id(const old_cluster_data &prm_old_cluster_data, ///< The old_cluster_data to query
	                                                      const cluster_id_t     &prm_cluster_id        ///< The ID of the cluster of interest
	                                                      ) {
		return prm_old_cluster_data.get_clust_info().get_info_of_cluster_of_id( prm_cluster_id );
	}

	/// \brief Get the size of the cluster with the specified ID from the specified old_cluster_data
	///
	/// \relates old_cluster_data
	inline size_t get_size_of_cluster_of_id(const old_cluster_data &prm_old_cluster_data, ///< The old_cluster_data to query
	                                        const cluster_id_t     &prm_cluster_id        ///< The ID of the cluster of interest
	                                        ) {
		return get_info_of_cluster_of_id( prm_old_cluster_data, prm_cluster_id ).get_size();
	}

	/// \brief Get the total number of entries in the clusters of the specified old_cluster_data
	///
	/// \relates old_cluster_data
	inline size_t get_num_entries(const old_cluster_data &prm_old_cluster_data ///< The old_cluster_data to query
	                              ) {
		return common::accumulate_proj(
			common::indices( get_num_clusters( prm_old_cluster_data ) ),
			0_z,
			std::plus<>{},
			[&] (const size_t &x) { return get_size_of_cluster_of_id( prm_old_cluster_data, x ); }
		);
	}

	/// \brief Get the name of the cluster with the specified ID in the specified old_cluster_data
	///
	/// \relates old_cluster_data
	inline const std::string & get_name_of_cluster_of_id(const old_cluster_data &prm_old_cluster_data, ///< The old_cluster_data to query
	                                                     const cluster_id_t     &prm_cluster_id        ///< The ID of the cluster of interest
	                                                     ) {
		return prm_old_cluster_data.get_clust_info().get_name_of_cluster_of_id( prm_cluster_id );
	}

	std::string to_string(const old_cluster_data &);
	::std::optional<ptrdiff_t> largest_number_if_names_all_numeric_integers(const old_cluster_data &);
	::std::optional<ptrdiff_t> largest_number_if_names_all_numeric_integers_of_val_if_none(const old_cluster_data_opt &,
	                                                                                       const ptrdiff_t &);

} // namespace cath::clust

#endif // CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_OLD_CLUSTER_DATA_HPP
