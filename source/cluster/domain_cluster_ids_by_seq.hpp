/// \file
/// \brief The domain_cluster_ids_by_seq class header

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

#ifndef _CATH_TOOLS_SOURCE_CLUSTER_DOMAIN_CLUSTER_IDS_BY_SEQ_H
#define _CATH_TOOLS_SOURCE_CLUSTER_DOMAIN_CLUSTER_IDS_BY_SEQ_H

#include <boost/range/adaptor/transformed.hpp>
#include <boost/utility/string_ref.hpp>

#include "cluster/domain_cluster_ids.hpp"
#include "common/container/id_of_string.hpp"

// #include <string>

namespace cath {
	namespace clust {

		/// \brief Store domain_cluster_ids under the name of the entry (sequence) to which it belongs
		///
		/// This tallies domain_cluster_ids to entry (sequence) names via a reference_wrapper to an id_of_string
		class domain_cluster_ids_by_seq final {
		private:
			/// \brief A vector of domain_cluster_ids objects
			///
			/// The numbering tallies with that of the id_of_string
			std::vector<domain_cluster_ids> domain_cluster_ids_of_seq_id;

			/// \brief A reference_wrapper to an id_of_string object
			std::reference_wrapper<common::id_of_string> id_of_seq_name;

		public:
			domain_cluster_ids_by_seq(common::id_of_string &);

			/// \brief Disallow construction from a temporary id_of_string
			///        because this stores a reference
			domain_cluster_ids_by_seq(const common::id_of_string &&) = delete;

			domain_cluster_ids_by_seq & add(const boost::string_ref &,
			                                domain_cluster_id);

			const common::id_of_string & get_id_of_seq_name() const;

			const domain_cluster_ids & operator[](const size_t &) const;
		};

		/// \brief Ctor from an id_of_string
		inline domain_cluster_ids_by_seq::domain_cluster_ids_by_seq(common::id_of_string &arg_id_of_seq_name ///< The id_of_string to use to map from entry (sequence) names to numeric IDs
		                                                            ) : id_of_seq_name{ arg_id_of_seq_name } {
		}

		/// \brief Add the specified domain_cluster_id under the sequence with the specified name
		///
		/// \TODO try using a id_of_string_view instead of the id_of_string because
		///       that might not be (much) more expensive when storing a new string
		///       and can avoid making a new string when it doesn't need to store.
		///       That said, there's probably less reuse of IDs than for clusters
		///       so it may make little difference.
		inline domain_cluster_ids_by_seq & domain_cluster_ids_by_seq::add(const boost::string_ref &arg_seq_id,           ///< The name of the sequence under which to store the domain_cluster_id
		                                                                  domain_cluster_id        arg_domain_cluster_id ///< The domain_cluster_id to store
		                                                                  ) {
			const auto &id = id_of_seq_name.get().emplace( arg_seq_id.to_string() ).second;
			if ( id >= domain_cluster_ids_of_seq_id.size() ) {
				domain_cluster_ids_of_seq_id.resize( id + 1 );
			}
			domain_cluster_ids_of_seq_id[ id ].emplace_back( std::move( arg_domain_cluster_id ) );
			return *this;
		}

		/// \brief Get the id_of_seq_name lookup from sequence name to sequence ID
		inline const common::id_of_string & domain_cluster_ids_by_seq::get_id_of_seq_name() const {
			return id_of_seq_name.get();
		}

		/// \brief Get the domain_cluster_ids associated with the sequence with the specified ID
		inline const domain_cluster_ids & domain_cluster_ids_by_seq::operator[](const size_t &arg_seq_id ///< The ID of the sequence to query
		                                                                        ) const {
			return domain_cluster_ids_of_seq_id[ arg_seq_id ];
		}

		/// \brief Get the cluster ID of the sequence with the specified name from the specified domain_cluster_ids_by_seq
		///
		/// \relates domain_cluster_ids_by_seq
		inline common::id_of_string::id_type get_cluster_id_of_seq_name(const domain_cluster_ids_by_seq &arg_dom_cluster_ids_by_seq, ///< The domain_cluster_ids_by_seq to query
		                                                                const std::string               &arg_seq_name                ///< The name of the sequence of interest
		                                                                ) {
			return arg_dom_cluster_ids_by_seq.get_id_of_seq_name()[ arg_seq_name ];
		}

		/// \brief Get the domain_cluster_ids of the sequence with the specified name from the specified domain_cluster_ids_by_seq
		///
		/// \relates domain_cluster_ids_by_seq
		inline const domain_cluster_ids & get_domain_cluster_ids_of_seq_name(const domain_cluster_ids_by_seq &arg_dom_cluster_ids_by_seq, ///< The domain_cluster_ids_by_seq to query
		                                                                     const std::string               &arg_seq_name                ///< The name of the sequence of interest
		                                                                     ) {
			return arg_dom_cluster_ids_by_seq[ get_cluster_id_of_seq_name( arg_dom_cluster_ids_by_seq, arg_seq_name ) ];
		}

	} // namespace clust
} // namespace cath

#endif
