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

#ifndef _CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_DOMAIN_CLUSTER_IDS_BY_SEQ_HPP
#define _CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_DOMAIN_CLUSTER_IDS_BY_SEQ_HPP

#include <boost/range/adaptor/transformed.hpp>
#include <boost/utility/string_ref.hpp>

#include "cath/cluster/domain_cluster_ids.hpp"
#include "cath/common/container/id_of_str_bidirnl.hpp"

namespace cath {
	namespace clust {

		/// \brief Store domain_cluster_ids under the name of the entry (sequence) to which it belongs
		///
		/// This tallies domain_cluster_ids to entry (sequence) names via a reference_wrapper to an id_of_str_bidirnl
		class domain_cluster_ids_by_seq final {
		private:
			/// \brief A vector of domain_cluster_ids objects
			///
			/// The numbering tallies with that of the id_of_str_bidirnl
			std::vector<domain_cluster_ids> domain_cluster_ids_of_seq_id;

			/// \brief A reference_wrapper to an id_of_str_bidirnl object
			std::reference_wrapper<common::id_of_str_bidirnl> id_of_seq_name;

		public:
			explicit domain_cluster_ids_by_seq(common::id_of_str_bidirnl &);

			/// \brief Disallow construction from a temporary id_of_str_bidirnl
			///        because this stores a reference
			domain_cluster_ids_by_seq(const common::id_of_str_bidirnl &&) = delete;

			clust_entry_problem add(const boost::string_ref &,
			                        domain_cluster_id);

			bool empty() const;
			size_t size() const;

			const common::id_of_str_bidirnl & get_id_of_seq_name() const;

			const domain_cluster_ids & operator[](const size_t &) const;
		};

		/// \brief Ctor from an id_of_str_bidirnl
		inline domain_cluster_ids_by_seq::domain_cluster_ids_by_seq(common::id_of_str_bidirnl &prm_id_of_seq_name ///< The id_of_str_bidirnl to use to map from entry (sequence) names to numeric IDs
		                                                            ) : id_of_seq_name{ prm_id_of_seq_name } {
		}

		/// \brief Add the specified domain_cluster_id under the sequence with the specified name
		///
		/// \TODO try using a id_of_str_bidirnl_view instead of the id_of_str_bidirnl because
		///       that might not be (much) more expensive when storing a new string
		///       and can avoid making a new string when it doesn't need to store.
		///       That said, there's probably less reuse of IDs than for clusters
		///       so it may make little difference.
		inline clust_entry_problem domain_cluster_ids_by_seq::add(const boost::string_ref &prm_seq_id,           ///< The name of the sequence under which to store the domain_cluster_id
		                                                          domain_cluster_id        prm_domain_cluster_id ///< The domain_cluster_id to store
		                                                          ) {
			const auto &id = id_of_seq_name.get().add_name( prm_seq_id.to_string() );
			if ( id >= domain_cluster_ids_of_seq_id.size() ) {
				domain_cluster_ids_of_seq_id.resize( id + 1 );
			}
			domain_cluster_ids &the_domain_cluster_ids = domain_cluster_ids_of_seq_id[ id ];

			for (const auto &the_domain_cluster_id : the_domain_cluster_ids) {
				const auto intrcn = interaction( prm_domain_cluster_id, the_domain_cluster_id );
				if ( intrcn != clust_entry_problem::NONE ) {
					return intrcn;
				}
			}

			the_domain_cluster_ids.emplace_back( std::move( prm_domain_cluster_id ) );
			return clust_entry_problem::NONE;
		}

		/// \brief Return whether this is empty
		inline bool domain_cluster_ids_by_seq::empty() const {
			return domain_cluster_ids_of_seq_id.empty();
		}

		/// \brief Return the number of domain_cluster_ids entries (including possibly empty ones)
		inline size_t domain_cluster_ids_by_seq::size() const {
			return domain_cluster_ids_of_seq_id.size();
		}

		/// \brief Get the id_of_seq_name lookup from sequence name to sequence ID
		inline const common::id_of_str_bidirnl & domain_cluster_ids_by_seq::get_id_of_seq_name() const {
			return id_of_seq_name.get();
		}

		/// \brief Get the domain_cluster_ids associated with the sequence with the specified ID
		inline const domain_cluster_ids & domain_cluster_ids_by_seq::operator[](const size_t &prm_seq_id ///< The ID of the sequence to query
		                                                                        ) const {
			return domain_cluster_ids_of_seq_id[ prm_seq_id ];
		}

		/// \brief Get the cluster ID of the sequence with the specified name from the specified domain_cluster_ids_by_seq
		///
		/// \relates domain_cluster_ids_by_seq
		inline cluster_id_t get_cluster_id_of_seq_name(const domain_cluster_ids_by_seq &prm_dom_cluster_ids_by_seq, ///< The domain_cluster_ids_by_seq to query
		                                               const std::string               &prm_seq_name                ///< The name of the sequence of interest
		                                               ) {
			return prm_dom_cluster_ids_by_seq.get_id_of_seq_name().get_id_of_name( prm_seq_name );
		}

		/// \brief Get whether there is a non-empty domain_cluster_ids for the sequence with the specified id in the specified domain_cluster_ids_by_seq
		///
		/// \relates domain_cluster_ids_by_seq
		inline bool has_domain_cluster_ids_of_seq_id(const domain_cluster_ids_by_seq &prm_dom_cluster_ids_by_seq, ///< The domain_cluster_ids_by_seq to query
		                                             const cluster_id_t              &prm_seq_id                  ///< The id of the sequence of interest
		                                             ) {
			return (
				prm_seq_id < prm_dom_cluster_ids_by_seq.size()
				&&
				! prm_dom_cluster_ids_by_seq[ prm_seq_id ].empty()
			);
		}

		/// \brief Get whether there is a non-empty domain_cluster_ids for the sequence with the specified name in the specified domain_cluster_ids_by_seq
		///
		/// \relates domain_cluster_ids_by_seq
		inline bool has_domain_cluster_ids_of_seq_name(const domain_cluster_ids_by_seq &prm_dom_cluster_ids_by_seq, ///< The domain_cluster_ids_by_seq to query
		                                               const std::string               &prm_seq_name                ///< The name of the sequence of interest
		                                               ) {
			return has_domain_cluster_ids_of_seq_id(
				prm_dom_cluster_ids_by_seq,
				get_cluster_id_of_seq_name(
					prm_dom_cluster_ids_by_seq,
					prm_seq_name
				)
			);
		}

		/// \brief Get the domain_cluster_ids of the sequence with the specified name from the specified domain_cluster_ids_by_seq
		///
		/// \relates domain_cluster_ids_by_seq
		inline const domain_cluster_ids & get_domain_cluster_ids_of_seq_name(const domain_cluster_ids_by_seq &prm_dom_cluster_ids_by_seq, ///< The domain_cluster_ids_by_seq to query
		                                                                     const std::string               &prm_seq_name                ///< The name of the sequence of interest
		                                                                     ) {
			return prm_dom_cluster_ids_by_seq[ get_cluster_id_of_seq_name( prm_dom_cluster_ids_by_seq, prm_seq_name ) ];
		}

	} // namespace clust
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_DOMAIN_CLUSTER_IDS_BY_SEQ_HPP
