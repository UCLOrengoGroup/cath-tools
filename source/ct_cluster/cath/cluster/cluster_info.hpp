/// \file
/// \brief The cluster_info class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_CLUSTER_INFO_HPP
#define _CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_CLUSTER_INFO_HPP

#include <boost/operators.hpp>
#include <boost/utility/string_ref.hpp>

#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/optional/make_optional_if.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/seq/seq_seg_run.hpp"

#include <string>

namespace cath {
	namespace clust {

		/// \brief Store the summary information on cluster, which can be used for determining
		///        the mapping threshold and for ordering unmapped, new clusters
		class cluster_info final : private boost::less_than_comparable<cluster_info> {
		private:
			/// \brief The size of the cluster
			size_t size;

			/// \brief The total of the square-roots of the domains' lengths
			doub_opt total_sqrt_length;

			/// \brief The total of the domains' middle indices (via middle_index(const seq_seg_run &))
			doub_opt total_mid_point_index;

			/// \brief The lowest domain ID in the cluster
			std::string lowest_domain_id;

		public:
			/// \brief Default ctor
			cluster_info() = default;

			cluster_info & add_entry(const boost::string_ref &,
			                         const seq::seq_seg_run_opt &);

			[[nodiscard]] const size_t &     get_size() const;
			[[nodiscard]] const doub_opt &   get_total_sqrt_length() const;
			[[nodiscard]] const doub_opt &   get_total_mid_point_index() const;
			[[nodiscard]] const std::string &get_lowest_domain_id() const;
		};

		bool operator<(const cluster_info &,
		               const cluster_info &);

		/// \brief Update the cluster to account for a new entry with the specified name and (optional) segments
		inline cluster_info & cluster_info::add_entry(const boost::string_ref    &prm_name,    ///< The name of the new entry
		                                              const seq::seq_seg_run_opt &prm_segments ///< The (optional) segments of the new entry
		                                              ) {
			if ( size == 0 ) {
				lowest_domain_id = prm_name.to_string();
				if ( prm_segments ) {
					total_sqrt_length     = 0.0;
					total_mid_point_index = 0_z;
				}
			}
			else {
				if ( static_cast<bool>( prm_segments) != static_cast<bool>( total_sqrt_length ) ) {
					BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot mix specifying domain regions within a cluster"));
				}
			}
			++size;
			if ( prm_name < lowest_domain_id ) {
				lowest_domain_id = prm_name.to_string();
			}
			if ( prm_segments ) {
				*total_sqrt_length     += sqrt( static_cast<double>( get_total_length( *prm_segments ) ) );
				*total_mid_point_index += middle_index( *prm_segments );
			}
			return *this;
		}

		/// \brief Get the number of entries
		inline const size_t & cluster_info::get_size() const {
			return size;
		}

		/// \brief Getter for the total sqrt length
		///
		/// \pre `get_size() > 0` else an invalid_argument_exception will be thrown
		inline const doub_opt & cluster_info::get_total_sqrt_length() const {
			if ( get_size() <= 0 ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot get_total_sqrt_length() from cluster_info that's empty"));
			}
			return total_sqrt_length;
		}

		/// \brief Getter for the total mid-point index
		///
		/// \pre `get_size() > 0` else an invalid_argument_exception will be thrown
		inline const doub_opt & cluster_info::get_total_mid_point_index() const {
			if ( get_size() <= 0 ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot get_total_mid_point_index() from cluster_info that's empty"));
			}
			return total_mid_point_index;
		}

		/// \brief Getter for the lowest domain ID
		///
		/// \pre `get_size() > 0` else an invalid_argument_exception will be thrown
		inline const std::string & cluster_info::get_lowest_domain_id() const {
			if ( get_size() <= 0 ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot get_lowest_domain_id() from cluster_info that's empty"));
			}
			return lowest_domain_id;
		}

		/// \brief Get the average mid point index of the specified cluster
		inline doub_opt get_average_mid_point_index(const cluster_info &prm_cluster_info ///< The cluster_info to query
		                                            ) {
			return if_then_optional(
				prm_cluster_info.get_total_sqrt_length()
				&&
				prm_cluster_info.get_total_mid_point_index(),
				(
					static_cast<double>( *prm_cluster_info.get_total_mid_point_index() )
					/
					static_cast<double>(  prm_cluster_info.get_size                 () )
				)
			);
		}

		/// \brief Return whether the first specified cluster_info should be placed before the second when ordering
		///        unmapped new clusters
		///
		///  * descending on sum over domains of sqrt(total_dom_length) (ie earlier FunFams have more/longer sequences, with more emphasis on having more sequence)
		///  * descending on number of sequences
		///  * ascending on average mid-point index
		///  * ascending on first domain ID
		inline bool operator<(const cluster_info &prm_lhs, ///< The first  cluster to compare
		                      const cluster_info &prm_rhs  ///< The second cluster to compare
		                      ) {
			const auto average_mid_point_index_lhs = get_average_mid_point_index( prm_lhs );
			const auto average_mid_point_index_rhs = get_average_mid_point_index( prm_rhs );
			return (
				std::tie(
					prm_rhs.get_total_sqrt_length(),
					prm_rhs.get_size(),
					average_mid_point_index_lhs,
					prm_lhs.get_lowest_domain_id()
				)
				<
				std::tie(
					prm_lhs.get_total_sqrt_length(),
					prm_lhs.get_size(),
					average_mid_point_index_rhs,
					prm_rhs.get_lowest_domain_id()
				)
			);
		}

	} // namespace clust
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_CLUSTER_INFO_HPP
