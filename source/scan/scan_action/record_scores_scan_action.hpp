/// \file
/// \brief The record_scores_scan_action class header

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

#ifndef _CATH_TOOLS_SOURCE_SCAN_SCAN_ACTION_RECORD_SCORES_SCAN_ACTION_H
#define _CATH_TOOLS_SOURCE_SCAN_SCAN_ACTION_RECORD_SCORES_SCAN_ACTION_H

#include "scan/detail/res_pair/single_struc_res_pair.hpp"

#include "scan/scan_query_set.hpp"
namespace cath {
	namespace scan {

		/// \brief TODOCUMENT
		///
		/// \todo Consider changing to a doub_vec_vec
		///       (begause, eg, an all-vs-all within a set of 12,000 would require >1Gb)
		class record_scores_scan_action final {
		private:
			/// \brief TODOCUMENT
			const size_t num_queries;

			/// \brief TODOCUMENT
			const size_t num_matches;

			/// \brief TODOCUMENT
			doub_vec scores;

			double & get_entry(const size_t &,
			                   const size_t &);
			const double & get_entry(const size_t &,
			                         const size_t &) const;

		public:
			record_scores_scan_action(const size_t &,
			                          const size_t &);

			void operator()(const detail::single_struc_res_pair &,
			                const detail::single_struc_res_pair &,
			                const index_type &,
			                const index_type &);

			const double & get_score(const size_t &,
			                         const size_t &) const;
		};

		/// \brief TODOCUMENT
		inline double & record_scores_scan_action::get_entry(const size_t &arg_query_index, ///< TODOCUMENT
		                                                     const size_t &arg_match_index  ///< TODOCUMENT
		                                                     ) {
			return const_cast<double &>(
				static_cast<const record_scores_scan_action &>( *this ).get_entry(
					arg_query_index,
					arg_match_index
				)
			);
		}

		/// \brief TODOCUMENT
		inline const double & record_scores_scan_action::get_entry(const size_t &arg_query_index, ///< TODOCUMENT
		                                                           const size_t &arg_match_index  ///< TODOCUMENT
		                                                           ) const {
			return scores[ arg_query_index * num_matches + arg_match_index ];
		}

		/// \brief TODOCUMENT
		inline record_scores_scan_action::record_scores_scan_action(const size_t &arg_num_queries, ///< TODOCUMENT
		                                                            const size_t &arg_num_matches  ///< TODOCUMENT
		                                                            ) : num_queries ( arg_num_queries ),
		                                                                num_matches ( arg_num_matches ),
		                                                                scores      ( num_queries * num_matches, 0.0 ) {
		}

		/// \brief TODOCUMENT
		inline void record_scores_scan_action::operator()(const detail::single_struc_res_pair &arg_res_pair_a,  ///< TODOCUMENT
		                                                  const detail::single_struc_res_pair &arg_res_pair_b,  ///< TODOCUMENT
		                                                  const index_type                    &arg_structure_a, ///< TODOCUMENT
		                                                  const index_type                    &arg_structure_b  ///< TODOCUMENT
		                                                  ) {
			const auto distance = sqrt( detail::squared_distance( arg_res_pair_a, arg_res_pair_b ) );
			get_entry( arg_structure_a, arg_structure_b ) += 1.0 - ( distance / 7.0 );
		}

		/// \brief TODOCUMENT
		inline const double & record_scores_scan_action::get_score(const size_t &arg_query_index, ///< TODOCUMENT
		                                                           const size_t &arg_match_index  ///< TODOCUMENT
		                                                           ) const {
			return get_entry( arg_query_index, arg_match_index );
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		record_scores_scan_action make_record_scores_scan_action(const scan_query_set<KPs...> &arg_query_set, ///< TODOCUMENT
		                                                         const scan_index<KPs...>     &arg_index      ///< TODOCUMENT
		                                                         ) {
			return {
				arg_query_set.get_num_structures(),
				arg_index.get_num_structures()
			};
		}
	} // namespace scan
} // namespace cath

#endif
