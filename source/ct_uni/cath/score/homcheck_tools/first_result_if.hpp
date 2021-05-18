/// \file
/// \brief The first_result_if header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_HOMCHECK_TOOLS_FIRST_RESULT_IF_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_HOMCHECK_TOOLS_FIRST_RESULT_IF_HPP

#include <optional>
#include <vector>

#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/algorithm/sort.hpp>

#include "cath/common/algorithm/copy_build.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/range/adaptor/limited.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/file/file_type_aliases.hpp"
#include "cath/score/homcheck_tools/ssaps_and_prcs_of_query.hpp"

namespace cath { namespace homcheck { class ssaps_and_prcs_of_query; } }

namespace cath {
	namespace homcheck {
		namespace detail {


			/// \brief Return the first result ordered by the specified less-than function that meets the specified predicate
			template <typename RES, typename LT_FN, typename PRED_FN>
			auto first_result_if_impl(const RES &prm_results,   ///< The results to query
			                          LT_FN      prm_less_than, ///< The less-than function to determine which results should come first ("lower" comes earlier)
			                          PRED_FN    prm_pred       ///< The predicate to specify acceptable results
			                          ) -> decltype( ::std::make_optional( std::cref( prm_results[ 0 ] ) ) ) {
				auto indices = common::copy_build<size_vec>( common::indices( prm_results.size() ) );
				boost::range::sort(
					indices,
					[&] (const size_t &x, const size_t &y) {
						return prm_less_than( prm_results[ x ], prm_results[ y ] );
					}
				);
				const auto find_itr = boost::range::find_if(
					indices,
					[&] (const size_t &x) {
						return prm_pred( prm_results[ x ] );
					}
				);
				return ( find_itr != ::std::cend( indices ) ) ? ::std::make_optional( std::cref( prm_results[ *find_itr ] ) )
				                                               : ::std::nullopt;
			}

			/// \brief Return the first N results ordered by the specified less-than function that meet the specified predicate
			template <typename RES, typename LT_FN, typename PRED_FN>
			auto first_n_results_if_impl(const RES    &prm_results,   ///< The results to query
			                             LT_FN         prm_less_than, ///< The less-than function to determine which results should come first ("lower" comes earlier)
			                             PRED_FN       prm_pred,      ///< The predicate to specify acceptable results
			                             const size_t &prm_n          ///< The maximum number of results to return
			                             ) -> std::vector< std::remove_const_t< std::remove_reference_t< decltype( prm_results[ 0 ] ) > > > {
				using result_type        = std::remove_const_t< std::remove_reference_t< decltype( prm_results[ 0 ] ) > >;
				using result_vector_type = std::vector<result_type>;
				auto indices = common::copy_build<size_vec>( common::indices( prm_results.size() ) );
				boost::range::sort(
					indices,
					[&] (const size_t &x, const size_t &y) {
						return prm_less_than( prm_results[ x ], prm_results[ y ] );
					}
				);

				return common::transform_build<result_vector_type> (
					indices
						| boost::adaptors::filtered( [&] (const size_t &x) { return prm_pred( prm_results[ x ] ); } )
						| common::limited( prm_n ),
					[&] (const size_t &x) {
						return prm_results[ x ];
					}
				);
			}
		} // namespace detail

		/// \brief Return the first SSAP and PRC result ordered by the specified less-than function that meets the specified predicate
		template <typename LT_FN, typename PRED_FN>
		ssap_and_prc_cref_opt first_result_if(const ssaps_and_prcs_of_query &prm_ssaps_and_prcs, ///< The SSAP and PRC results to query
		                                      LT_FN                          prm_less_than,      ///< The less-than function to determine which results should come first ("lower" comes earlier)
		                                      PRED_FN                        prm_pred            ///< The predicate to specify acceptable results
		                                      ) {
			return detail::first_result_if_impl(
				prm_ssaps_and_prcs,
				prm_less_than,
				prm_pred
			);
		}
	} // namespace homcheck
} // namespace cath


namespace cath {
	namespace file {


		/// \brief Return the first SSAP result ordered by the specified less-than function that meets the specified predicate
		template <typename LT_FN, typename PRED_FN>
		ssap_scores_entry_cref_opt first_result_if(const ssap_scores_entry_vec &prm_ssaps,     ///< The SSAP results to query
		                                           LT_FN                        prm_less_than, ///< The less-than function to determine which results should come first ("lower" comes earlier)
		                                           PRED_FN                      prm_pred       ///< The predicate to specify acceptable results
		                                           ) {
			return homcheck::detail::first_result_if_impl(
				prm_ssaps,
				prm_less_than,
				prm_pred
			);
		}

		/// \brief Return the first N PRC results ordered by the specified less-than function that meet the specified predicate
		template <typename LT_FN, typename PRED_FN>
		prc_scores_entry_vec first_n_results_if(const prc_scores_entry_vec &prm_prcs,      ///< The PRC results to query
		                                        LT_FN                       prm_less_than, ///< The less-than function to determine which results should come first ("lower" comes earlier)
		                                        PRED_FN                     prm_pred,      ///< The predicate to specify acceptable results
		                                        const size_t               &prm_n          ///< The maximum number of results to return
		                                        ) {
			return homcheck::detail::first_n_results_if_impl(
				prm_prcs,
				prm_less_than,
				prm_pred,
				prm_n
			);
		}

		/// \brief Return the first N SSAP results ordered by the specified less-than function that meet the specified predicate
		template <typename LT_FN, typename PRED_FN>
		ssap_scores_entry_vec first_n_results_if(const ssap_scores_entry_vec &prm_ssaps,     ///< The SSAP results to query
		                                         LT_FN                        prm_less_than, ///< The less-than function to determine which results should come first ("lower" comes earlier)
		                                         PRED_FN                      prm_pred,      ///< The predicate to specify acceptable results
		                                         const size_t                &prm_n          ///< The maximum number of results to return
		                                         ) {

			return homcheck::detail::first_n_results_if_impl(
				prm_ssaps,
				prm_less_than,
				prm_pred,
				prm_n
			);
		}


	} // namespace file
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_HOMCHECK_TOOLS_FIRST_RESULT_IF_HPP
