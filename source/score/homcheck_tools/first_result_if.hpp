/// \file
/// \brief The first_result_if header

#ifndef _CATH_TOOLS_SOURCE_SCORE_HOMCHECK_TOOLS_FIRST_RESULT_IF_H
#define _CATH_TOOLS_SOURCE_SCORE_HOMCHECK_TOOLS_FIRST_RESULT_IF_H

#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/irange.hpp>
#include <boost/optional.hpp>

#include "common/algorithm/copy_build.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/adaptor/limited.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "common/size_t_literal.hpp"
#include "common/type_aliases.hpp"
#include "file/file_type_aliases.hpp"
#include "score/homcheck_tools/ssaps_and_prcs_of_query.hpp"

#include <vector>

using namespace cath::common::literals;

namespace cath { namespace homcheck { class ssaps_and_prcs_of_query; } }

namespace cath {
	namespace homcheck {
		namespace detail {


			/// \brief Return the first result ordered by the specified less-than function that meets the specified predicate
			template <typename RES, typename LT_FN, typename PRED_FN>
			auto first_result_if_impl(const RES &arg_results,   ///< The results to query
			                          LT_FN      arg_less_than, ///< The less-than function to determine which results should come first ("lower" comes earlier)
			                          PRED_FN    arg_pred       ///< The predicate to specify acceptable results
			                          ) -> decltype( boost::make_optional( std::cref( arg_results[ 0 ] ) ) ) {
				auto indices = common::copy_build<size_vec>( boost::irange( 0_z, arg_results.size() ) );
				boost::range::sort(
					indices,
					[&] (const size_t &x, const size_t &y) {
						return arg_less_than( arg_results[ x ], arg_results[ y ] );
					}
				);
				const auto find_itr = boost::range::find_if(
					indices,
					[&] (const size_t &x) {
						return arg_pred( arg_results[ x ] );
					}
				);
				return ( find_itr != common::cend( indices ) ) ? boost::make_optional( std::cref( arg_results[ *find_itr ] ) )
				                                               : boost::none;
			}

			/// \brief Return the first N results ordered by the specified less-than function that meet the specified predicate
			template <typename RES, typename LT_FN, typename PRED_FN>
			auto first_n_results_if_impl(const RES    &arg_results,   ///< The results to query
			                             LT_FN         arg_less_than, ///< The less-than function to determine which results should come first ("lower" comes earlier)
			                             PRED_FN       arg_pred,      ///< The predicate to specify acceptable results
			                             const size_t &arg_n          ///< The maximum number of results to return
			                             ) -> std::vector< std::remove_const_t< std::remove_reference_t< decltype( arg_results[ 0 ] ) > > > {
				using result_type        = std::remove_const_t< std::remove_reference_t< decltype( arg_results[ 0 ] ) > >;
				using result_vector_type = std::vector<result_type>;
				auto indices = common::copy_build<size_vec>( boost::irange( 0_z, arg_results.size() ) );
				boost::range::sort(
					indices,
					[&] (const size_t &x, const size_t &y) {
						return arg_less_than( arg_results[ x ], arg_results[ y ] );
					}
				);

				return common::transform_build<result_vector_type> (
					indices
						| boost::adaptors::filtered( [&] (const size_t &x) { return arg_pred( arg_results[ x ] ); } )
						| common::limited( arg_n ),
					[&] (const size_t &x) {
						return arg_results[ x ];
					}
				);
			}
		} // namespace detail

		/// \brief Return the first SSAP and PRC result ordered by the specified less-than function that meets the specified predicate
		template <typename LT_FN, typename PRED_FN>
		ssap_and_prc_cref_opt first_result_if(const ssaps_and_prcs_of_query &arg_ssaps_and_prcs, ///< The SSAP and PRC results to query
		                                      LT_FN                          arg_less_than,      ///< The less-than function to determine which results should come first ("lower" comes earlier)
		                                      PRED_FN                        arg_pred            ///< The predicate to specify acceptable results
		                                      ) {
			return detail::first_result_if_impl(
				arg_ssaps_and_prcs,
				arg_less_than,
				arg_pred
			);
		}
	} // namespace homcheck
} // namespace cath


namespace cath {
	namespace file {


		/// \brief Return the first SSAP result ordered by the specified less-than function that meets the specified predicate
		template <typename LT_FN, typename PRED_FN>
		ssap_scores_entry_cref_opt first_result_if(const ssap_scores_entry_vec &arg_ssaps,     ///< The SSAP results to query
		                                           LT_FN                        arg_less_than, ///< The less-than function to determine which results should come first ("lower" comes earlier)
		                                           PRED_FN                      arg_pred       ///< The predicate to specify acceptable results
		                                           ) {
			return homcheck::detail::first_result_if_impl(
				arg_ssaps,
				arg_less_than,
				arg_pred
			);
		}

		/// \brief Return the first N PRC results ordered by the specified less-than function that meet the specified predicate
		template <typename LT_FN, typename PRED_FN>
		prc_scores_entry_vec first_n_results_if(const prc_scores_entry_vec &arg_prcs,      ///< The PRC results to query
		                                        LT_FN                       arg_less_than, ///< The less-than function to determine which results should come first ("lower" comes earlier)
		                                        PRED_FN                     arg_pred,      ///< The predicate to specify acceptable results
		                                        const size_t               &arg_n          ///< The maximum number of results to return
		                                        ) {
			return homcheck::detail::first_n_results_if_impl(
				arg_prcs,
				arg_less_than,
				arg_pred,
				arg_n
			);
		}

		/// \brief Return the first N SSAP results ordered by the specified less-than function that meet the specified predicate
		template <typename LT_FN, typename PRED_FN>
		ssap_scores_entry_vec first_n_results_if(const ssap_scores_entry_vec &arg_ssaps,     ///< The SSAP results to query
		                                         LT_FN                        arg_less_than, ///< The less-than function to determine which results should come first ("lower" comes earlier)
		                                         PRED_FN                      arg_pred,      ///< The predicate to specify acceptable results
		                                         const size_t                &arg_n          ///< The maximum number of results to return
		                                         ) {

			return homcheck::detail::first_n_results_if_impl(
				arg_ssaps,
				arg_less_than,
				arg_pred,
				arg_n
			);
		}


	} // namespace file
} // namespace cath

#endif
