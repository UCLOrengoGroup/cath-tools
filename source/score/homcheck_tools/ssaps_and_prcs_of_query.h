/// \file
/// \brief The ssaps_and_prcs_of_query class header

#ifndef SSAPS_AND_PRCS_OF_QUERY_H_INCLUDED
#define SSAPS_AND_PRCS_OF_QUERY_H_INCLUDED

#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/irange.hpp>
#include <boost/optional.hpp>

#include "common/algorithm/copy_build.h"
#include "common/c++14/cbegin_cend.h"
#include "common/size_t_literal.h"
#include "common/type_aliases.h"
#include "score/homcheck_tools/ssap_and_prc.h"
#include "score/score_type_aliases.h"

#include <iosfwd>
#include <vector>

namespace cath { namespace file { class prc_scores_entry; } }
namespace cath { namespace file { using prc_scores_entry_vec = std::vector<prc_scores_entry>; } }
namespace cath { namespace file { class ssap_scores_entry; } }
namespace cath { namespace file { using ssap_scores_entry_vec = std::vector<ssap_scores_entry>; } }
namespace cath { namespace homcheck { class ssap_and_prc; } }

using namespace cath::common::literals;

namespace cath {
	namespace homcheck {

		/// \brief Represent the ssap_and_prc results associated with a single query
		///
		/// \invariant All ssap_and_prc_entries have the same query_id
		class ssaps_and_prcs_of_query final {
		private:
			/// \brief The SSAP and PRC entries
			ssap_and_prc_vec ssap_and_prc_entries;

			void sanity_check() const;

		public:
			using const_iterator = ssap_and_prc_vec::const_iterator;

			ssaps_and_prcs_of_query() = default;
			ssaps_and_prcs_of_query(const ssap_and_prc_vec &);

			size_t size() const;

			const ssap_and_prc & operator[](const size_t &) const;

			const_iterator begin() const;
			const_iterator end() const;
		};

		boost::optional<std::reference_wrapper<const ssap_and_prc>> best_magic_function(const ssaps_and_prcs_of_query &);

		/// \brief TODOCUMENT
		template <typename FN>
		boost::optional<std::reference_wrapper<const ssap_and_prc>> best_magic_function_if(const ssaps_and_prcs_of_query &arg_ssaps_and_prcs, ///< TODOCUMENT
		                                                                                   FN                             arg_pred            ///< TODOCUMENT
		                                                                                   ) {
			auto indices = common::copy_build<size_vec>( boost::irange( 0_z, arg_ssaps_and_prcs.size() ) );
			boost::range::sort(
				indices,
				[&] (const size_t &x, const size_t &y) {
					// Reverse inequality to put the highest magic_function values to the start
					return arg_ssaps_and_prcs[ x ].get_magic_function_score() > arg_ssaps_and_prcs[ y ].get_magic_function_score();
				}
			);
			const auto find_itr = boost::range::find_if(
				indices,
				[&] (const size_t &x) {
					return arg_pred( arg_ssaps_and_prcs[ x ] );
				}
			);
			return ( find_itr != common::cend( indices ) ) ? boost::make_optional( std::cref( arg_ssaps_and_prcs[ *find_itr ] ) )
			                                               : boost::none;
		}

		ssaps_and_prcs_of_query make_ssaps_and_prcs_of_query(const file::ssap_scores_entry_vec &,
		                                                     const file::prc_scores_entry_vec &);

	}
}

#endif
