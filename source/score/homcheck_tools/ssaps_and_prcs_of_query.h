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
#include "file/file_type_aliases.h"
#include "score/homcheck_tools/ssap_and_prc.h"
#include "score/score_type_aliases.h"

#include <iosfwd>
#include <vector>

namespace cath { namespace file { class prc_scores_entry; } }
namespace cath { namespace file { class ssap_scores_entry; } }
namespace cath { namespace file { using prc_scores_entry_vec = std::vector<prc_scores_entry>; } }
namespace cath { namespace file { using ssap_scores_entry_vec = std::vector<ssap_scores_entry>; } }
namespace cath { namespace homcheck { class ssap_and_prc; } }
namespace cath { namespace homcheck { class superfamily_of_domain; } }

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

			bool empty() const;
			size_t size() const;

			const ssap_and_prc & operator[](const size_t &) const;

			const_iterator begin() const;
			const_iterator end() const;
		};

		const std::string & get_query_id(const ssaps_and_prcs_of_query &);

		ssap_and_prc_cref_opt best_magic_function_assignable(const ssaps_and_prcs_of_query &,
		                                                     const superfamily_of_domain &);

		file::ssap_scores_entry_cref_opt best_fold_level_match(const file::ssap_scores_entry_vec &,
		                                                       const superfamily_of_domain &);

		ssaps_and_prcs_of_query make_ssaps_and_prcs_of_query(const file::ssap_scores_entry_vec &,
		                                                     const file::prc_scores_entry_vec &);

	}
}

#endif
