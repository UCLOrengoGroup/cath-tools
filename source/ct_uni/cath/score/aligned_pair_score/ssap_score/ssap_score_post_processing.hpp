/// \file
/// \brief The ssap_score_post_processing class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_SCORE_ALIGNED_PAIR_SCORE_SSAP_SCORE_SSAP_SCORE_POST_PROCESSING_HPP
#define _CATH_TOOLS_SOURCE_UNI_SCORE_ALIGNED_PAIR_SCORE_SSAP_SCORE_SSAP_SCORE_POST_PROCESSING_HPP

#include "cath/common/algorithm/constexpr_is_uniq.hpp"
#include "cath/common/cpp20/make_array.hpp"
#include "cath/common/detail/maybe_unused_namespace_scope_constexpr.hpp"

#include <array>
#include <iosfwd>
#include <map>

namespace cath {
	namespace score {

		/// \brief How to post process the basic score once it has been calculated
		///
		/// The normalisation steps all use the length getter, which is optionally specified in ssap_score's ctor.
		///
		/// All these strategies are followed by a multiplication by 100 to get a score out of 100
		enum class ssap_score_post_processing : char {
			SIMPLE_NORMLS,          ///< Perform a simple normalisation
			COMPLX_NORMLS,          ///< Perform a more intricate normalisation (which considers the number of residues not compared)
			SIMPLE_NORMLS_THEN_LOG, ///< Perform a simple normalisation and then (multiply by 1000 and) log the result
			COMPLX_NORMLS_THEN_LOG, ///< Perform a more intricate normalisation and then (multiply by 1000 and) log the result
			LOG_THEN_SIMPLE_NORMLS  ///< (Multiply by 1000 and) log the result and then perform a simple normalisation
		};

		/// \brief TODOCUMENT
		static constexpr auto all_ssap_score_post_processings = common::make_array(
			ssap_score_post_processing::SIMPLE_NORMLS,
			ssap_score_post_processing::COMPLX_NORMLS,
			ssap_score_post_processing::SIMPLE_NORMLS_THEN_LOG,
			ssap_score_post_processing::COMPLX_NORMLS_THEN_LOG,
			ssap_score_post_processing::LOG_THEN_SIMPLE_NORMLS
		);

		static_assert( common::constexpr_is_uniq( all_ssap_score_post_processings ), "all_ssap_score_post_processings shouldn't contain repeated values" );

		/// \brief TODOCUMENT
		static constexpr size_t num_ssap_score_post_processings = std::tuple_size< decltype( all_ssap_score_post_processings ) >::value;
		MAYBE_UNUSED_NAMESPACE_SCOPE_CONSTEXPR( num_ssap_score_post_processings )

		/// \brief TODOCUMENT
		struct name_of_ssap_score_post_processing final {
			static std::map<ssap_score_post_processing, std::string> get();
		};

		std::ostream & operator<<(std::ostream &,
		                          const ssap_score_post_processing &);


		bool has_post_log(const ssap_score_post_processing &);
		bool has_pre_log(const ssap_score_post_processing &);
		bool normalisation_is_simple(const ssap_score_post_processing &);

	} // namespace score
} // namespace cath

#endif
