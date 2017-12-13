/// \file
/// \brief The crh_score_spec class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_OPTIONS_SPEC_CRH_SCORE_SPEC_HPP
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_OPTIONS_SPEC_CRH_SCORE_SPEC_HPP

#include "resolve_hits/resolve_hits_type_aliases.hpp"

namespace cath {
	namespace rslv {

		/// \brief Specify how the scores should be handled in cath-resolve-hits
		class crh_score_spec final {
		private:
			/// \brief The degree to which long domains are preferred
			resscr_t long_domains_preference = DEFAULT_LONG_DOMAINS_PREFERENCE;

			/// \brief The degree to which high scores are preferred
			resscr_t high_scores_preference  = DEFAULT_HIGH_SCORES_PREFERENCE;

			/// \brief Whether to apply rules specific to CATH-Gene3D
			bool     apply_cath_rules        = DEFAULT_APPLY_CATH_RULES;

			/// \brief Whether to use a naive, greedy approach to resolving
			bool     naive_greedy            = DEFAULT_NAIVE_GREEDY;


		public:
			/// \brief The default value for the degree to which long domains are preferred
			static constexpr resscr_t DEFAULT_LONG_DOMAINS_PREFERENCE    =  0.0;

			/// \brief The default value for the degree to which high scores are preferred
			//
			// // For Jon, want x such that `( sqrt( 2 ) ) ^ x = 3`, which gives `x = ln( 3 ) / ln ( sqrt( 2 ) ) = 3.16992500144`
			// // \todo Come constexpr sqrt() and ln() or log(), make this be calculated at compile time
			// //       based on a constexpr power_base value elsewhere in the code
			// static constexpr resscr_t DEFAULT_HIGH_SCORES_PREFERENCE     = static_cast<resscr_t>( 3.17 );
			static constexpr resscr_t DEFAULT_HIGH_SCORES_PREFERENCE     = 0.0;

			/// \brief The default value for whether to apply rules specific to CATH-Gene3D
			static constexpr bool     DEFAULT_APPLY_CATH_RULES           = false;

			/// \brief The default value for whether to use a naive, greedy approach to resolving
			static constexpr bool     DEFAULT_NAIVE_GREEDY               = false;

			crh_score_spec() noexcept = default;

			explicit crh_score_spec(const bool &,
			                        const resscr_t & = DEFAULT_LONG_DOMAINS_PREFERENCE,
			                        const resscr_t & = DEFAULT_HIGH_SCORES_PREFERENCE,
			                        const bool & = DEFAULT_NAIVE_GREEDY);

			const resscr_t & get_long_domains_preference() const;
			const resscr_t & get_high_scores_preference() const;
			const bool & get_apply_cath_rules() const;
			const bool & get_naive_greedy() const;

			crh_score_spec & set_long_domains_preference(const resscr_t &);
			crh_score_spec & set_high_scores_preference(const resscr_t &);
			crh_score_spec & set_apply_cath_rules(const bool &);
			crh_score_spec & set_naive_greedy(const bool &);
		};

		crh_score_spec make_neutral_score_spec();

	} // namespace rslv
} // namespace cath

#endif
