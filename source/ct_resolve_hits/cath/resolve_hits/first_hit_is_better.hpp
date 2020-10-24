/// \file
/// \brief The first_hit_is_better header

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

#ifndef _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_FIRST_HIT_IS_BETTER_HPP
#define _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_FIRST_HIT_IS_BETTER_HPP

#include <boost/logic/tribool.hpp>

#include "cath/common/boost_addenda/tribool/tribool.hpp"
#include "cath/resolve_hits/calc_hit.hpp"
#include "cath/resolve_hits/full_hit.hpp"
#include "cath/resolve_hits/full_hit_list.hpp"

namespace cath {
	namespace rslv {

		/// \brief Return whether the first specified hit is better-than or worse-than the second specified hit
		///
		/// "a is better than b" means b can be dropped from resolving calculations involving a
		/// because the a is a better choice. In practice this means that the a's score is ≥ b's
		/// and a's segments are a subset of b's.
		///
		/// When scores and segments are equal, the labels and then label indices are compared
		/// because one of the two may as well be discarded to save later computation.
		///
		/// If the labels and label indices are equal, then the result is indeterminate.
		///
		/// \returns tribool{true} for better than; tribool{false} for worse-than; indeterminate otherwise
		///
		/// This induces a partial ordering because:
		///  * x not-better-than x
		///  * x better-than y => y not-better-than x
		///  * x better-than y && y better-than z => x better-than z
		/// ...but not a strict weak ordering because, for example:
		///  * [10-19, 1.0] equivalent-to [ 5-15, 1.0] and
		///  * [ 5-15, 1.0] equivalent-to [10-20, 1.0] but
		///  * [10-19, 1.0] better-than   [10-20, 1.0]
		///
		/// \todo Would make sense to return spaceship numbers, a la Herb Sutter's comparison proposal p0515r0,
		///       rather than a tribool
		///
		/// \relates calc_hit
		///
		/// \relatesalso full_hit_list
		inline boost::logic::tribool first_hit_is_better(const calc_hit      &prm_lhs,      ///< The first hit to compare
		                                                 const calc_hit      &prm_rhs,      ///< The second hit to compare
		                                                 const full_hit_list &prm_full_hits ///< The full hits, which are need for comparing labels for otherwise very similar hits
		                                                 ) {
			// If the neither of the hits covers the other, than neither can be better than the other
			if ( ! one_covers_other( prm_lhs, prm_rhs ) ) {
				return boost::logic::indeterminate;
			}

			// If the first's score is better...
			if ( prm_lhs.get_score() > prm_rhs.get_score() ) {
				// If the first is (non-strictly) within the second, it's better; else neither is better
				return first_is_not_outside_second( prm_lhs, prm_rhs )
					? boost::logic::tribool{ true }
					: boost::logic::indeterminate;
			}

			// If the second's score is better...
			if ( prm_lhs.get_score() < prm_rhs.get_score() ) {
				// If the second is (non-strictly) within the first, it's better; else neither is better
				return first_is_not_outside_second( prm_rhs, prm_lhs )
					? boost::logic::tribool{ false }
					: boost::logic::indeterminate;
			}

			// Otherwise, scores are equal so...

			// If the first hit is strictly within the second, it's better; else...
			if ( first_is_shorter_and_within_second( prm_lhs, prm_rhs ) ) {
				return boost::logic::tribool{ true };
			}
			// If the second hit is strictly within the first, it's better; else...
			if ( first_is_shorter_and_within_second( prm_rhs, prm_lhs ) ) {
				return boost::logic::tribool{ false };
			}

			// Otherwise, both score and segments are equal so...

			/// Compare labels and then label indices
			const hitidx_t    &idx_lhs   = prm_lhs.get_label_idx();
			const hitidx_t    &idx_rhs   = prm_rhs.get_label_idx();
			const std::string &label_lhs = prm_full_hits[ idx_lhs ].get_label();
			const std::string &label_rhs = prm_full_hits[ idx_rhs ].get_label();
			return ( std::tie( label_lhs, idx_lhs ) < std::tie( label_rhs, idx_rhs ) ) ? boost::logic::tribool{ true  } :
			       ( std::tie( label_lhs, idx_lhs ) > std::tie( label_rhs, idx_rhs ) ) ? boost::logic::tribool{ false } :
			                                                                             boost::logic::indeterminate;
		}

		/// \brief Return whether the first specified hit is better-than or worse-than the second specified hit
		///
		/// "a is better than b" means b can be dropped from resolving calculations involving a
		/// because the a is a better choice. In practice this means that the a's score is ≥ b's
		/// and a's segments are a subset of b's.
		///
		/// When scores and segments are equal, the labels and then label indices are compared
		/// because one of the two may as well be discarded to save later computation.
		///
		/// If the labels and label indices are equal, then the result is indeterminate.
		///
		/// \returns tribool{true} for better than; tribool{false} for worse-than; indeterminate otherwise
		///
		/// This induces a partial ordering because:
		///  * x not-better-than x
		///  * x better-than y => y not-better-than x
		///  * x better-than y && y better-than z => x better-than z
		/// ...but not a strict weak ordering because, for example:
		///  * [10-19, 1.0] equivalent-to [ 5-15, 1.0] and
		///  * [ 5-15, 1.0] equivalent-to [10-20, 1.0] but
		///  * [10-19, 1.0] better-than   [10-20, 1.0]
		///
		/// \todo Would make sense to return spaceship numbers, a la Herb Sutter's comparison proposal p0515r0,
		///       rather than a tribool
		///
		/// \relates calc_hit
		///
		/// \relatesalso full_hit_list
		inline boost::logic::tribool first_hit_is_better(const full_hit &prm_lhs, ///< The first hit to compare
		                                                 const full_hit &prm_rhs  ///< The second hit to compare
		                                                 ) {
			// If the neither of the hits covers the other, than neither can be better than the other
			if ( ! one_covers_other( prm_lhs, prm_rhs ) ) {
				return boost::logic::indeterminate;
			}

			// If the first's score is better...
			if ( prm_lhs.get_score() > prm_rhs.get_score() ) {
				// If the first is (non-strictly) within the second, it's better; else neither is better
				return first_is_not_outside_second( prm_lhs, prm_rhs )
					? boost::logic::tribool{ true }
					: boost::logic::indeterminate;
			}

			// If the second's score is better...
			if ( prm_lhs.get_score() < prm_rhs.get_score() ) {
				// If the second is (non-strictly) within the first, it's better; else neither is better
				return first_is_not_outside_second( prm_rhs, prm_lhs )
					? boost::logic::tribool{ false }
					: boost::logic::indeterminate;
			}

			// Otherwise, scores are equal so...

			// If the first hit is strictly within the second, it's better; else...
			if ( first_is_shorter_and_within_second( prm_lhs, prm_rhs ) ) {
				return boost::logic::tribool{ true };
			}
			// If the second hit is strictly within the first, it's better; else...
			if ( first_is_shorter_and_within_second( prm_rhs, prm_lhs ) ) {
				return boost::logic::tribool{ false };
			}

			// Otherwise, both score and segments are equal so...

			/// Compare labels and then label indices
			const std::string &label_lhs = prm_lhs.get_label();
			const std::string &label_rhs = prm_rhs.get_label();
			return ( label_lhs < label_rhs ) ? boost::logic::tribool{ true  } :
			       ( label_lhs > label_rhs ) ? boost::logic::tribool{ false } :
			                                   boost::logic::indeterminate;
		}

	} // namespace rslv
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_FIRST_HIT_IS_BETTER_HPP
