/// \file
/// \brief The scored_hit_arch class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_SCORED_HIT_ARCH_HPP
#define _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_SCORED_HIT_ARCH_HPP

#include "cath/resolve_hits/algo/scored_arch_proxy.hpp"
#include "cath/resolve_hits/hit_arch.hpp"

namespace cath { namespace rslv { class calc_hit_list; } }

namespace cath {
	namespace rslv {

		/// \brief A hit architecture (ie hit_arch) along with the score associated with it
		class scored_hit_arch final {
		private:
			/// \brief The score associated with the architecture
			resscr_t the_score = INIT_SCORE;

			/// \brief The architecture
			hit_arch the_arch;

		public:
			scored_hit_arch(const resscr_t &,
			                hit_arch);
			scored_hit_arch() = default;

			const resscr_t & get_score() const noexcept;
			const hit_arch & get_arch() const noexcept;

			scored_hit_arch & operator+=(const calc_hit &);
			scored_hit_arch & operator+=(const scored_hit_arch &);

			scored_hit_arch & operator-=(const calc_hit &);
			scored_hit_arch & operator-=(const calc_hit_vec &);
		};

		/// \brief Ctor
		inline scored_hit_arch::scored_hit_arch(const resscr_t &prm_score, ///< The score associated with the architecture
		                                        hit_arch        prm_arch   ///< The architecture
		                                        ) : the_score ( prm_score             ),
		                                            the_arch  ( std::move( prm_arch ) ) {
		}

		/// \brief Getter for the score
		inline const resscr_t & scored_hit_arch::get_score() const noexcept {
			return the_score;
		}

		/// \brief Getter for the architecture
		inline const hit_arch & scored_hit_arch::get_arch() const noexcept {
			return the_arch;
		}

		/// \brief Add the specified calc_hit including its score) to this scored_hit_arch
		inline scored_hit_arch & scored_hit_arch::operator+=(const calc_hit &prm_hit ///< The calc_hit to add
		                                                     ) {
			// Do score second so that this can propagate any hit_arch::operator+=(const calc_hit &) exception guarantee
			the_arch  += prm_hit;
			the_score += prm_hit.get_score();
			return *this;
		}

		/// \brief Add the specified scored_hit_arch (including its score) to this scored_hit_arch
		inline scored_hit_arch & scored_hit_arch::operator+=(const scored_hit_arch &prm_scored_hit_arch ///< The scored_hit_arch to add
		                                                     ) {
			// Do score second so that this can propagate any hit_arch::operator+=(const hit_arch &) exception guarantee 
			the_arch  += prm_scored_hit_arch.get_arch();
			the_score += prm_scored_hit_arch.get_score();
			return *this;
		}

		/// \brief Remove the specified calc_hit (including its score) from this scored_hit_arch
		inline scored_hit_arch & scored_hit_arch::operator-=(const calc_hit &prm_hit ///< The calc_hit to remove
		                                                     ) {
			if ( the_arch.remove( prm_hit ) ) {
				the_score -= prm_hit.get_score();
			}
			return *this;
		}

		/// \brief Remove the specified list of hits (including their scores) from this scored_hit_arch
		inline scored_hit_arch & scored_hit_arch::operator-=(const calc_hit_vec &prm_hit_vec ///< The hits to remove
		                                                     ) {
			for (const calc_hit &the_hit : prm_hit_vec) {
				(*this) -= the_hit;
			}
			return *this;
		}

		/// \brief Add the specified calc_hit (including its score) to a copy of the specified scored_hit_arch
		///
		/// \relates scored_hit_arch
		inline scored_hit_arch operator+(scored_hit_arch  prm_scored_hit_arch, ///< The scored_hit_arch to copy and then add the specified hit
		                                 const calc_hit  &prm_hit              ///< The calc_hit to add
		                                 ) {
			prm_scored_hit_arch += prm_hit;
			return prm_scored_hit_arch;
		}

		/// \brief Add the specified calc_hit (including its score) to a copy of the specified scored_hit_arch
		///
		/// \relates scored_hit_arch
		inline scored_hit_arch operator+(const calc_hit  &prm_hit,            ///< The calc_hit to add
		                                 scored_hit_arch  prm_scored_hit_arch ///< The scored_hit_arch to copy and then add the specified hit
		                                 ) {
			prm_scored_hit_arch += prm_hit;
			return prm_scored_hit_arch;
		}

		/// \brief Add the second specified scored_hit_arch (including its scores) to a copy of the first
		///
		/// \relates scored_hit_arch
		inline scored_hit_arch operator+(scored_hit_arch        prm_scored_hit_arch_lhs, ///< The scored_hit_arch to copy and then add the specified scored_hit_arch
		                                 const scored_hit_arch &prm_scored_hit_arch_rhs  ///< The scored_hit_arch to add
		                                 ) {
			prm_scored_hit_arch_lhs += prm_scored_hit_arch_rhs;
			return prm_scored_hit_arch_lhs;
		}

		/// \brief Remove the specified calc_hit (including its score) from a copy of the specified scored_hit_arch
		///
		/// \relates scored_hit_arch
		inline scored_hit_arch operator-(scored_hit_arch  prm_scored_hit_arch, ///< The scored_hit_arch to copy and then remove the specified hit
		                                 const calc_hit  &prm_hit              ///< The calc_hit to remove
		                                 ) {
			prm_scored_hit_arch -= prm_hit;
			return prm_scored_hit_arch;
		}

		/// \brief Remove the specified list of hits (including their scores) from a copy of the specified scored_hit_arch
		///
		/// \relates scored_hit_arch
		inline scored_hit_arch operator-(scored_hit_arch     prm_scored_hit_arch, ///< The scored_hit_arch to copy and then remove the specified hits
		                                 const calc_hit_vec &prm_hit_vec          ///< The hits to remove
		                                 ) {
			prm_scored_hit_arch -= prm_hit_vec;
			return prm_scored_hit_arch;
		}

		scored_hit_arch make_scored_hit_arch(const scored_arch_proxy &,
		                                     const calc_hit_list &);

		full_hit_list get_full_hits_of_hit_arch(const scored_hit_arch &,
		                                        const full_hit_list &);

	} // namespace rslv
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_SCORED_HIT_ARCH_HPP
