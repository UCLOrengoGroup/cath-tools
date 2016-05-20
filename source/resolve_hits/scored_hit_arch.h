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

#ifndef SCORED_HIT_ARCH_H_INCLUDED
#define SCORED_HIT_ARCH_H_INCLUDED

#include "resolve_hits/hit_arch.h"
#include "resolve_hits/scored_arch_proxy.h"

namespace cath { namespace rslv { class hit_list; } }

namespace cath {
	namespace rslv {

		/// \brief TODOCUMENT
		class scored_hit_arch final {
		private:
			/// \brief TODOCUMENT
			resscr_t the_score = INIT_SCORE;

			/// \brief TODOCUMENT
			hit_arch the_arch;

		public:
			scored_hit_arch(const resscr_t &,
			                const hit_arch &);
			scored_hit_arch() = default;

			const resscr_t & get_score() const noexcept;
			const hit_arch & get_arch() const noexcept;

			scored_hit_arch & operator+=(const hit &);
			scored_hit_arch & operator+=(const scored_hit_arch &);

			scored_hit_arch & operator-=(const hit &);
			scored_hit_arch & operator-=(const hit_vec &);
		};

		/// \brief TODOCUMENT
		inline scored_hit_arch::scored_hit_arch(const resscr_t &arg_score, ///< TODOCUMENT
		                                        const hit_arch &arg_arch   ///< TODOCUMENT
		                                        ) : the_score ( arg_score ),
		                                            the_arch  ( arg_arch  ) {
		}

		/// \brief TODOCUMENT
		inline const resscr_t & scored_hit_arch::get_score() const noexcept {
			return the_score;
		}

		/// \brief TODOCUMENT
		inline const hit_arch & scored_hit_arch::get_arch() const noexcept {
			return the_arch;
		}

		/// \brief TODOCUMENT
		inline scored_hit_arch & scored_hit_arch::operator+=(const hit &arg_hit ///< TODOCUMENT
		                                                     ) {
			// Do score second so that this can propagate any hit_arch::operator+=(const hit &) exception guarantee 
			the_arch  += arg_hit;
			the_score += arg_hit.get_score();
			return *this;
		}

		/// \brief TODOCUMENT
		inline scored_hit_arch & scored_hit_arch::operator+=(const scored_hit_arch &arg_scored_hit_arch ///< TODOCUMENT
		                                                     ) {
			// Do score second so that this can propagate any hit_arch::operator+=(const hit_arch &) exception guarantee 
			the_arch  += arg_scored_hit_arch.get_arch();
			the_score += arg_scored_hit_arch.get_score();
			return *this;
		}

		/// \brief TODOCUMENT
		inline scored_hit_arch & scored_hit_arch::operator-=(const hit &arg_hit ///< TODOCUMENT
		                                                     ) {
			if ( the_arch.remove( arg_hit ) ) {
				the_score -= arg_hit.get_score();
			}
			return *this;
		}

		/// \brief TODOCUMENT
		inline scored_hit_arch & scored_hit_arch::operator-=(const hit_vec &arg_hit_vec ///< TODOCUMENT
		                                                     ) {
			for (const hit &the_hit : arg_hit_vec) {
				(*this) -= the_hit;
			}
			return *this;
		}

		/// \brief TODOCUMENT
		///
		/// \relates scored_hit_arch
		inline scored_hit_arch operator+(scored_hit_arch  arg_scored_hit_arch, ///< TODOCUMENT
		                                 const hit       &arg_hit              ///< TODOCUMENT
		                                 ) {
			arg_scored_hit_arch += arg_hit;
			return arg_scored_hit_arch;
		}

		/// \brief TODOCUMENT
		///
		/// \relates scored_hit_arch
		inline scored_hit_arch operator+(const hit       &arg_hit,            ///< TODOCUMENT
		                                 scored_hit_arch  arg_scored_hit_arch ///< TODOCUMENT
		                                 ) {
			arg_scored_hit_arch += arg_hit;
			return arg_scored_hit_arch;
		}

		/// \brief TODOCUMENT
		///
		/// \relates scored_hit_arch
		inline scored_hit_arch operator+(scored_hit_arch        arg_scored_hit_arch_lhs, ///< TODOCUMENT
		                                 const scored_hit_arch &arg_scored_hit_arch_rhs  ///< TODOCUMENT
		                                 ) {
			arg_scored_hit_arch_lhs += arg_scored_hit_arch_rhs;
			return arg_scored_hit_arch_lhs;
		}

		/// \brief TODOCUMENT
		///
		/// \relates scored_hit_arch
		inline scored_hit_arch operator-(scored_hit_arch  arg_scored_hit_arch, ///< TODOCUMENT
		                                 const hit       &arg_hit              ///< TODOCUMENT
		                                 ) {
			arg_scored_hit_arch -= arg_hit;
			return arg_scored_hit_arch;
		}

		/// \brief TODOCUMENT
		///
		/// \relates scored_hit_arch
		inline scored_hit_arch operator-(scored_hit_arch  arg_scored_hit_arch, ///< TODOCUMENT
		                                 const hit_vec   &arg_hit_vec          ///< TODOCUMENT
		                                 ) {
			arg_scored_hit_arch -= arg_hit_vec;
			return arg_scored_hit_arch;
		}

		/// \brief TODOCUMENT
		///
		/// \relates scored_hit_arch
		scored_hit_arch make_scored_hit_arch(const scored_arch_proxy &,
		                                     const hit_list &);

		// std::string to_string(const scored_hit_arch &);
		// std::ostream & operator<<(std::ostream &,
		//                           const scored_hit_arch &);

	}
}

#endif
