/// \file
/// \brief The scored_arch_proxy class header

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

#ifndef SCORED_ARCH_PROXY_H_INCLUDED
#define SCORED_ARCH_PROXY_H_INCLUDED

#include "common/c++14/cbegin_cend.h"
#include "resolve_hits/resolve_hits_type_aliases.h"

namespace cath {
	namespace rslv {

		/// \brief TODOCUMENT
		class scored_arch_proxy final {
		private:
			/// \brief TODOCUMENT
			resscr_t the_score = INIT_SCORE;

			/// \brief TODOCUMENT
			hitidx_vec hit_indices;

		public:
			using const_iterator = hitidx_vec::const_iterator;

			scored_arch_proxy() = default;

			const resscr_t & get_score() const;

			bool empty() const;
			size_t size() const;
			const hitidx_t & operator[](const size_t &) const;

			const_iterator begin() const;
			const_iterator end() const;

			scored_arch_proxy & add_hit(const resscr_t &,
			                            const hitidx_t &);
			// scored_arch_proxy & operator+=(const resscr_t &,
			//                                const hitidx_vec);
			// scored_arch_proxy & operator-=(const resscr_t &,
			//                                const hitidx_t);
			// scored_arch_proxy & operator-=(const resscr_t &,
			//                                const hitidx_vec);
			// scored_arch_proxy & operator+=(const scored_arch_proxy &);

			// scored_arch_proxy & operator-=(const hit &);
			// scored_arch_proxy & operator-=(const hit_vec &);
		};

		/// \brief TODOCUMENT
		inline const resscr_t & scored_arch_proxy::get_score() const {
			return the_score;
		}

		/// \brief TODOCUMENT
		inline bool scored_arch_proxy::empty() const {
			return hit_indices.empty();
		}

		/// \brief TODOCUMENT
		inline size_t scored_arch_proxy::size() const {
			return hit_indices.size();
		}

		/// \brief TODOCUMENT
		inline const hitidx_t & scored_arch_proxy::operator[](const size_t &arg_index ///< TODOCUMENT
		                                                      ) const {
			return hit_indices[ arg_index ];
		}

		/// \brief TODOCUMENT
		inline auto scored_arch_proxy::begin() const -> const_iterator {
			return common::cbegin( hit_indices );
		}

		/// \brief TODOCUMENT
		inline auto scored_arch_proxy::end() const -> const_iterator {
			return common::cend( hit_indices );
		}

		/// \brief TODOCUMENT
		inline scored_arch_proxy & scored_arch_proxy::add_hit(const resscr_t &arg_score,    ///< TODOCUMENT
		                                                      const hitidx_t &arg_hit_index ///< TODOCUMENT
		                                                      ) {
			the_score += arg_score;
			hit_indices.push_back( arg_hit_index );
			return *this;
		}

		/// \brief TODOCUMENT
		inline scored_arch_proxy add_hit_copy(scored_arch_proxy  arg_scored_arch_proxy,///< TODOCUMENT
		                                      const resscr_t    &arg_score,             ///< TODOCUMENT
		                                      const hitidx_t    &arg_hit_index          ///< TODOCUMENT
		                                      ) {
			arg_scored_arch_proxy.add_hit( arg_score, arg_hit_index );
			return arg_scored_arch_proxy;
		}

		// /// \brief TODOCUMENT
		// ///
		// /// \relates scored_arch_proxy
		// inline scored_arch_proxy operator+(scored_arch_proxy  arg_scored_arch_proxy, ///< TODOCUMENT
		//                                  const hit       &arg_hit              ///< TODOCUMENT
		//                                  ) {
		// 	arg_scored_arch_proxy += arg_hit;
		// 	return arg_scored_arch_proxy;
		// }

		// /// \brief TODOCUMENT
		// ///
		// /// \relates scored_arch_proxy
		// inline scored_arch_proxy operator+(const hit       &arg_hit,            ///< TODOCUMENT
		//                                  scored_arch_proxy  arg_scored_arch_proxy ///< TODOCUMENT
		//                                  ) {
		// 	arg_scored_arch_proxy += arg_hit;
		// 	return arg_scored_arch_proxy;
		// }

		// /// \brief TODOCUMENT
		// ///
		// /// \relates scored_arch_proxy
		// inline scored_arch_proxy operator+(scored_arch_proxy        arg_scored_arch_proxy_lhs, ///< TODOCUMENT
		//                                  const scored_arch_proxy &arg_scored_arch_proxy_rhs  ///< TODOCUMENT
		//                                  ) {
		// 	arg_scored_arch_proxy_lhs += arg_scored_arch_proxy_rhs;
		// 	return arg_scored_arch_proxy_lhs;
		// }

		// /// \brief TODOCUMENT
		// ///
		// /// \relates scored_arch_proxy
		// inline scored_arch_proxy operator-(scored_arch_proxy  arg_scored_arch_proxy, ///< TODOCUMENT
		//                                  const hit       &arg_hit              ///< TODOCUMENT
		//                                  ) {
		// 	arg_scored_arch_proxy -= arg_hit;
		// 	return arg_scored_arch_proxy;
		// }

		// /// \brief TODOCUMENT
		// ///
		// /// \relates scored_arch_proxy
		// inline scored_arch_proxy operator-(scored_arch_proxy  arg_scored_arch_proxy, ///< TODOCUMENT
		//                                  const hit_vec   &arg_hit_vec          ///< TODOCUMENT
		//                                  ) {
		// 	arg_scored_arch_proxy -= arg_hit_vec;
		// 	return arg_scored_arch_proxy;
		// }


	}
}

#endif
