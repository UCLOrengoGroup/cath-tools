/// \file
/// \brief The hit_arch class header

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

#ifndef HIT_ARCH_H_INCLUDED
#define HIT_ARCH_H_INCLUDED

#include <boost/log/trivial.hpp>
#include <boost/range/algorithm/adjacent_find.hpp>
#include <boost/range/algorithm/sort.hpp>

#include "common/c++14/cbegin_cend.h"
#include "common/chrono/duration_to_seconds_string.h"
#include "resolve_hits/hit.h"

#include <tuple>

namespace cath {
	namespace rslv {

		/// \brief TODOCUMENT
		class hit_arch final {
		private:
			/// \brief TODOCUMENT
			hit_vec the_hits;

			static void sort_hit_vec(hit_vec &);

			void sanity_check() const;

		public:
			using iterator       = hit_vec::iterator;
			using const_iterator = hit_vec::const_iterator;

			explicit hit_arch(const hit_vec &);
			hit_arch() = default;

			hit_arch(const hit_arch &arg_rhs
			         ) : the_hits( arg_rhs.the_hits ) {
			}

			hit_arch(hit_arch &&) noexcept = default;
			hit_arch & operator=(const hit_arch &) = default;
			hit_arch & operator=(hit_arch &&) noexcept = default;

			size_t size() const;
			bool empty() const;

			const hit & operator[](const size_t &) const;

			const_iterator begin() const;
			const_iterator end() const;

			bool remove(const hit &);
			hit_arch & operator+=(const hit &);
			hit_arch & operator+=(const hit_arch &);
		};

		std::string to_string(const hit_arch &,
		                      const hit_output_format & = hit_output_format::CLASS);
		std::ostream & operator<<(std::ostream &,
		                          const hit_arch &);
		hit_arch operator+(hit_arch,
		                   const hit &);
		hit_arch operator+(hit_arch,
		                   const hit_arch &);

		/// \brief TODOCUMENT
		inline void hit_arch::sort_hit_vec(hit_vec &arg_hit_vec ///< TODOCUMENT
		                                   ) {
			boost::range::sort(
				arg_hit_vec,
				hit::get_hit_start_less()
			);
		}

		inline void hit_arch::sanity_check() const {
			const auto overlap_itr = boost::range::adjacent_find(
				the_hits,
				&hits_overlap
			);
			if ( overlap_itr != common::cend( the_hits ) ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot create hit_arch with overlapping domains"));
			}
		}

		/// \brief TODOCUMENT
		inline hit_arch::hit_arch(const hit_vec &arg_hit_arch ///< TODOCUMENT
		                          ) : the_hits( arg_hit_arch ) {
			sort_hit_vec( the_hits );
			sanity_check();
		}

		/// \brief TODOCUMENT
		inline size_t hit_arch::size() const {
			return the_hits.size();
		}

		/// \brief TODOCUMENT
		inline bool hit_arch::empty() const {
			return the_hits.empty();
		}

		/// \brief TODOCUMENT
		inline const hit & hit_arch::operator[](const size_t &arg_index ///< TODOCUMENT
		                                        ) const {
			return the_hits[ arg_index ];
		}

		/// \brief TODOCUMENT
		inline auto hit_arch::begin() const -> const_iterator {
			return common::cbegin( the_hits );
		}

		/// \brief TODOCUMENT
		inline auto hit_arch::end() const -> const_iterator {
			return common::cend( the_hits );
		}

		/// \brief TODOCUMENT
		inline bool hit_arch::remove(const hit &arg_hit ///< TODOCUMENT
		                             ) {
			// Do score second so that this can propagate any hit_arch::operator-=(const hit_arch &) exception guarantee
			const auto find_itr = boost::range::find( the_hits, arg_hit );
			if ( find_itr != common::cend( the_hits ) ) {
				the_hits.erase( find_itr );
				return true;
			}
			return false;
		}

		/// \brief TODOCUMENT
		inline hit_arch & hit_arch::operator+=(const hit &arg_hit ///< TODOCUMENT
		                                       ) {
			the_hits.push_back( arg_hit );
			sort_hit_vec( the_hits );
			sanity_check();
			return *this;
		}

		/// \brief TODOCUMENT
		inline hit_arch & hit_arch::operator+=(const hit_arch &arg_hit_arch ///< TODOCUMENT
		                                       ) {
			for (const hit &the_hit : arg_hit_arch) {
				( *this ) += the_hit;
			}
			return *this;
		}

		/// \brief TODOCUMENT
		inline hit_arch operator+(hit_arch   arg_hit_arch, ///< TODOCUMENT
		                          const hit &arg_hit       ///< TODOCUMENT
		                          ) {
			arg_hit_arch += arg_hit;
			return arg_hit_arch;
		}

		/// \brief TODOCUMENT
		inline hit_arch operator+(hit_arch        arg_hit_arch_lhs, ///< TODOCUMENT
		                          const hit_arch &arg_hit_arch_rhs  ///< TODOCUMENT
		                          ) {
			arg_hit_arch_lhs += arg_hit_arch_rhs;
			return arg_hit_arch_lhs;
		}


	}
}

#endif
