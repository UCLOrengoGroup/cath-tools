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

#include "common/chrono/duration_to_seconds_string.h"
#include "common/cpp14/cbegin_cend.h"
#include "resolve_hits/hit.h"

#include <tuple>

namespace cath {
	namespace rslv {

		/// \brief An architecture of non-overlapping hits
		///
		/// \invariant The hits are kept sorted in ascending order of start residue
		class hit_arch final {
		private:
			/// \brief The (non-overlapping) hits that make up the architecture
			hit_vec the_hits;

			static void sort_hit_vec(hit_vec &);

			void sanity_check() const;

		public:
			/// \brief An iterator type alias to make this a range
			using iterator       = hit_vec::iterator;

			/// \brief A const_iterator type alias to make this a range
			using const_iterator = hit_vec::const_iterator;

			hit_arch() = default; ///< \brief Default ctor
			explicit hit_arch(const hit_vec &);

			hit_arch            (const hit_arch &)          = default; ///< \brief Default copy-ctor
			hit_arch            (hit_arch &&     ) noexcept = default; ///< \brief Default move-ctor
			hit_arch & operator=(const hit_arch &)          = default; ///< \brief Default copy-assignment-operator
			hit_arch & operator=(hit_arch &&     ) noexcept = default; ///< \brief Default move-assignment-operator

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
		                      const str_vec &,
		                      const hit_output_format & = hit_output_format::CLASS,
		                      const std::string & = std::string{});

		hit_arch operator+(hit_arch,
		                   const hit &);
		hit_arch operator+(hit_arch,
		                   const hit_arch &);

		/// \brief In-place sort the specified vector of hits by hit::get_hit_start_less() (ie, by their starts)
		inline void hit_arch::sort_hit_vec(hit_vec &arg_hit_vec ///< The vector of hits to in-place sort
		                                   ) {
			boost::range::sort(
				arg_hit_vec,
				hit::get_hit_start_less()
			);
		}

		/// \brief Check that there are no overlaps between any hits, and throw if any are found
		///
		/// \pre The hits must be sorted before any calls to this method because it assumes
		///      that any overlaps will be detectable in neighbours.
		///
		/// \todo Make this check more comprehensive so that it would no longer miss cases like:
		///       ~~~~~
		///       **   **
		///          **
		///            **
		///       ~~~~~
		///       This is low-priority because it seems very unlikely that the calling code creates
		///       overlaps that have never yet been detected
		inline void hit_arch::sanity_check() const {
			const auto overlap_itr = boost::range::adjacent_find(
				the_hits,
				&hits_overlap
			);
			if ( overlap_itr != common::cend( the_hits ) ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot create hit_arch with overlapping domains"));
			}
		}

		/// \brief Ctor from a vector of hits
		inline hit_arch::hit_arch(const hit_vec &arg_hit_arch ///< The vector of hits from which to construct the hit_arch. Must have no mutually overlapping hits. Need not be pre-sorted.
		                          ) : the_hits( arg_hit_arch ) {
			sort_hit_vec( the_hits );
			sanity_check();
		}

		/// \brief Return the number of hits in the architecture
		inline size_t hit_arch::size() const {
			return the_hits.size();
		}

		/// \brief Return whether there are zero hits in this architecture
		inline bool hit_arch::empty() const {
			return the_hits.empty();
		}

		/// \brief Const subscript operator for accessing the hit at the specified index (after sorting in ascending order of start residue)
		inline const hit & hit_arch::operator[](const size_t &arg_index ///< The index of the hit to return
		                                        ) const {
			return the_hits[ arg_index ];
		}

		/// \brief Standard const begin() operator to make hit_arch into a range
		inline auto hit_arch::begin() const -> const_iterator {
			return common::cbegin( the_hits );
		}

		/// \brief Standard const end() operator to make hit_arch into a range
		inline auto hit_arch::end() const -> const_iterator {
			return common::cend( the_hits );
		}

		/// \brief If the hit_arch contains the specified hit, remove it and return true; otherwise, return false
		inline bool hit_arch::remove(const hit &arg_hit ///< The hit to remove from the hit_arch
		                             ) {
			// Do score second so that this can propagate any hit_arch::operator-=(const hit_arch &) exception guarantee
			const auto find_itr = boost::range::find( the_hits, arg_hit );
			if ( find_itr != common::cend( the_hits ) ) {
				the_hits.erase( find_itr );
				return true;
			}
			return false;
		}

		/// \brief Add the specified hit to this hit_arch
		///
		/// \pre The specified hit may not overlap with any of the hits contained within the hit_arch
		inline hit_arch & hit_arch::operator+=(const hit &arg_hit ///< The hit to add to this hit_arch
		                                       ) {
			the_hits.push_back( arg_hit );
			sort_hit_vec( the_hits );
			sanity_check();
			return *this;
		}

		/// \brief Add the specified hit_arch's hits to this hit_arch
		///
		/// \pre The specified hit_arch's hits may not overlap with any of the hits contained within the hit_arch
		inline hit_arch & hit_arch::operator+=(const hit_arch &arg_hit_arch ///< The hit_arch whose hits should be added to this hit_arch
		                                       ) {
			for (const hit &the_hit : arg_hit_arch) {
				( *this ) += the_hit;
			}
			return *this;
		}

		/// \brief Add the specified hit to the specified hit_arch
		///
		/// \pre The specified hit may not overlap with any of the hits contained within the hit_arch
		///
		/// The first hit_arch is taken by non-const value to avoid copying an rvalue argument
		///
		/// \relates hit_arch
		inline hit_arch operator+(hit_arch   arg_hit_arch, ///< The hit_arch to which the hit should be added
		                          const hit &arg_hit       ///< The hit to add to the hit_arch
		                          ) {
			arg_hit_arch += arg_hit;
			return arg_hit_arch;
		}

		/// \brief Add the second specified hit_arch's hits to the first specified hit_arch
		///
		/// The first hit_arch is taken by non-const value to avoid copying an rvalue argument
		///
		/// \relates hit_arch
		inline hit_arch operator+(hit_arch        arg_hit_arch_lhs, ///< The hit_arch to which the hits should be added
		                          const hit_arch &arg_hit_arch_rhs  ///< The hit_arch whose hits should be added to the other hit_arch
		                          ) {
			arg_hit_arch_lhs += arg_hit_arch_rhs;
			return arg_hit_arch_lhs;
		}


	}
}

#endif
