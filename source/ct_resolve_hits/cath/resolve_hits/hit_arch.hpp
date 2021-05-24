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

#ifndef _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_HIT_ARCH_HPP
#define _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_HIT_ARCH_HPP

#include <boost/range/algorithm/adjacent_find.hpp>
#include <boost/range/algorithm/sort.hpp>

#include "cath/common/chrono/duration_to_seconds_string.hpp"
#include "cath/resolve_hits/calc_hit.hpp"
#include "cath/resolve_hits/hit_output_format.hpp"
#include "cath/resolve_hits/options/spec/hit_boundary_output.hpp"
#include "cath/resolve_hits/trim/trim_spec.hpp"

#include <tuple>

// clang-format off
namespace cath::rslv { class full_hit_list; }
// clang-format on

namespace cath::rslv {

	/// \brief An architecture of non-overlapping calc_hits
	///
	/// \invariant The calc_hits are kept sorted in ascending order of start residue
	class hit_arch final {
	private:
		/// \brief The (non-overlapping) calc_hits that make up the architecture
		calc_hit_vec the_hits;

		static void sort_hit_vec(calc_hit_vec &);

		void sanity_check() const;

	public:
		/// \brief An iterator type alias to make this a range
		using iterator       = calc_hit_vec::iterator;

		/// \brief A const_iterator type alias to make this a range
		using const_iterator = calc_hit_vec::const_iterator;

		hit_arch() = default; ///< \brief Default ctor
		explicit hit_arch(calc_hit_vec);

		hit_arch            (const hit_arch &)          = default; ///< \brief Default copy-ctor
		hit_arch            (hit_arch &&     ) noexcept = default; ///< \brief Default move-ctor
		hit_arch & operator=(const hit_arch &)          = default; ///< \brief Default copy-assignment-operator
		hit_arch & operator=(hit_arch &&     ) noexcept = default; ///< \brief Default move-assignment-operator

		[[nodiscard]] size_t size() const;
		[[nodiscard]] bool   empty() const;

		const calc_hit & operator[](const size_t &) const;

		[[nodiscard]] const_iterator begin() const;
		[[nodiscard]] const_iterator end() const;

		bool remove(const calc_hit &);
		hit_arch & operator+=(const calc_hit &);
		hit_arch & operator+=(const hit_arch &);
	};

	full_hit_list get_full_hits_of_hit_arch(const hit_arch &,
	                                        const full_hit_list &);

	std::string to_output_string(const hit_arch &,
	                             const full_hit_list &,
	                             const crh_segment_spec &,
	                             const hit_output_format & = hit_output_format::CLASS,
	                             const std::string & = std::string{},
	                             const hit_boundary_output & = hit_boundary_output::ORIG);

	hit_arch operator+(hit_arch,
	                   const calc_hit &);
	hit_arch operator+(hit_arch,
	                   const hit_arch &);

	/// \brief In-place sort the specified vector of calc_hits by calc_hit::get_hit_start_less() (ie, by their starts)
	inline void hit_arch::sort_hit_vec(calc_hit_vec &prm_hit_vec ///< The vector of calc_hits to in-place sort
	                                   ) {
		boost::range::sort(
			prm_hit_vec,
			calc_hit::get_hit_start_less()
		);
	}

	/// \brief Check that there are no overlaps between any calc_hits, and throw if any are found
	///
	/// \pre The calc_hits must be sorted before any calls to this method because it assumes
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
			&are_overlapping
		);
		if ( overlap_itr != ::std::cend( the_hits ) ) {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot create hit_arch with overlapping domains"));
		}
	}

	/// \brief Ctor from a vector of calc_hits
	inline hit_arch::hit_arch(calc_hit_vec prm_hit_arch ///< The vector of calc_hits from which to construct the hit_arch. Must have no mutually overlapping calc_hits. Need not be pre-sorted.
	                          ) : the_hits( std::move( prm_hit_arch ) ) {
		sort_hit_vec( the_hits );
		sanity_check();
	}

	/// \brief Return the number of calc_hits in the architecture
	inline size_t hit_arch::size() const {
		return the_hits.size();
	}

	/// \brief Return whether there are zero calc_hits in this architecture
	inline bool hit_arch::empty() const {
		return the_hits.empty();
	}

	/// \brief Const subscript operator for accessing the calc_hit at the specified index (after sorting in ascending order of start residue)
	inline const calc_hit & hit_arch::operator[](const size_t &prm_index ///< The index of the calc_hit to return
	                                             ) const {
		return the_hits[ prm_index ];
	}

	/// \brief Standard const begin() operator to make hit_arch into a range
	inline auto hit_arch::begin() const -> const_iterator {
		return ::std::cbegin( the_hits );
	}

	/// \brief Standard const end() operator to make hit_arch into a range
	inline auto hit_arch::end() const -> const_iterator {
		return ::std::cend( the_hits );
	}

	/// \brief If the hit_arch contains the specified calc_hit, remove it and return true; otherwise, return false
	inline bool hit_arch::remove(const calc_hit &prm_hit ///< The calc_hit to remove from the hit_arch
	                             ) {
		// Do score second so that this can propagate any hit_arch::operator-=(const hit_arch &) exception guarantee
		const auto find_itr = boost::range::find( the_hits, prm_hit );
		if ( find_itr != ::std::cend( the_hits ) ) {
			the_hits.erase( find_itr );
			return true;
		}
		return false;
	}

	/// \brief Add the specified calc_hit to this hit_arch
	///
	/// \pre The specified calc_hit may not overlap with any of the calc_hits contained within the hit_arch
	inline hit_arch & hit_arch::operator+=(const calc_hit &prm_hit ///< The calc_hit to add to this hit_arch
	                                       ) {
		the_hits.push_back( prm_hit );
		sort_hit_vec( the_hits );
		sanity_check();
		return *this;
	}

	/// \brief Add the specified hit_arch's calc_hits to this hit_arch
	///
	/// \pre The specified hit_arch's calc_hits may not overlap with any of the calc_hits contained within the hit_arch
	inline hit_arch & hit_arch::operator+=(const hit_arch &prm_hit_arch ///< The hit_arch whose calc_hits should be added to this hit_arch
	                                       ) {
		for (const calc_hit &the_hit : prm_hit_arch) {
			( *this ) += the_hit;
		}
		return *this;
	}

	/// \brief Add the specified calc_hit to the specified hit_arch
	///
	/// \pre The specified calc_hit may not overlap with any of the calc_hits contained within the hit_arch
	///
	/// The first hit_arch is taken by non-const value to avoid copying an rvalue argument
	///
	/// \relates hit_arch
	inline hit_arch operator+(hit_arch        prm_hit_arch, ///< The hit_arch to which the calc_hit should be added
	                          const calc_hit &prm_hit       ///< The calc_hit to add to the hit_arch
	                          ) {
		prm_hit_arch += prm_hit;
		return prm_hit_arch;
	}

	/// \brief Add the second specified hit_arch's calc_hits to the first specified hit_arch
	///
	/// The first hit_arch is taken by non-const value to avoid copying an rvalue argument
	///
	/// \relates hit_arch
	inline hit_arch operator+(hit_arch        prm_hit_arch_lhs, ///< The hit_arch to which the calc_hits should be added
	                          const hit_arch &prm_hit_arch_rhs  ///< The hit_arch whose calc_hits should be added to the other hit_arch
	                          ) {
		prm_hit_arch_lhs += prm_hit_arch_rhs;
		return prm_hit_arch_lhs;
	}

} // namespace cath::rslv

#endif // _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_HIT_ARCH_HPP
