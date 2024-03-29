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

#ifndef CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO_SCORED_ARCH_PROXY_HPP
#define CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO_SCORED_ARCH_PROXY_HPP

#include "cath/resolve_hits/resolve_hits_type_aliases.hpp"

namespace cath::rslv {

	/// \brief Represent an architecture by the indices of its hits in the calc_hit_list, along with the associated score
	///
	/// It can be much quicker in places to manipulate a scored_arch_proxy
	/// than a full scored_hit_arch (because, for example, adding an index
	/// involves copying 4 bytes; adding a full hit involves copying 60 bytes
	/// plus any memory allocated for the fragments).
	///
	/// A scored_arch_proxy can be easily made back into a scored_hit_arch
	/// using make_scored_hit_arch().
	class scored_arch_proxy final {
	private:
		/// \brief The score associated with the architecture
		resscr_t the_score = INIT_SCORE;

		/// \brief The architecture, represented by the indices of its hits in the calc_hit_list
		hitidx_vec hit_indices;

	public:
		/// \brief A const_iterator type alias as part of making this a range over hit indices
		using const_iterator = hitidx_vec::const_iterator;

		scored_arch_proxy() = default;

		[[nodiscard]] const resscr_t &get_score() const;

		[[nodiscard]] bool   empty() const;
		[[nodiscard]] size_t size() const;
		const hitidx_t & operator[](const size_t &) const;

		[[nodiscard]] const_iterator begin() const;
		[[nodiscard]] const_iterator end() const;

		scored_arch_proxy & add_hit(const resscr_t &,
		                            const hitidx_t &);
	};

	/// \brief Getter for the score associated with the architecture
	inline const resscr_t & scored_arch_proxy::get_score() const {
		return the_score;
	}

	/// \brief Get whether this architecture currently contains zero entries
	inline bool scored_arch_proxy::empty() const {
		return hit_indices.empty();
	}

	/// \brief Get the number of entries in this architecture
	inline size_t scored_arch_proxy::size() const {
		return hit_indices.size();
	}

	/// \brief Get the hit index stored at the specified index
	inline const hitidx_t & scored_arch_proxy::operator[](const size_t &prm_index ///< The index of the hit index to return
	                                                      ) const {
		return hit_indices[ prm_index ];
	}

	/// \brief Standard const begin() method, as part of making this a range over hit indices
	inline auto scored_arch_proxy::begin() const -> const_iterator {
		return ::std::cbegin( hit_indices );
	}

	/// \brief Standard const end() method, as part of making this a range over hit indices
	inline auto scored_arch_proxy::end() const -> const_iterator {
		return ::std::cend( hit_indices );
	}

	/// \brief Add the specified hit index and associated score to this scored_arch_proxy
	inline scored_arch_proxy & scored_arch_proxy::add_hit(const resscr_t &prm_score,    ///< The score associated with the hit to add
	                                                      const hitidx_t &prm_hit_index ///< The index of the hit to add
	                                                      ) {
		the_score += prm_score;
		hit_indices.push_back( prm_hit_index );
		return *this;
	}

	/// \brief Add the specified hit index and associated score to a copy of the specified scored_arch_proxy
	///
	/// \relates scored_arch_proxy
	inline scored_arch_proxy add_hit_copy(scored_arch_proxy  prm_scored_arch_proxy, ///< The scored_arch_proxy to copy and then add the hit to the copy of
	                                      const resscr_t    &prm_score,             ///< The score associated with the hit to add
	                                      const hitidx_t    &prm_hit_index          ///< The index of the hit to add
	                                      ) {
		prm_scored_arch_proxy.add_hit( prm_score, prm_hit_index );
		return prm_scored_arch_proxy;
	}

	std::string to_string(const scored_arch_proxy &);
	std::ostream & operator<<(std::ostream &,
	                          const scored_arch_proxy &);

} // namespace cath::rslv

#endif // CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO_SCORED_ARCH_PROXY_HPP
