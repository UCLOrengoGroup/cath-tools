/// \file
/// \brief The full_hit_list class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_FULL_HIT_LIST_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_FULL_HIT_LIST_H

#include "common/cpp14/cbegin_cend.hpp"
#include "resolve_hits/resolve_hits_type_aliases.hpp"

namespace cath {
	namespace rslv {

		/// \brief Represent a list of full_hits (which can then be resolved)
		///
		/// This is contained within calc_hit_list
		///
		/// \invariant The full_hits kept sorted by get_less_than_fn() (roughly, by stop, then start, then score)
		class full_hit_list final {
		private:
			/// \brief The list of full_hits
			full_hit_vec the_full_hits;

		public:
			/// \brief A const_iterator type alias as part of making this a range over full_hits
			using iterator       = full_hit_vec::iterator;

			/// \brief A const_iterator type alias as part of making this a range over full_hits
			using const_iterator = full_hit_vec::const_iterator;

			/// \todo GCC 4.9.2 (on the UCL CS CentOS 6 machines) errors out if this is labelled
			///       noexcept, but it should be so change that when a compiler that old no longer
			///       needs to be supported.
			full_hit_list() = default;
			explicit full_hit_list(full_hit_vec);

			size_t size() const;
			bool empty() const;

			template <typename... Ts>
			void emplace_back(Ts &&...args);

			const full_hit & operator[](const size_t &) const;

			iterator begin();
			iterator end();
			const_iterator begin() const;
			const_iterator end() const;
		};

		std::string to_string(const full_hit_list &);
		std::ostream & operator<<(std::ostream &,
		                          const full_hit_list &);
		residx_opt get_max_stop(const full_hit_list &);


		/// \brief Ctor
		inline full_hit_list::full_hit_list(full_hit_vec arg_full_hit_list ///< The full_hits
		                                    ) : the_full_hits( std::move( arg_full_hit_list ) ) {
			// sort_full_hit_vec( the_full_hits, full_hit_labels );
		}

		/// \brief Return the number of full_hits
		inline size_t full_hit_list::size() const {
			return the_full_hits.size();
		}

		/// \brief Return whether there are zero full_hits
		inline bool full_hit_list::empty() const {
			return the_full_hits.empty();
		}

		/// \brief Emplace_back a hit (by perfect forwarding the arguments to the new full_hit's ctor)
		template <typename... Ts>
		inline void full_hit_list::emplace_back(Ts &&...args ///< The arguments to perfect-forward
		                                        ) {
			the_full_hits.emplace_back( std::forward<Ts>( args )... );
		}

		/// \brief Return the full_hit stored at the specified index
		inline const full_hit & full_hit_list::operator[](const size_t &arg_index ///< The index of the full_hit to return
		                                                  ) const {
			return the_full_hits[ arg_index ];
		}

		/// \brief Standard non-const begin() method, as part of making this into a range over the full_hits
		inline auto full_hit_list::begin() -> iterator {
			return std::begin( the_full_hits );
		}

		/// \brief Standard non-const end() method, as part of making this into a range over the full_hits
		inline auto full_hit_list::end() -> iterator {
			return std::end( the_full_hits );
		}

		/// \brief Standard const begin() method, as part of making this into a range over the full_hits
		inline auto full_hit_list::begin() const -> const_iterator {
			return common::cbegin( the_full_hits );
		}

		/// \brief Standard const end() method, as part of making this into a range over the full_hits
		inline auto full_hit_list::end() const -> const_iterator {
			return common::cend( the_full_hits );
		}

	} // namespace rslv
} // namespace cath

#endif
