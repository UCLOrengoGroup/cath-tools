/// \file
/// \brief The hit_list class header

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

#ifndef HIT_LIST_H_INCLUDED
#define HIT_LIST_H_INCLUDED

#include <boost/log/trivial.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/range/algorithm/sort.hpp>

#include "common/chrono/duration_to_seconds_string.h"
#include "common/cpp14/cbegin_cend.h"
#include "resolve_hits/hit.h"

#include <tuple>

namespace cath { namespace rslv { class read_and_resolve_mgr; } }

namespace cath {
	namespace rslv {

		/// \brief Represent a list of hits (which can then be resolved)
		///
		/// \invariant The hits kept sorted by get_less_than_fn() (roughly, by stop, then start, then score)
		class hit_list final {
		private:
			/// \brief The list of hits
			hit_vec the_hits;

			/// \brief The list of corresponding labels
			///
			/// Note that the list of labels may not be in the same order as the
			/// list of hits; each hit has an index that indicates which is its corresponding label
			str_vec hit_labels;

			static void sort_hit_vec(hit_vec &,
			                         const str_vec &);

		public:
			/// \brief A const_iterator type alias as part of making this a range over hits
			using iterator       = hit_vec::iterator;

			/// \brief A const_iterator type alias as part of making this a range over hits
			using const_iterator = hit_vec::const_iterator;

			/// \brief Less-than fuction as used for keeping the hits sorted in the hit_list
			static auto get_less_than_fn(const str_vec &arg_hit_labels ///< The list of corresponding labels for the hits
			                             ) {
				return [&] (const hit &x, const hit &y) {
					const auto lhs_tie = std::tie( x.get_stop_arrow(), x.get_start_arrow(), x.get_score() );
					const auto rhs_tie = std::tie( y.get_stop_arrow(), y.get_start_arrow(), y.get_score() );
					if ( lhs_tie < rhs_tie) {
						return true;
					}
					else if ( rhs_tie < lhs_tie ) {
						return false;
					}

					const auto lhs_segs_str = get_segments_string( x );
					const auto rhs_segs_str = get_segments_string( y );
					return (
						std::tie( lhs_segs_str, x.get_label( arg_hit_labels ) )
						<
						std::tie( rhs_segs_str, y.get_label( arg_hit_labels ) )
					);
				};
			}

			explicit hit_list(const hit_vec &,
			                  const str_vec &);
			explicit hit_list(hit_vec &&,
			                  str_vec &&);

			size_t size() const;
			bool empty() const;

			const hit & operator[](const size_t &) const;

			const str_vec & get_labels() const;

			iterator begin();
			iterator end();
			const_iterator begin() const;
			const_iterator end() const;
		};

		void read_hit_list_from_file(read_and_resolve_mgr &,
		                             const boost::filesystem::path &);
		void read_hit_list_from_istream(read_and_resolve_mgr &,
		                                std::istream &);
		std::string to_string(const hit_list &);
		std::ostream & operator<<(std::ostream &,
		                          const hit_list &);
		residx_opt get_max_stop(const hit_list &);
		resscr_opt get_best_score(const hit_list &);
		hit_list::const_iterator find_first_hit_stopping_at_or_after(const hit_list  &,
		                                                             const res_arrow &);
		hit_list::const_iterator find_first_hit_stopping_after(const hit_list  &,
		                                                       const res_arrow &);
		boost::integer_range<hitidx_t> indices_of_hits_that_stop_in_range(const hit_list &,
		                                                                  const res_arrow &,
		                                                                  const res_arrow &);

		/// \brief Private-static method for in-place sorting hits using get_less_than_fn()
		inline void hit_list::sort_hit_vec(hit_vec       &arg_hit_vec,   ///< The hits to in-place sort
		                                   const str_vec &arg_hit_labels ///< The labels corresponding to the hits (but not necessarily in the same order)
		                                   ) {
			boost::range::sort(
				arg_hit_vec,
				get_less_than_fn( arg_hit_labels )
			);
		}

		/// \brief Ctor from lvalues
		inline hit_list::hit_list(const hit_vec &arg_hit_list,  ///< The hits
		                          const str_vec &arg_hit_labels ///< The labels corresponding to the hits (but not necessarily in the same order)
		                          ) : the_hits   ( arg_hit_list   ),
		                              hit_labels ( arg_hit_labels ) {
			sort_hit_vec( the_hits, hit_labels );
		}

		/// \brief Ctor from rvalues
		inline hit_list::hit_list(hit_vec &&arg_hit_list,  ///< The hits
		                          str_vec &&arg_hit_labels ///< The labels corresponding to the hits (but not necessarily in the same order)
		                          ) : the_hits   ( std::move( arg_hit_list   ) ),
		                              hit_labels ( std::move( arg_hit_labels ) ) {
			sort_hit_vec( the_hits, hit_labels );
		}

		/// \brief Return the number of hits
		inline size_t hit_list::size() const {
			return the_hits.size();
		}

		/// \brief Return whether there are zero hits
		inline bool hit_list::empty() const {
			return the_hits.empty();
		}

		/// \brief Return the hit stored at the specified index
		inline const hit & hit_list::operator[](const size_t &arg_index ///< The index of the hit to return
		                                        ) const {
			return the_hits[ arg_index ];
		}

		/// \brief Get the list of labels corresponding to the hits (but not necessarily in the same order)
		inline const str_vec & hit_list::get_labels() const {
			return hit_labels;
		}

		/// \brief Standard non-const begin() method, as part of making this into a range over the hits
		inline auto hit_list::begin() -> iterator {
			return std::begin( the_hits );
		}

		/// \brief Standard non-const end() method, as part of making this into a range over the hits
		inline auto hit_list::end() -> iterator {
			return std::end( the_hits );
		}

		/// \brief Standard const begin() method, as part of making this into a range over the hits
		inline auto hit_list::begin() const -> const_iterator {
			return common::cbegin( the_hits );
		}

		/// \brief Standard const end() method, as part of making this into a range over the hits
		inline auto hit_list::end() const -> const_iterator {
			return common::cend( the_hits );
		}

	}
}

#endif
