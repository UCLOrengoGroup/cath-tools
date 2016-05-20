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

#include "common/c++14/cbegin_cend.h"
#include "common/chrono/duration_to_seconds_string.h"
#include "resolve_hits/hit.h"

#include <tuple>

namespace cath {
	namespace rslv {

		/// \brief TODOCUMENT
		class hit_list final {
		private:
			/// \brief TODOCUMENT
			hit_vec the_hits;

			static void sort_hit_vec(hit_vec &);

		public:
			using iterator       = hit_vec::iterator;
			using const_iterator = hit_vec::const_iterator;

			static auto get_less_than_fn() {
				return [] (const hit &x, const hit &y) {
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
						std::tie( lhs_segs_str, x.get_label() )
						<
						std::tie( rhs_segs_str, y.get_label() )
					);

//					// Not using tie() trick because sorting on (unsigned) stop descending and start ascending
//					const auto &x_stop = x.get_stop();
//					const auto &y_stop = y.get_stop();
//					if      ( x_stop > y_stop ) {
//						return true;
//					}
//					else if ( x_stop < y_stop ) {
//						return false;
//					}
//					else {
//						return x.get_start() < y.get_start();
//					}
				};
			}

			explicit hit_list(const hit_vec &);

			size_t size() const;
			bool empty() const;

			const hit & operator[](const size_t &) const;

			iterator begin();
			iterator end();
			const_iterator begin() const;
			const_iterator end() const;
		};

		hit_list read_hit_list_from_file(const boost::filesystem::path &);
		std::string to_string(const hit_list &);
		std::ostream & operator<<(std::ostream &,
		                          const hit_list &);
		residx_t get_max_stop(const hit_list &);
		hit_list::const_iterator find_first_hit_stopping_at_or_after(const hit_list  &,
		                                                             const res_arrow &);
		hit_list::const_iterator find_first_hit_stopping_after(const hit_list  &,
		                                                       const res_arrow &);
		boost::integer_range<hitidx_t> indices_of_hits_that_stop_in_range(const hit_list &,
		                                                                  const res_arrow &,
		                                                                  const res_arrow &);

		/// \brief TODOCUMENT
		inline void hit_list::sort_hit_vec(hit_vec &arg_hit_vec ///< TODOCUMENT
		                                   ) {
			// const auto sort_start_time = std::chrono::high_resolution_clock::now();
			boost::range::sort(
				arg_hit_vec,
				get_less_than_fn()
			);
			// BOOST_LOG_TRIVIAL( warning ) << "Sorting takes   : " << common::durn_to_seconds_string( std::chrono::high_resolution_clock::now() - sort_start_time ) << "\n";
		}

		/// \brief TODOCUMENT
		inline hit_list::hit_list(const hit_vec &arg_hit_list ///< TODOCUMENT
		                          ) : the_hits( arg_hit_list ) {
			sort_hit_vec( the_hits );
		}

		/// \brief TODOCUMENT
		inline size_t hit_list::size() const {
			return the_hits.size();
		}

		/// \brief TODOCUMENT
		inline bool hit_list::empty() const {
			return the_hits.empty();
		}

		/// \brief TODOCUMENT
		inline const hit & hit_list::operator[](const size_t &arg_index ///< TODOCUMENT
		                                        ) const {
			return the_hits[ arg_index ];
		}

		/// \brief TODOCUMENT
		inline auto hit_list::begin() -> iterator {
			return std::begin( the_hits );
		}

		/// \brief TODOCUMENT
		inline auto hit_list::end() -> iterator {
			return std::end( the_hits );
		}

		/// \brief TODOCUMENT
		inline auto hit_list::begin() const -> const_iterator {
			return common::cbegin( the_hits );
		}

		/// \brief TODOCUMENT
		inline auto hit_list::end() const -> const_iterator {
			return common::cend( the_hits );
		}

	}
}

#endif
