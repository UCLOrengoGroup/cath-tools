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

#include <boost/filesystem/path.hpp>

#include "common/cpp14/cbegin_cend.hpp"
#include "common/rapidjson_addenda/to_rapidjson_string.hpp"
#include "common/type_aliases.hpp"
#include "resolve_hits/full_hit.hpp"
#include "resolve_hits/options/spec/hit_boundary_output.hpp"
#include "resolve_hits/resolve_hits_type_aliases.hpp"

namespace cath { namespace rslv { class read_and_process_mgr; } }

namespace cath {
	namespace rslv {

		using doub_opt = boost::optional<double>;

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

		void read_full_hit_list_from_file(read_and_process_mgr &,
		                                  const boost::filesystem::path &);
		void read_full_hit_list_from_istream(read_and_process_mgr &,
		                                     std::istream &);
		std::string to_output_string(const full_hit_list &,
		                             const crh_segment_spec &,
		                             const hit_output_format & = hit_output_format::CLASS,
		                             const std::string & = std::string{},
		                             const hit_boundary_output & = hit_boundary_output::ORIG);
		std::string to_string(const full_hit_list &);
		std::ostream & operator<<(std::ostream &,
		                          const full_hit_list &);
		residx_opt get_max_stop(const full_hit_list &);
		resscr_opt get_best_crh_score(const full_hit_list &,
		                              const crh_score_spec &);

		seg_boundary_pair_vec resolved_boundaries(const full_hit &,
		                                          const full_hit_list &,
		                                          const crh_segment_spec &);

		hit_seg_opt_vec merge_boundaries(const hit_seg_vec &,
		                                 const seg_boundary_pair_vec &,
		                                 const crh_segment_spec &);

		hit_seg_opt_vec resolve_all_boundaries(const full_hit &,
		                                       const full_hit_list &,
		                                       const crh_segment_spec &);

		std::string get_all_resolved_segments_string(const full_hit &,
		                                             const full_hit_list &,
		                                             const crh_segment_spec &);

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

		using full_hit_tpl     = std::tuple<std::string, residx_residx_pair_vec, double>;
		using full_hit_tpl_vec = std::vector<full_hit_tpl>;

		/// \brief Write the specified full_hit_list to the specified rapidjson_writer
		template <common::json_style Style>
		void write_to_rapidjson_with_compact_fullhits(common::rapidjson_writer<Style> &arg_writer,       ///< The rapidjson_writer to which the full_hit_list should be written
		                                              const full_hit_list             &arg_full_hit_list ///< The full_hit_list to write
		                                              ) {
			arg_writer.start_array();
			for (const auto &the_full_hit : arg_full_hit_list) {
				arg_writer.write_raw_string( common::to_rapidjson_string<common::json_style::COMPACT>( the_full_hit ) );
			}
			arg_writer.end_array();
		}

		std::string to_json_string_with_compact_fullhits(const full_hit_list &);

	} // namespace rslv
} // namespace cath

#endif
