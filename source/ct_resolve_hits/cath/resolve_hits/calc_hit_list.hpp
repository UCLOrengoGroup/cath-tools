/// \file
/// \brief The calc_hit_list class header

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

#ifndef CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_CALC_HIT_LIST_HPP
#define CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_CALC_HIT_LIST_HPP

#include <filesystem>
#include <tuple>

#include <boost/range/algorithm/sort.hpp>

#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/chrono/duration_to_seconds_string.hpp"
#include "cath/resolve_hits/calc_hit.hpp"
#include "cath/resolve_hits/full_hit.hpp"
#include "cath/resolve_hits/full_hit_list.hpp"
#include "cath/resolve_hits/options/spec/crh_filter_spec.hpp"
#include "cath/resolve_hits/score_functions.hpp"
#include "cath/resolve_hits/seg_dupl_hit_policy.hpp"

// clang-format off
namespace cath::rslv { class read_and_process_mgr; }
namespace cath::rslv { class crh_segment_spec; }
// clang-format on

namespace cath::rslv {

	/// \brief Represent a list of hits (which can then be resolved)
	///
	/// This contains a full full_hit_list inside
	///
	/// \invariant The hits kept sorted by get_less_than_fn() (roughly, by stop, then start, then score)
	class calc_hit_list final {
	private:
		/// \brief The list of full_hits from which these hits were drawn
		///
		/// Note that the list may not be in the same order as the
		/// list of hits; each calc_hit has an index that indicates which is its corresponding full_hit
		full_hit_list full_hits;

		/// \brief The list of hits
		calc_hit_vec the_hits;

		static void sort_hit_vec(calc_hit_vec &,
		                         const full_hit_list &);

	public:
		/// \brief A const_iterator type alias as part of making this a range over hits
		using iterator       = calc_hit_vec::iterator;

		/// \brief A const_iterator type alias as part of making this a range over hits
		using const_iterator = calc_hit_vec::const_iterator;

		/// \brief Less-than function as used for keeping the hits sorted in the calc_hit_list
		static auto get_less_than_fn(const full_hit_list &prm_full_hits ///< The full_hits from which these hits were drawn
		                             ) {
			return [&] (const calc_hit &x, const calc_hit &y) {
				const auto num_segs_lhs = get_num_segments( x );
				const auto num_segs_rhs = get_num_segments( y );
				const auto lhs_tie = std::tie( get_stop_arrow( x ), get_start_arrow( x ), x.get_score(), num_segs_lhs );
				const auto rhs_tie = std::tie( get_stop_arrow( y ), get_start_arrow( y ), y.get_score(), num_segs_rhs );
				if ( lhs_tie < rhs_tie) {
					return true;
				}
				if ( rhs_tie < lhs_tie ) {
					return false;
				}

				for (const size_t &seg_ctr : common::indices( get_num_segments( x ) ) ) {
					const auto lhs_seg_tie = std::tie( get_start_arrow_of_segment( x, seg_ctr ), get_stop_arrow_of_segment( x, seg_ctr ) );
					const auto rhs_seg_tie = std::tie( get_start_arrow_of_segment( y, seg_ctr ), get_stop_arrow_of_segment( y, seg_ctr ) );
					if ( lhs_seg_tie < rhs_seg_tie) {
						return true;
					}
					if ( rhs_seg_tie < lhs_seg_tie ) {
						return false;
					}
				}

				return (
					prm_full_hits[ x.get_label_idx() ].get_label()
					<
					prm_full_hits[ y.get_label_idx() ].get_label()
				);
			};
		}

		explicit calc_hit_list(full_hit_list,
		                       const crh_score_spec &,
		                       const crh_segment_spec &,
		                       const crh_filter_spec & = make_accept_all_filter_spec(),
		                       const seg_dupl_hit_policy & = seg_dupl_hit_policy::PRESERVE);

		[[nodiscard]] size_t size() const;
		[[nodiscard]] bool   empty() const;

		const calc_hit & operator[](const size_t &) const;

		[[nodiscard]] const full_hit_list &get_full_hits() const;

		iterator begin();
		iterator end();
		[[nodiscard]] const_iterator begin() const;
		[[nodiscard]] const_iterator end() const;
	};

	calc_hit_vec make_hit_list_from_full_hit_list(const full_hit_list &,
	                                              const crh_score_spec &,
	                                              const crh_segment_spec &,
	                                              const crh_filter_spec & = make_accept_all_filter_spec());
	calc_hit_vec make_sorted_pruned_calc_hit_vec(const full_hit_list &,
	                                             const crh_score_spec &,
	                                             const crh_segment_spec &,
	                                             const crh_filter_spec &,
	                                             const seg_dupl_hit_policy &);

	void read_hit_list_from_file(read_and_process_mgr &,
	                             const ::std::filesystem::path &,
	                             const hit_score_type &);
	void read_hit_list_from_istream(read_and_process_mgr &,
	                                std::istream &,
	                                const hit_score_type &);
	std::string to_string(const calc_hit_list &);
	std::ostream & operator<<(std::ostream &,
	                          const calc_hit_list &);
	seq::residx_opt get_max_stop(const calc_hit_list &);
	resscr_opt get_best_score(const calc_hit_list &);
	calc_hit_list::const_iterator find_first_hit_stopping_at_or_after(const calc_hit_list  &,
	                                                                  const seq::seq_arrow &);
	calc_hit_list::const_iterator find_first_hit_stopping_after(const calc_hit_list  &,
	                                                            const seq::seq_arrow &);
	boost::integer_range<hitidx_t> indices_of_hits_that_stop_in_range(const calc_hit_list &,
	                                                                  const seq::seq_arrow &,
	                                                                  const seq::seq_arrow &);

	calc_hit_vec_citr_vec identify_redundant_hits(const calc_hit_vec &,
	                                              const full_hit_list &);

	void remove_redundant_hits(calc_hit_vec &,
	                           const full_hit_list &);


	/// \brief Private-static method for in-place sorting hits using get_less_than_fn()
	inline void calc_hit_list::sort_hit_vec(calc_hit_vec        &prm_hit_vec,  ///< The hits to in-place sort
	                                        const full_hit_list &prm_full_hits ///< The full_hits from which these hits were drawn
	                                        ) {
		boost::range::sort(
			prm_hit_vec,
			get_less_than_fn( prm_full_hits )
		);
	}

	/// \brief Ctor
	inline calc_hit_list::calc_hit_list(full_hit_list              prm_full_hits,        ///< The full_hits from which these hits are to be drawn
	                                    const crh_score_spec      &prm_score_spec,       ///< The crh_score_spec to specify how the crh-scores are to be calculated from the full-hits
	                                    const crh_segment_spec    &prm_crh_segment_spec, ///< The crh_segment_spec to specify how the segments are to be handled before being put into the hits for calculation
	                                    const crh_filter_spec     &prm_filter_spec,      ///< The crh_filter_spec specifying how hits should be filtered
	                                    const seg_dupl_hit_policy &prm_policy            ///< Whether the strictly-worse hits should be pruned
	                                    ) : full_hits { std::move( prm_full_hits ) },
	                                        the_hits  { make_sorted_pruned_calc_hit_vec(
	                                        	full_hits,
	                                        	prm_score_spec,
	                                        	prm_crh_segment_spec,
	                                        	prm_filter_spec,
	                                        	prm_policy
	                                        ) } {
		remove_redundant_hits( the_hits, full_hits );
	}

	/// \brief Return the number of hits
	inline size_t calc_hit_list::size() const {
		return the_hits.size();
	}

	/// \brief Return whether there are zero hits
	inline bool calc_hit_list::empty() const {
		return the_hits.empty();
	}

	/// \brief Return the calc_hit stored at the specified index
	inline const calc_hit & calc_hit_list::operator[](const size_t &prm_index ///< The index of the calc_hit to return
	                                                  ) const {
		return the_hits[ prm_index ];
	}

	/// \brief Get the list of labels corresponding to the hits (but not necessarily in the same order)
	inline const full_hit_list & calc_hit_list::get_full_hits() const {
		return full_hits;
	}

	/// \brief Standard non-const begin() method, as part of making this into a range over the hits
	inline auto calc_hit_list::begin() -> iterator {
		return std::begin( the_hits );
	}

	/// \brief Standard non-const end() method, as part of making this into a range over the hits
	inline auto calc_hit_list::end() -> iterator {
		return std::end( the_hits );
	}

	/// \brief Standard const begin() method, as part of making this into a range over the hits
	inline auto calc_hit_list::begin() const -> const_iterator {
		return ::std::cbegin( the_hits );
	}

	/// \brief Standard const end() method, as part of making this into a range over the hits
	inline auto calc_hit_list::end() const -> const_iterator {
		return ::std::cend( the_hits );
	}

	using hit_tpl     = std::tuple<std::string, seq::residx_residx_pair_vec, double>;
	using hit_tpl_vec = std::vector<hit_tpl>;

} // namespace cath::rslv

#endif // CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_CALC_HIT_LIST_HPP
