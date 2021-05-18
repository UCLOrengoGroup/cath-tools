/// \file
/// \brief The best_scan_arches class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO_BEST_SCAN_ARCHES_HPP
#define _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO_BEST_SCAN_ARCHES_HPP

#include <boost/core/ignore_unused.hpp>

#include "cath/common/boost_addenda/range/back.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/resolve_hits/algo/scored_arch_proxy.hpp"
#include "cath/resolve_hits/hit_arch.hpp"
#include "cath/resolve_hits/resolve_hits_type_aliases.hpp"
#include "cath/seq/seq_arrow.hpp"

#include <vector>

namespace cath {
	namespace rslv {

		/// \brief The best architectures (and their scores) that have been seen in a straight
		///        scan over various res_arrows
		///
		/// This is the key part of a single layer of dynamic-programming scanning
		///
		/// This works by storing the unique best architectures seen (as a vector of scored_arch_proxy)
		/// along with a lookup table from boundary indices to the index of the corresponding best architecture
		///
		/// This design means that little work is required for extending an architecture's region of best-ness
		/// and that lookup is very quick
		class best_scan_arches final {
		private:
			/// \brief The best architectures seen up to each of the points
			scored_arch_proxy_vec best_arches;

			/// \brief The index (in best_arches) of the best architecture seen for all residue positions
			std::vector<seq::resarw_t> bests;

		public:
			explicit best_scan_arches(const seq::residx_t &);

			[[nodiscard]] const scored_arch_proxy &get_best_scored_arch_up_to_arrow( const seq::seq_arrow & ) const;

			[[nodiscard]] const scored_arch_proxy &get_best_scored_arch_so_far() const;

			resscr_t extend_up_to_arrow(const seq::seq_arrow &);

			void add_best_up_to_arrow(const seq::seq_arrow &,
			                          const scored_arch_proxy &);
		};

		inline const resscr_t & get_best_score_up_to_arrow(const best_scan_arches &,
		                                                   const seq::seq_arrow &);
		inline const hit_arch & get_best_arch_up_to_arrow(const best_scan_arches &,
		                                                  const seq::seq_arrow &);
		inline const resscr_t & get_best_score_so_far(const best_scan_arches &);
		inline const hit_arch & get_best_arch_so_far(const best_scan_arches &);

		/// \brief Ctor from the number of residues expected
		inline best_scan_arches::best_scan_arches(const seq::residx_t &prm_num_residues ///< The number of residues expected
		                                          ) {
			bests.reserve( prm_num_residues + 1 ); // 1 more arrow than there are residues they surround

			// Initialise that the best we have seen up to now, is an empty architecture with zero score up to the start

			best_arches.emplace_back();
			bests.emplace_back( 0 );
		}

		/// \brief Get the best architecture (scored_arch_proxy) seen so far up to and including the specified arrow
		inline const scored_arch_proxy & best_scan_arches::get_best_scored_arch_up_to_arrow(const seq::seq_arrow &prm_arrow ///< The point at which we want to know the optimum solution
		                                                                                    ) const {
			const auto &index = prm_arrow.get_index();
#ifndef NDEBUG
			if ( index >= bests.size() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Index cannot go further than what has been seen"));
			}
#endif

			return best_arches[ bests[ index ] ];
		}

		/// \brief Get the best seen architecture so far
		inline const scored_arch_proxy & best_scan_arches::get_best_scored_arch_so_far() const {
			return best_arches[ bests.back() ];
		}

		/// \brief Extend the previous best seen architecture as the best up to the specified arrow
		inline resscr_t best_scan_arches::extend_up_to_arrow(const seq::seq_arrow &prm_arrow ///< The point up to which the previous best seen architecture is now known to still be the best
		                                                     ) {
			const auto &index = prm_arrow.get_index();
#ifndef NDEBUG
			if ( index + 1 < bests.size() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception(
					"Arrow index "
					+ std::to_string( index )
					+ " should be at least as high as the last seen ("
					+ std::to_string( bests.size() - 1 )
					+ ")"
				));
			}
#endif
			bests.resize( index + 1, bests.back() );

			return get_best_score_so_far( *this );
		}

		/// \brief Add the specified architecture as the best seen architecture up to the specified boundary
		///
		/// Be sure to have previously extended the previous architecture up to the boundary before
		///
		/// \pre prm_arrow must be one greater than the most recent addition / extension
		///      (else, in a debug build, an exception will be thrown)
		///
		/// \relates best_scan_arches
		inline void best_scan_arches::add_best_up_to_arrow(const seq::seq_arrow    &prm_arrow,      ///< The boundary associated with the new best
		                                                   const scored_arch_proxy &prm_scored_arch ///< The new best architecture
		                                                   ) {
#ifndef NDEBUG
			if ( prm_arrow.get_index() != bests.size() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception(
					"Arrow index "
					+ ::std::to_string( prm_arrow.get_index() )
					+ " doesn't go exactly one further than the last seen ("
					+ ::std::to_string( bests.size() - 1 )
					+ ")"
				));
			}
#else
			boost::ignore_unused( prm_arrow );
#endif
			best_arches.push_back( prm_scored_arch );
			bests.emplace_back( best_arches.size() - 1 );
		}

		/// \brief Get the best score seen so far
		///
		/// \relates best_scan_arches
		inline const resscr_t & get_best_score_so_far(const best_scan_arches &prm_best_scan_arches ///< The best_scan_arches to query
		                                              ) {
			return prm_best_scan_arches.get_best_scored_arch_so_far().get_score();
		}

	} // namespace rslv
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO_BEST_SCAN_ARCHES_HPP
