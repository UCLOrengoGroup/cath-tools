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

#ifndef BEST_SCAN_ARCHES_H_INCLUDED
#define BEST_SCAN_ARCHES_H_INCLUDED

#include <boost/core/ignore_unused.hpp>

#include "common/boost_addenda/range/back.h"
#include "exception/invalid_argument_exception.h"
#include "resolve_hits/hit_arch.h"
#include "resolve_hits/res_arrow.h"
#include "resolve_hits/resolve_hits_type_aliases.h"
#include "resolve_hits/scored_arch_proxy.h"

#include <vector>

namespace cath {
	namespace rslv {

		/// \brief TODOCUMENT
		class best_scan_arches final {
		private:
			/// \brief TODOCUMENT
			std::vector<resarw_t> bests;

			/// \brief TODOCUMENT
			scored_arch_proxy_vec best_arches;

		public:
			explicit best_scan_arches(const residx_t &);

			const scored_arch_proxy & get_best_scored_arch_up_to_arrow(const res_arrow &) const;

			const scored_arch_proxy & get_best_scored_arch_so_far() const;

			resscr_t extend_up_to_arrow(const res_arrow &);

			void add_best_up_to_arrow(const res_arrow &,
			                          const scored_arch_proxy &);
		};

		inline const resscr_t & get_best_score_up_to_arrow(const best_scan_arches &,
		                                                   const res_arrow &);
		inline const hit_arch & get_best_arch_up_to_arrow(const best_scan_arches &,
		                                                  const res_arrow &);
		inline const resscr_t & get_best_score_so_far(const best_scan_arches &);
		inline const hit_arch & get_best_arch_so_far(const best_scan_arches &);

		/// \brief TODOCUMENT
		inline best_scan_arches::best_scan_arches(const residx_t &arg_num_residues ///< TODOCUMENT
		                                          ) {
			bests.reserve( arg_num_residues + 1 ); // 1 more arrow than there are residues they surround
			bests.emplace_back( 0 );
			best_arches.emplace_back();
		}

		/// \brief TODOCUMENT
		inline const scored_arch_proxy & best_scan_arches::get_best_scored_arch_up_to_arrow(const res_arrow &arg_arrow ///< TODOCUMENT
		                                                                                    ) const {
			const auto &index = arg_arrow.get_index();
#ifndef NDEBUG
			if ( index >= bests.size() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Index cannot go further than what has been seen"));
			}
#endif

			return best_arches[ bests[ index ] ];
		}

		/// \brief TODOCUMENT
		inline const scored_arch_proxy & best_scan_arches::get_best_scored_arch_so_far() const {
			return best_arches[ bests.back() ];
		}

		/// \brief TODOCUMENT
		inline resscr_t best_scan_arches::extend_up_to_arrow(const res_arrow &arg_arrow ///< TODOCUMENT
		                                                     ) {
			const auto &index = arg_arrow.get_index();
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

		/// \brief TODOCUMENT
		///
		/// \relates best_scan_arches
		inline void best_scan_arches::add_best_up_to_arrow(const res_arrow         &arg_arrow,      ///< TODOCUMENT
		                                                   const scored_arch_proxy &arg_scored_arch ///< TODOCUMENT
		                                                   ) {
#ifndef NDEBUG
			if ( arg_arrow.get_index() != bests.size() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception(
					"Arrow index "
					+ ::std::to_string( arg_arrow.get_index() )
					+ " doesn't go exactly one further than the last seen ("
					+ ::std::to_string( bests.size() - 1 )
					+ ")"
				));
			}
#else
			boost::ignore_unused( arg_arrow );
#endif
			best_arches.push_back( arg_scored_arch );
			bests.emplace_back( best_arches.size() - 1 );
		}

		// /// \brief TODOCUMENT
		// inline const resscr_t & get_best_score_up_to_arrow(const best_scan_arches &arg_best_scan_arches, ///< TODOCUMENT
		//                                                    const res_arrow        &arg_arrow             ///< TODOCUMENT
		//                                                    ) {
		// 	return arg_best_scan_arches.get_best_scored_arch_up_to_arrow( arg_arrow ).get_score();
		// }

		// /// \brief TODOCUMENT
		// ///
		// /// \relates best_scan_arches
		// inline const hit_arch & get_best_arch_up_to_arrow(const best_scan_arches &arg_best_scan_arches, ///< TODOCUMENT
		//                                                   const res_arrow        &arg_arrow             ///< TODOCUMENT
		//                                                   ) {
		// 	return arg_best_scan_arches.get_best_scored_arch_up_to_arrow( arg_arrow ).get_arch();
		// }

		/// \brief TODOCUMENT
		///
		/// \relates best_scan_arches
		inline const resscr_t & get_best_score_so_far(const best_scan_arches &arg_best_scan_arches ///< TODOCUMENT
		                                              ) {
			return arg_best_scan_arches.get_best_scored_arch_so_far().get_score();
		}

		// /// \brief TODOCUMENT
		// ///
		// /// \relates best_scan_arches
		// inline const hit_arch & get_best_arch_so_far(const best_scan_arches &arg_best_scan_arches ///< TODOCUMENT
		//                                              ) {
		// 	return arg_best_scan_arches.get_best_scored_arch_so_far().get_arch();
		// }

		// /// \brief TODOCUMENT
		// ///
		// /// \relates best_scan_arches
		// inline void add_best_up_to_arrow(best_scan_arches &arg_best_scan_arches, ///< TODOCUMENT
		//                                  const res_arrow  &arg_arrow,      ///< TODOCUMENT
		//                                  const resscr_t   &arg_best_score, ///< TODOCUMENT
		//                                  const hit_arch   &arg_arch        ///< TODOCUMENT
		//                                  ) {
		// 	return arg_best_scan_arches.add_best_up_to_arrow(
		// 		arg_arrow,
		// 		scored_arch_proxy{
		// 			arg_best_score,
		// 			arg_arch
		// 		}
		// 	);
		// }

	}
}

#endif
