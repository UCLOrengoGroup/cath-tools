/// \file
/// \brief The pair_scan_action class header

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

#ifndef PAIR_SCAN_ACTION_H_INCLUDED
#define PAIR_SCAN_ACTION_H_INCLUDED

namespace cath { namespace index { namespace detail { class vcie_match_criteria; } } }
namespace cath { namespace index { class quad_find_action; } }
namespace cath { namespace index { class view_cache_index_entry; } }

namespace cath {
	namespace index {
		namespace detail {

			/// \brief TODOCUMENT
			template <typename IDX>
			class pair_scan_action final {
			private:
				/// \brief TODOCUMENT
				const IDX &index;

				/// \brief TODOCUMENT
				const vcie_match_criteria &criteria;

				/// \brief TODOCUMENT
				quad_find_action &quad_action;

			public:
				pair_scan_action(const IDX &,
				                 const vcie_match_criteria &,
				                 quad_find_action &);
				void operator()(const view_cache_index_entry &) const;
			};

			/// \brief TODOCUMENT
			template <typename IDX>
			pair_scan_action<IDX>::pair_scan_action(const IDX                 &arg_index,      ///< TODOCUMENT
			                                        const vcie_match_criteria &arg_criteria,   ///< TODOCUMENT
			                                        quad_find_action          &arg_quad_action ///< TODOCUMENT
			                                        ) : index      ( arg_index       ),
			                                            criteria   ( arg_criteria    ),
			                                            quad_action( arg_quad_action ) {
			}

			/// \brief TODOCUMENT
			template <typename IDX>
			inline void pair_scan_action<IDX>::operator()(const view_cache_index_entry &arg_entry ///< TODOCUMENT
			                                              ) const {
				index.perform_action_on_matches( arg_entry, criteria, quad_action );
			}

		}
	}
}

#endif
