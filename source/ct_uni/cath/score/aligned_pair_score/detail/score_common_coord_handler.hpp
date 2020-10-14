/// \file
/// \brief The score_common_coord_handler class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_SCORE_ALIGNED_PAIR_SCORE_DETAIL_SCORE_COMMON_COORD_HANDLER_HPP
#define _CATH_TOOLS_SOURCE_UNI_SCORE_ALIGNED_PAIR_SCORE_DETAIL_SCORE_COMMON_COORD_HANDLER_HPP

#include "cath/alignment/common_atom_selection_policy/common_atom_select_ca_policy.hpp"
#include "cath/alignment/common_residue_selection_policy/common_residue_select_all_policy.hpp"
#include "cath/common/clone/clone_ptr.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/score/score_type_aliases.hpp"
#include "cath/structure/geometry/coord_list.hpp"
#include "cath/structure/structure_type_aliases.hpp"

#include <memory>
#include <string>

namespace cath { namespace align { class alignment; } }
namespace cath { namespace align { class common_atom_selection_policy; } }
namespace cath { namespace align { class common_residue_selection_policy; } }
namespace cath { namespace geom { class coord_list; } }
namespace cath { class protein; }

namespace cath {
	namespace score {
		namespace detail {

			class score_common_coord_handler;

			bool operator<(const score_common_coord_handler &,
			               const score_common_coord_handler &);

			/// \brief A convenience class to be used in the implementation of aligned_pair_score classes
			///        for handling the selection of common coordinates from an alignment and associated
			///        proteins.
			///
			/// \todo Does this really need to be hidden away in cath::score::detail as if embarrassing?
			///       Consider making this class more general and more generally accessible.
			///
			/// This can be parameterised by the common_residue_selection_policy to use to select the
			/// common coordinates.
			///
			/// \todo Also parameterise this by the policy determining which common atoms to select
			///       this should include options for:
			///         * all backbone atoms
			///         * nitrogen
			///         * carbon_alpha
			///         * carbon
			///         * (carbon_beta?)
			///         * (oxygen?)
			///
			/// \todo Provide methods for providing strings that can used to provide more details
			///       in short names, long names and descriptions.
			class score_common_coord_handler final : private boost::equivalent<score_common_coord_handler,
			                                                 boost::totally_ordered<score_common_coord_handler> > {
			private:
				friend bool operator<(const score_common_coord_handler &,
				                      const score_common_coord_handler &);

				friend class boost::serialization::access;

				template<class archive> void serialize(archive &ar,
				                                       const size_t /*version*/
				                                       ) {
					ar & BOOST_SERIALIZATION_NVP( comm_res_seln_pol_ptr  );
					ar & BOOST_SERIALIZATION_NVP( comm_atom_seln_pol_ptr );
				}

				/// \brief The common_residue_selection_policy that should be used to select the common
				///        coordinates from the alignment and associated proteins.
				///
				/// The default (as defined by the default ctor) is common_residue_select_all_policy.
				common::clone_ptr<const align::common_residue_selection_policy> comm_res_seln_pol_ptr { align::common_residue_select_all_policy().clone() };

				/// \brief The common_atom_selection_policy that should be used to select the common
				///        atoms for the equivalent residues.
				///
				/// The default (as defined by the default ctor) is common_atom_select_ca_policy.
				common::clone_ptr<const align::common_atom_selection_policy> comm_atom_seln_pol_ptr   { align::common_atom_select_ca_policy().clone()     };

				/// \todo Consider making these public
				const align::common_residue_selection_policy & get_comm_res_seln_pol()  const;

				/// \todo Consider making these public
				const align::common_atom_selection_policy    & get_comm_atom_seln_pol() const;

				str_str_pair get_policy_description_strings() const;
				
			public:
				score_common_coord_handler() = default;
				score_common_coord_handler(const align::common_residue_selection_policy &,
				                           const align::common_atom_selection_policy &);
				score_common_coord_handler(const score_common_coord_handler &) = default;

				std::string short_suffix_string() const;
				std::string long_suffix_string() const;
				std::string description_brackets_string() const;

				str_bool_pair_vec short_name_suffixes() const;

				geom::coord_list_coord_list_pair get_common_coords(const align::alignment &,
				                                                   const protein &,
				                                                   const protein &) const;

				std::pair<cath::geom::coord_list_vec, cath::geom::coord_list_vec> get_common_coords_by_residue(const align::alignment &,
				                                                                                               const protein &,
				                                                                                               const protein &) const;
			};

			score_common_coord_handler_vec get_all_score_common_coord_handlers();

			bool operator<(const score_common_coord_handler &,
			               const score_common_coord_handler &);

			std::ostream & operator<<(std::ostream &,
			                          const score_common_coord_handler &);

		} // namespace detail
	} // namespace score
} // namespace cath
#endif
