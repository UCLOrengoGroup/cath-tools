/// \file
/// \brief The dssp_ss_calc class header

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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_SEC_STRUC_CALC_DSSP_SS_CALC_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_SEC_STRUC_CALC_DSSP_SS_CALC_H

#include <boost/optional/optional.hpp>

#include "structure/structure_type_aliases.hpp"

#include <vector>

namespace cath { class protein; }
namespace cath { namespace file { class pdb; } }
namespace cath { namespace sec { class bifur_hbond; } }
namespace cath { namespace sec { class bifur_hbond_list; } }
namespace cath { namespace sec { struct hbond_half; } }
namespace cath { namespace sec { using hbond_half_opt = boost::optional<hbond_half>; } }
namespace cath { namespace sec { using hbond_half_opt_pair = std::pair<hbond_half_opt, hbond_half_opt>; } }

namespace cath {
	namespace sec {
		namespace detail {

			/// \brief Represent the type of a beta-bridge
			enum class beta_bridge_type : bool {
				PARALLEL,     ///< Parallel beta-bridge
				ANTI_PARALLEL ///< Anti-parallel beta-bridge
			};

			std::string to_string(const beta_bridge_type &);
			std::ostream & operator<<(std::ostream &,
			                          const beta_bridge_type &);

			/// \brief Hold the information associated with a beta-bridge from a particular source residue
			struct beta_bridge {
				/// \brief The index of the residue on the other side of this beta-bridge
				size_t partner_idx;

				/// \brief The type of this beta-bridge (parallel or anti-parallel)
				beta_bridge_type type;
			};

			bool operator==(const beta_bridge &,
			                const beta_bridge &);

			std::string to_string(const beta_bridge &);
			std::ostream & operator<<(std::ostream &,
			                          const beta_bridge &);

			/// \brief Type alias for an optional beta_bridge
			using beta_bridge_opt = boost::optional<beta_bridge>;

			/// \brief Type alias for a vector of beta_bridge values
			using beta_bridge_vec = std::vector<beta_bridge>;

			/// \brief Represent whether a residue is 4-helix bonded to the residue 4 before this one,
			///        4 after this one, or both
			enum class helix_category : char {
				BONDED_TO_LATER_ONLY,   ///< The residue is 4-helix h-bonded to the residue 4 places after this one
				BONDED_TO_EARLIER_ONLY, ///< The residue is 4-helix h-bonded to the residue 4 places before this one
				BONDED_TO_BOTH,         ///< The residue is 4-helix h-bonded to the residue 4 places before this one and to the residue 4 places after this one
			};

			bool is_bonded_to_earlier(const helix_category &);
			bool is_bonded_to_later(const helix_category &);

			bool is_bonded_to(const hbond_half_opt_pair &,
			                  const size_t &);

			boost::optional<helix_category> four_helix_cat(const bifur_hbond &,
			                                               const size_t &);
			boost::optional<helix_category> four_helix_cat(const bifur_hbond_list &,
			                                               const size_t &);
			bool is_four_helix(const bifur_hbond_list &,
			                   const size_t &);
			bool beta_index_in_range(const bifur_hbond_list &,
			                         const size_t &);
			beta_bridge_opt has_parallel_beta_bridge_bonds_to_src(const bifur_hbond_list &,
			                                                      const size_t &,
			                                                      const size_t &);
			beta_bridge_opt has_parallel_beta_bridge_bonds_straddling_src(const bifur_hbond_list &,
			                                                              const size_t &,
			                                                              const size_t &);
			beta_bridge_opt has_antiparallel_beta_bridge_bonds_to_src(const bifur_hbond_list &,
			                                                          const size_t &,
			                                                          const size_t &);
			beta_bridge_opt has_antiparallel_beta_bridge_bonds_straddling_src(const bifur_hbond_list &,
			                                                                  const size_t &,
			                                                                  const size_t &);
			beta_bridge_vec has_parallel_beta_bridge(const bifur_hbond_list &,
			                                         const size_t &);
			beta_bridge_vec has_antiparallel_beta_bridge(const bifur_hbond_list &,
			                                             const size_t &);
			beta_bridge_vec has_beta_bridge(const bifur_hbond_list &,
			                                const size_t &);

			beta_bridge_vec has_beta_bridge(const bifur_hbond_list &,
			                                const size_t &);
			bool is_beta_bulge(const beta_bridge &,
			                   const size_t &,
			                   const beta_bridge &,
			                   const size_t &);

			void add_beta_bulge(sec_struc_type_vec &,
			                    const size_t &,
			                    const size_t &);

		} // namespace detail

		sec_struc_type_vec calc_sec_strucs(const bifur_hbond_list &);

		sec_struc_type_vec calc_sec_strucs_of_pdb__recalc_backbone_residues(const file::pdb &,
		                                                                    const ostream_ref_opt & = boost::none);

		sec_struc_type_vec calc_sec_strucs_of_backbone_complete_pdb(const file::pdb &);

		sec_struc_type_vec get_sec_strucs(const protein &);

	} // namespace sec
} // namespace cath

#endif
