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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_SEC_STRUC_CALC_DSSP_DSSP_SS_CALC_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_SEC_STRUC_CALC_DSSP_DSSP_SS_CALC_HPP

#include <optional>
#include <vector>

#include "cath/structure/structure_type_aliases.hpp"

namespace cath { class protein; }
namespace cath { namespace file { class pdb; } }
namespace cath { namespace sec { class bifur_hbond; } }
namespace cath { namespace sec { class bifur_hbond_list; } }
namespace cath { namespace sec { struct hbond_half; } }
namespace cath { namespace sec { using hbond_half_opt = ::std::optional<hbond_half>; } }
namespace cath { namespace sec { using hbond_half_opt_pair = std::pair<hbond_half_opt, hbond_half_opt>; } }

namespace cath {
	namespace sec {
		namespace detail {

			/// \brief A struct for containing some useful constants
			struct sec_struc_consts final {

				/// \brief The minimum allowable difference between residues for them to be ends of a beta bridge
				static constexpr size_t MIN_ALLOWABLE_RES_DIFF_FOR_BETA_BRIDGE = 3;

				/// \brief The maximum allowed difference between source residues for two beta bridges to
				///        be considered to enclose a beta-bulge
				static constexpr size_t BETA_BULGE_MAX_DIFF_SOURCE             = 2;

				/// \brief The maximum allowed difference between destination residues for two beta bridges to
				///        be considered to enclose a beta-bulge
				static constexpr size_t BETA_BULGE_MAX_DIFF_DEST               = 5;

				/// \brief The default type of helix
				static constexpr size_t DEFAULT_HELIX_N                        = 4;

			};

			/// \brief Represent the type of a beta-bridge
			enum class beta_bridge_type : bool {
				PARALLEL,     ///< Parallel beta-bridge
				ANTI_PARALLEL ///< Anti-parallel beta-bridge
			};

			std::string to_string(const beta_bridge_type &);
			std::ostream & operator<<(std::ostream &,
			                          const beta_bridge_type &);

			/// \brief Represent the context of a beta-bridge
			enum class beta_bridge_context : bool {
				LONE_BRIDGE, ///< The beta-bridge isn't in a beta-sheet (or that hasn't been calculated yet)
				IN_SHEET     ///< The bridge appears in a beta-sheet
			};

			std::string to_string(const beta_bridge_context &);
			std::ostream & operator<<(std::ostream &,
			                          const beta_bridge_context &);

			/// \brief Hold the information associated with a beta-bridge from a particular source residue
			struct beta_bridge {
				/// \brief The index of the residue on the other side of this beta-bridge
				size_t partner_idx;

				/// \brief The type of this beta-bridge (parallel or anti-parallel)
				beta_bridge_type type;

				/// \brief The context of this beta-bridge (alone or in a sheet)
				///
				/// Note that this is initially set to LONE_BRIDGE and the is updated
				/// to IN_SHEET depending on calculations
				beta_bridge_context context = DEFAULT_CONTEXT;

				static constexpr beta_bridge_context DEFAULT_CONTEXT = beta_bridge_context::LONE_BRIDGE;

				beta_bridge(const size_t              &prm_partner_idx,
				            const beta_bridge_type    &prm_type,
				            const beta_bridge_context &prm_context = DEFAULT_CONTEXT
				            ) : partner_idx { prm_partner_idx },
				                type        { prm_type        },
				                context     { prm_context     } {
				}

				beta_bridge(const beta_bridge &) noexcept = default;
				beta_bridge(beta_bridge &&) noexcept = default;
				beta_bridge & operator=(const beta_bridge &) noexcept = default;
				beta_bridge & operator=(beta_bridge &&) noexcept = default;
			};

			bool operator==(const beta_bridge &,
			                const beta_bridge &);

			std::string to_string(const beta_bridge &);
			std::ostream & operator<<(std::ostream &,
			                          const beta_bridge &);

			/// \brief Type alias for an optional beta_bridge
			using beta_bridge_opt     = ::std::optional<beta_bridge>;

			/// \brief Type alias for a vector of beta_bridge values
			using beta_bridge_vec     = std::vector<beta_bridge>;

			/// \brief Type alias for a vector of beta_bridge_vecs
			using beta_bridge_vec_vec = std::vector<beta_bridge_vec>;

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

			bool are_strictly_nh_to_co_bonded(const bifur_hbond_list &,
			                                  const size_t &,
			                                  const size_t &);

			bool are_nh_to_co_bonded(const bifur_hbond_list &,
			                         const size_t &,
			                         const size_t &);

			bool are_co_to_nh_bonded(const bifur_hbond_list &,
			                         const size_t &,
			                         const size_t &);

			::std::optional<helix_category> n_helix_cat(const bifur_hbond_list &,
			                                            const size_vec &,
			                                            const size_t &,
			                                            const size_t & = sec_struc_consts::DEFAULT_HELIX_N);

			bool is_n_helix_bonded_to_later(const bifur_hbond_list &,
			                                const size_vec &,
			                                const size_t &,
			                                const size_t &);

			bool could_start_n_helix(const bifur_hbond_list &,
			                         const size_vec &,
			                         const size_t &,
			                         const size_t & = sec_struc_consts::DEFAULT_HELIX_N);

			bool starts_3_helix(const bifur_hbond_list &,
			                    const size_vec &,
			                    const size_t &);

			bool starts_5_helix(const bifur_hbond_list &,
			                    const size_vec &,
			                    const size_t &);

			bool in_5_helix(const bifur_hbond_list &,
			                const size_vec &,
			                const size_t &);


			bool is_in_4_helix_not_conflicting_with_5_helix(const bifur_hbond_list &,
			                                                const size_vec &,
			                                                const size_t &);

			bool is_in_n_helix(const bifur_hbond_list &,
			                   const size_vec &,
			                   const size_t &,
			                   const size_t & = sec_struc_consts::DEFAULT_HELIX_N);

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

			bool is_beta_bulge(const beta_bridge &,
			                   const size_t &,
			                   const beta_bridge &,
			                   const size_t &);

			void add_beta_bulge(sec_struc_type_vec &,
			                    const size_t &,
			                    const size_t &);


			bool are_consecutive_bridges_in_sheet(const beta_bridge &,
			                                      const beta_bridge &);

			void set_bridges_contexts(beta_bridge_vec_vec &);

			beta_bridge_vec_vec set_bridges_contexts_copy(beta_bridge_vec_vec);

			void remove_bridges_to_chain_break_residues(beta_bridge_vec_vec &,
			                                            const size_vec &);

			beta_bridge_vec_vec remove_bridges_to_chain_break_residues_copy(beta_bridge_vec_vec,
			                                                                const size_vec &);

			bool indices_straddle_break(const size_vec &,
			                            const size_t &,
			                            const size_t &);

		} // namespace detail

		sec_struc_type_vec calc_sec_strucs(const bifur_hbond_list &,
		                                   const size_vec &);

		sec_struc_type_vec calc_sec_strucs_of_pdb__recalc_backbone_residues(const file::pdb &,
		                                                                    const ostream_ref_opt & = ::std::nullopt);

		sec_struc_type_vec calc_sec_strucs_of_backbone_complete_pdb(const file::pdb &);

		sec_struc_type_vec get_sec_strucs(const protein &);

	} // namespace sec
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_SEC_STRUC_CALC_DSSP_DSSP_SS_CALC_HPP
