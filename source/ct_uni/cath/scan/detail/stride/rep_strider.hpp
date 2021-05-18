/// \file
/// \brief The rep_strider class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_STRIDE_REP_STRIDER_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_STRIDE_REP_STRIDER_HPP

#include <optional>

#include <boost/operators.hpp>
#include <boost/range/irange.hpp>

#include "cath/scan/detail/scan_type_aliases.hpp"
#include "cath/scan/detail/stride/co_stride.hpp"

namespace cath {
	namespace scan {
		namespace detail {

			/// \brief Represent the stride between residues that are selected as reps
			///
			/// A stride of 0 means that all residues are reps.
			/// A stride of 1 means that every other residue is a rep.
			///
			/// This class and related non-member, non-friend (NMNF) functions are constexpr.
			class rep_strider final : private boost::equality_comparable<rep_strider>  {
			private:
				/// \brief The stride value (number of residues in between one rep and another)
				index_type stride = 0;

			public:
				/// \brief Compiler-generated default ctor
				constexpr explicit rep_strider() = default;

				/// \brief Ctor from stride value
				constexpr explicit rep_strider(const index_type &prm_stride ///< The stride value to use
				                               ) : stride( prm_stride ) {
				}

				/// \brief Getter for the stride value
				[[nodiscard]] constexpr const index_type &get_stride() const {
					return stride;
				}
			};

			/// \brief TODOCUMENT
			constexpr bool operator==(const rep_strider &prm_rep_strider_a, ///< TODOCUMENT
			                          const rep_strider &prm_rep_strider_b  ///< TODOCUMENT
			                          ) {
				return ( prm_rep_strider_a.get_stride() == prm_rep_strider_b.get_stride() );
			}

			/// \brief Get the number of residues that the specified rep_strider
			///        allocates to the specified number of known residues.
			///
			/// \relates rep_strider
			constexpr res_rep_index_type get_num_reps_of_num_residues(const rep_strider &prm_rep_strider, ///< The rep_strider to be queried
			                                                          const index_type  &prm_num_residues ///< The number of residues for the rep_strider to handle
			                                                          ) {
				return static_cast<res_rep_index_type>(
					( prm_rep_strider.get_stride() + prm_num_residues )
					/
					( prm_rep_strider.get_stride() + 1                )
				);
			}

			/// \brief Get the index of the residue that the specified rep_strider
			///        would use as the rep with the specified index.
			///
			/// \relates rep_strider
			constexpr index_type get_index_of_rep_index(const rep_strider        &prm_rep_strider, ///< The rep_strider to be queried
			                                            const res_rep_index_type &prm_rep_index    ///< The rep index of the residue for which the absolute index is required
			                                            ) {
				return ( prm_rep_index * ( prm_rep_strider.get_stride() + 1 ) );
			}

//			constexpr res_rep_index_type get_num_reps_of_num_residues(prm_rep_strider, prm_num_residues );
//			constexpr index_type get_index_of_rep_index( prm_rep_strider, prm_rep_index );

			/// \brief TODOCUMENT
			///
			/// \relates rep_strider
			inline auto get_indices_range(const rep_strider &prm_rep_strider, ///< TODOCUMENT
			                              const index_type  &prm_num_residues ///< TODOCUMENT
			                              ) {
				const auto num_reps           = get_num_reps_of_num_residues( prm_rep_strider, prm_num_residues );
				const auto one_past_end_index = get_index_of_rep_index      ( prm_rep_strider, num_reps         );
				return boost::irange<index_type>( 0, one_past_end_index, prm_rep_strider.get_stride() + 1 );
			}

			/// \brief TODOCUMENT
			///
			/// \relates rep_strider
			inline auto get_rep_indices_range(const rep_strider &prm_rep_strider, ///< TODOCUMENT
			                                  const index_type  &prm_num_residues ///< TODOCUMENT
			                                  ) {
				const auto num_reps = get_num_reps_of_num_residues( prm_rep_strider, prm_num_residues );
				return boost::irange<res_rep_index_type>( 0, num_reps );
			}

			/// \brief TODOCUMENT
			///
			/// \relates rep_strider
			inline constexpr index_type co_stride(const rep_strider &prm_rep_strider_a, ///< TODOCUMENT
			                                      const rep_strider &prm_rep_strider_b  ///< TODOCUMENT
			                                      ) {
				return co_stride(
					prm_rep_strider_a.get_stride(),
					prm_rep_strider_b.get_stride()
				);
			}

			::std::optional<index_type> centre_index_of_index_and_next_centre_index(const index_type &,
			                                                                        const index_type &,
			                                                                        const index_type &);

			rep_rep_pair_opt get_rep_of_indices(const rep_strider &,
			                                    const index_type  &,
			                                    const rep_strider &,
			                                    const index_type  &);

//
//			/// \brief TODOCUMENT
//			///
//			/// \relates scan_stride
//			inline constexpr rep_strider get_query_from_strider(const scan_stride &prm_stride ///< The scan_stride object containing the relevant stride values
//			                                                    ) {
//				/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
//				return rep_strider{ prm_stride.get_query_from_stride() };
//			}
//
//			/// \brief TODOCUMENT
//			///
//			/// \relates scan_stride
//			inline constexpr rep_strider get_query_to_strider(const scan_stride &prm_stride ///< The scan_stride object containing the relevant stride values
//			                                                  ) {
//				/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
//				return rep_strider{ prm_stride.get_query_to_stride() };
//			}
//
//			/// \brief TODOCUMENT
//			///
//			/// \relates scan_stride
//			inline constexpr rep_strider get_index_from_strider(const scan_stride &prm_stride ///< The scan_stride object containing the relevant stride values
//			                                                    ) {
//				/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
//				return rep_strider{ prm_stride.get_index_from_stride() };
//			}
//
//			/// \brief TODOCUMENT
//			///
//			/// \relates scan_stride
//			inline constexpr rep_strider get_index_to_strider(const scan_stride &prm_stride ///< The scan_stride object containing the relevant stride values
//			                                                  ) {
//				/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
//				return rep_strider{ prm_stride.get_index_to_stride() };
//			}


		} // namespace detail
	} // namespace scan
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_STRIDE_REP_STRIDER_HPP
