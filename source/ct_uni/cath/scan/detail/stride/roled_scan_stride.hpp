/// \file
/// \brief The roled_scan_stride class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_SCAN_DETAIL_STRIDE_ROLED_SCAN_STRIDE_HPP
#define _CATH_TOOLS_SOURCE_UNI_SCAN_DETAIL_STRIDE_ROLED_SCAN_STRIDE_HPP

#include "cath/scan/detail/scan_role.hpp"
#include "cath/scan/scan_stride.hpp"
#include "cath/scan/detail/stride/rep_strider.hpp"

namespace cath {
	namespace scan {
		namespace detail {

			/// \brief TODOCUMENT
			class roled_scan_stride final : private boost::equality_comparable<roled_scan_stride> {
			private:
				/// \brief TODOCUMENT
				scan_role the_role;

				/// \brief TODOCUMENT
				scan_stride the_stride;

			public:
				constexpr roled_scan_stride(const scan_role &,
				                            scan_stride);

				constexpr const scan_role & get_scan_role() const;
				constexpr const scan_stride & get_scan_stride() const;
			};

			/// \brief Ctor from scan_role and scan_stride
			inline constexpr roled_scan_stride::roled_scan_stride(const scan_role &prm_scan_role,  ///< TODOCUMENT
			                                                      scan_stride      prm_scan_stride ///< TODOCUMENT
			                                                      ) : the_role   { prm_scan_role                },
			                                                          the_stride { std::move( prm_scan_stride ) } {
			}

			/// \brief TODOCUMENT
			inline constexpr const scan_role & roled_scan_stride::get_scan_role() const {
				return the_role;
			}

			/// \brief TODOCUMENT
			inline constexpr const scan_stride & roled_scan_stride::get_scan_stride() const {
				return the_stride;
			}

			/// \brief TODOCUMENT
			inline constexpr bool operator==(const roled_scan_stride &prm_roled_scan_stride_a, ///< TODOCUMENT
			                                 const roled_scan_stride &prm_roled_scan_stride_b  ///< TODOCUMENT
			                                 ) {
				return (
					prm_roled_scan_stride_a.get_scan_role()   == prm_roled_scan_stride_b.get_scan_role()
					&&
					prm_roled_scan_stride_a.get_scan_stride() == prm_roled_scan_stride_b.get_scan_stride()
				);
			}

			/// \brief TODOCUMENT
			///
			/// \relates roled_scan_stride
			inline constexpr const rep_strider & get_query_from_strider(const roled_scan_stride &prm_roled_stride ///< TODOCUMENT
			                                                            ) {
				return prm_roled_stride.get_scan_stride().get_query_from_strider();
			}

			/// \brief TODOCUMENT
			///
			/// \relates roled_scan_stride
			inline constexpr const rep_strider & get_query_to_strider(const roled_scan_stride &prm_roled_stride ///< TODOCUMENT
			                                                          ) {
				return prm_roled_stride.get_scan_stride().get_query_to_strider();
			}

			/// \brief TODOCUMENT
			///
			/// \relates roled_scan_stride
			inline constexpr const rep_strider & get_index_from_strider(const roled_scan_stride &prm_roled_stride ///< TODOCUMENT
			                                                            ) {
				return prm_roled_stride.get_scan_stride().get_index_from_strider();
			}

			/// \brief TODOCUMENT
			///
			/// \relates roled_scan_stride
			inline constexpr const rep_strider & get_index_to_strider(const roled_scan_stride &prm_roled_stride ///< TODOCUMENT
			                                                          ) {
				return prm_roled_stride.get_scan_stride().get_index_to_strider();
			}



			/// \brief TODOCUMENT
			///
			/// \relates roled_scan_stride
			inline constexpr const rep_strider & get_this_from_strider(const roled_scan_stride &prm_roled_stride ///< TODOCUMENT
			                                                           ) {
				return ( prm_roled_stride.get_scan_role() == scan_role::QUERY ) ? get_query_from_strider( prm_roled_stride )
				                                                                : get_index_from_strider( prm_roled_stride );
			}

			/// \brief TODOCUMENT
			///
			/// \relates roled_scan_stride
			inline constexpr const rep_strider & get_this_to_strider(const roled_scan_stride &prm_roled_stride ///< TODOCUMENT
			                                                         ) {
				return ( prm_roled_stride.get_scan_role() == scan_role::QUERY ) ? get_query_to_strider( prm_roled_stride )
				                                                                : get_index_to_strider( prm_roled_stride );
			}

			/// \brief TODOCUMENT
			///
			/// \relates roled_scan_stride
			inline constexpr const rep_strider & get_other_from_strider(const roled_scan_stride &prm_roled_stride ///< TODOCUMENT
			                                                            ) {
				return ( prm_roled_stride.get_scan_role() == scan_role::QUERY ) ? get_index_from_strider( prm_roled_stride )
				                                                                : get_query_from_strider( prm_roled_stride );
			}

			/// \brief TODOCUMENT
			///
			/// \relates roled_scan_stride
			inline constexpr const rep_strider & get_other_to_strider(const roled_scan_stride &prm_roled_stride ///< TODOCUMENT
			                                                          ) {
				return ( prm_roled_stride.get_scan_role() == scan_role::QUERY ) ? get_index_to_strider( prm_roled_stride )
				                                                                : get_query_to_strider( prm_roled_stride );
			}



			/// \brief TODOCUMENT
			///
			/// \relates roled_scan_stride
			inline constexpr index_type get_this_from_stride(const roled_scan_stride &prm_roled_stride ///< TODOCUMENT
			                                                 ) {
				return get_this_from_strider( prm_roled_stride ).get_stride();
			}

			/// \brief TODOCUMENT
			///
			/// \relates roled_scan_stride
			inline constexpr index_type get_this_to_stride(const roled_scan_stride &prm_roled_stride ///< TODOCUMENT
			                                               ) {
				return get_this_to_strider( prm_roled_stride ).get_stride();
			}

			/// \brief TODOCUMENT
			///
			/// \relates roled_scan_stride
			inline constexpr index_type get_other_from_stride(const roled_scan_stride &prm_roled_stride ///< TODOCUMENT
			                                                  ) {
				return get_other_from_strider( prm_roled_stride ).get_stride();
			}

			/// \brief TODOCUMENT
			///
			/// \relates roled_scan_stride
			inline constexpr index_type get_other_to_stride(const roled_scan_stride &prm_roled_stride ///< TODOCUMENT
			                                                ) {
				return get_other_to_strider( prm_roled_stride ).get_stride();
			}


			/// \brief TODOCUMENT
			///
			/// \relates scan_stride
			inline constexpr index_type from_co_stride(const roled_scan_stride &prm_roled_scan_stride ///< TODOCUMENT
			                                           ) {
				return from_co_stride( prm_roled_scan_stride.get_scan_stride() );
			}

			/// \brief TODOCUMENT
			///
			/// \relates scan_stride
			inline constexpr index_type to_co_stride(const roled_scan_stride &prm_roled_scan_stride ///< TODOCUMENT
			                                         ) {
				return to_co_stride( prm_roled_scan_stride.get_scan_stride() );
			}

		} // namespace detail
	} // namespace scan
} // namespace cath

#endif
