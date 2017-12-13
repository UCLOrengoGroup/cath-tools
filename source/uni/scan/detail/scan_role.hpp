/// \file
/// \brief The scan_role header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_SCAN_DETAIL_SCAN_ROLE_HPP
#define _CATH_TOOLS_SOURCE_UNI_SCAN_DETAIL_SCAN_ROLE_HPP

namespace cath {
	namespace scan {
		namespace detail {

			/// \brief The role in a scan of being either a query or an index
			///
			/// Motivation:
			///  * scan_stride objects store the from/to stride for both the query and the index
			///  * some code needs to know not only its own from/to stride but also that in its counterparts
			///    (because two counterparts need to store the same pattern of neighbours around each rep)
			///  * such code can conviently store four strides in a single scan_stride object but then
			///    must also store whether their side is the query or index side
			///
			/// The scan_stride and scan_role are packaged together in the roled_scan_stride
			enum class scan_role : bool {
				QUERY, ///< The role of being a query that can be searched against indices
				INDEX  ///< The role of being an index, against which queries can be scanned
			};

		} // namespace detail
	} // namespace scan
} // namespace cath

#endif
