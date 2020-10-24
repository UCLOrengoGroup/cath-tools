/// \file
/// \brief The hierarchy_ref class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_HIERARCHY_HIERARCHY_REF_HPP
#define _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_HIERARCHY_HIERARCHY_REF_HPP

namespace cath {
	namespace clust {

		/// \brief Whether an hierarchy_value in a hierarchy refers to
		///        an entry or a cluster in the next-deepest layer
		enum class hierarchy_ref : bool {
			ENTRY,  ///< Refers to an entry (ie a leaf node)
			CLUSTER ///< Refers to a cluster in the next deepest layer
		};

	} // namespace clust
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_HIERARCHY_HIERARCHY_REF_HPP
