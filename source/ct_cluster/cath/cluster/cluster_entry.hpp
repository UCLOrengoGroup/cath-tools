/// \file
/// \brief The cluster_entry class header

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

#ifndef CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_CLUSTER_ENTRY_HPP
#define CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_CLUSTER_ENTRY_HPP

#include "cath/seq/seq_seg_run.hpp"

#include <string>

namespace cath::clust {

	/// \brief Store the data of a cluster entry
	class cluster_entry final {
	private:
		/// \brief The name of the cluster_entry
		std::string name;

		/// \brief The segments of the cluster_entry
		seq::seq_seg_run segments;

	public:
		/// \brief No default ctor for now (because seq_seg_run doesn't default construct)
		cluster_entry() = delete;

		/// \brief Ctor from name and segments
		explicit cluster_entry(std::string      prm_name, ///< The name of the cluster_entry
		                       seq::seq_seg_run prm_segs  ///< The segments of the cluster_entry
		                       ) : name    { std::move( prm_name ) },
		                           segments{ std::move( prm_segs ) } {
		}
	};

} // namespace cath::clust

#endif // CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_CLUSTER_ENTRY_HPP
