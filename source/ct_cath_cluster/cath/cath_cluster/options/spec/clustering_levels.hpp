/// \file
/// \brief The clustering_levels class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_SPEC_CLUSTERING_LEVELS_HPP
#define _CATH_TOOLS_SOURCE_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_SPEC_CLUSTERING_LEVELS_HPP

#include <boost/any.hpp>

#include "cath/clustagglom/clustagglom_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"

namespace cath {
	namespace clust {

		/// \brief Wrap a strength_vec of clustering levels in a distinctive type for Boost program_ptions validation
		struct clustering_levels final {
			/// \brief The wrapped levels
			strength_vec levels;
		};

		std::istream & operator>>(std::istream &,
		                          clustering_levels &);

		void validate(boost::any &,
		              const str_vec &,
		              clustering_levels *,
		              int);

	} // namespace clust
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_SPEC_CLUSTERING_LEVELS_HPP
