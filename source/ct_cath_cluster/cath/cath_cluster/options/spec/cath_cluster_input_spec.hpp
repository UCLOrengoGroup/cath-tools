/// \file
/// \brief The cath_cluster_input_spec class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_SPEC_CATH_CLUSTER_INPUT_SPEC_HPP
#define _CATH_TOOLS_SOURCE_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_SPEC_CATH_CLUSTER_INPUT_SPEC_HPP

#include <cstddef>

#include "cath/clustagglom/link_dirn.hpp"
#include "cath/common/path_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"

namespace cath::clust {

	/// \brief Specify the input for cath-cluster
	class cath_cluster_input_spec final {
	private:
		/// \brief An optional file from which links should be read
		path_opt links_infile;

		/// \brief The direction of links in the input
		///        (ie whether a high number represents strong link or a high dissimilarity)
		link_dirn the_link_dirn = link_dirn::STRENGTH;

		/// \brief The index of the column from which the link values are to be read
		///        (offset 0; must be >= 2)
		size_t column_idx       = 2;

		/// \brief An optional file from which names should be read
		path_opt names_infile;

	  public:
		[[nodiscard]] const path_opt & get_links_infile() const;
		[[nodiscard]] const link_dirn &get_link_dirn() const;
		[[nodiscard]] const size_t &   get_column_idx() const;
		[[nodiscard]] const path_opt & get_names_infile() const;

		cath_cluster_input_spec &set_links_infile( const path_opt & );
		cath_cluster_input_spec &set_link_dirn( const link_dirn & );
		cath_cluster_input_spec &set_column_idx( const size_t & );
		cath_cluster_input_spec &set_names_infile( const path_opt & );
	};

	str_opt get_invalid_description(const cath_cluster_input_spec &);

} // namespace cath::clust

#endif // _CATH_TOOLS_SOURCE_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_SPEC_CATH_CLUSTER_INPUT_SPEC_HPP
