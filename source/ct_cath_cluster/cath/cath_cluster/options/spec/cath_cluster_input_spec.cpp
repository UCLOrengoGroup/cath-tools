/// \file
/// \brief The cath_cluster_input_spec class definitions

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

#include "cath_cluster_input_spec.hpp"

#include "cath/common/exception/invalid_argument_exception.hpp"

using namespace ::cath;
using namespace ::cath::clust;
using namespace ::cath::common;
using namespace ::std::literals::string_literals;

using ::boost::none;

/// \brief Getter for an optional file from which links should be read
const path_opt & cath_cluster_input_spec::get_links_infile() const {
	return links_infile;
}

/// \brief Getter for the direction of links in the input
const link_dirn & cath_cluster_input_spec::get_link_dirn() const {
	return the_link_dirn;
}

/// \brief Getter for the index of the column from which the link values are to be read
const size_t & cath_cluster_input_spec::get_column_idx() const {
	return column_idx;
}

/// \brief Getter for an optional file from which names should be read
const path_opt & cath_cluster_input_spec::get_names_infile() const {
	return names_infile;
}

/// \brief Setter for an optional file from which links should be read
cath_cluster_input_spec & cath_cluster_input_spec::set_links_infile(const path_opt &prm_links_infile ///< An optional file from which links should be read
                                                                    ) {
	links_infile = prm_links_infile;
	return *this;
}

/// \brief Setter for the direction of links in the input
cath_cluster_input_spec & cath_cluster_input_spec::set_link_dirn(const link_dirn &prm_the_link_dirn ///< The direction of links in the input
                                                                 ) {
	the_link_dirn = prm_the_link_dirn;
	return *this;
}

/// \brief Setter for the index of the column from which the link values are to be read
cath_cluster_input_spec & cath_cluster_input_spec::set_column_idx(const size_t &prm_column_idx ///< The index of the column from which the link values are to be read
                                                                  ) {
	if ( prm_column_idx < 2 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot specify a clustering column index of 0 or 1"));
	}
	column_idx = prm_column_idx;
	return *this;
}

/// \brief Setter for an optional file from which names should be read
cath_cluster_input_spec & cath_cluster_input_spec::set_names_infile(const path_opt &prm_names_infile ///< An optional file from which names should be read
                                                                    ) {
	names_infile = prm_names_infile;
	return *this;
}

/// \brief Generate a description of any problem that makes the specified cath_cluster_input_spec invalid
///        or none otherwise
///
/// \relates cath_cluster_input_spec
str_opt cath::clust::get_invalid_description(const cath_cluster_input_spec &prm_input_spec ///< The cath_cluster_input_spec to query
                                             ) {
	if ( ! prm_input_spec.get_links_infile() ) {
		return "An input file of links must be specified"s;
	}

	return none;
}
