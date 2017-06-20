/// \file
/// \brief The cluster_membership_file class header

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

#ifndef _CATH_TOOLS_SOURCE_CLUSTER_FILE_CLUSTER_MEMBERSHIP_FILE_H
#define _CATH_TOOLS_SOURCE_CLUSTER_FILE_CLUSTER_MEMBERSHIP_FILE_H

#include <boost/filesystem/path.hpp>

#include "cluster/new_cluster_data.hpp" // Required for the deleted function definitions
#include "cluster/old_cluster_data.hpp" // Required for the deleted function definitions

#include <iostream>
#include <string>

namespace cath { namespace common { class id_of_string; } }

namespace cath {
	namespace clust {

		old_cluster_data parse_old_membership(std::istream &,
		                                      common::id_of_string &);

		/// \brief Prevent calling with an rvalue id_of_string
		old_cluster_data parse_old_membership(std::istream &,
		                                      const common::id_of_string &&) = delete;


		old_cluster_data parse_old_membership(const std::string &,
		                                      common::id_of_string &);

		/// \brief Prevent calling with an rvalue id_of_string
		old_cluster_data parse_old_membership(const std::string &,
		                                      const common::id_of_string &&) = delete;


		old_cluster_data parse_old_membership(const boost::filesystem::path &,
		                                      common::id_of_string &);

		/// \brief Prevent calling with an rvalue id_of_string
		old_cluster_data parse_old_membership(const boost::filesystem::path &,
		                                      const common::id_of_string &&) = delete;




		new_cluster_data parse_new_membership(std::istream &,
		                                      common::id_of_string &);

		/// \brief Prevent calling with an rvalue id_of_string
		new_cluster_data parse_new_membership(std::istream &,
		                                      const common::id_of_string &&) = delete;


		new_cluster_data parse_new_membership(const std::string &,
		                                      common::id_of_string &);

		/// \brief Prevent calling with an rvalue id_of_string
		new_cluster_data parse_new_membership(const std::string &,
		                                      const common::id_of_string &&) = delete;


		new_cluster_data parse_new_membership(const boost::filesystem::path &,
		                                      common::id_of_string &);

		/// \brief Prevent calling with an rvalue id_of_string
		new_cluster_data parse_new_membership(const boost::filesystem::path &,
		                                      const common::id_of_string &&) = delete;

	} // namespace clust
} // namespace cath

#endif
