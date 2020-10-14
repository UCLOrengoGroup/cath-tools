/// \file
/// \brief The cath_tools_git_version definitions

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

#include <string>

#include "cath/external_info/cath_tools_git_version.hpp"
#include "cath/external_info/cath_tools_git_version_impl.hpp"

using ::std::string;

/// The git version (`git describe --tags --long`) of
/// the source directory from which this is being built
///
/// TODO: Come C++17, consider returning string_view
string cath::cath_tools_git_version() {
	return CATH_TOOLS_GIT_VERSION;
}

/// The git version (`git log -1 --date=short --pretty=format:%cd`) of
/// the source directory from which this is being built
///
/// TODO: Come C++17, consider returning string_view
string cath::cath_tools_git_date() {
	return CATH_TOOLS_GIT_DATE;
}
