/// \file
/// \brief The clustagglom_fixture class definitions

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

#include "clustagglom_fixture.hpp"

#include <filesystem>

#include "cath/test/global_test_constants.hpp"

using namespace ::cath::clust;

using ::std::filesystem::path;

/// \brief The directory in which clustagglom data is held
path clustagglom_fixture::CLUSTAGGLOM_DIR() {
	return global_test_constants::TEST_SOURCE_DATA_DIR() / "clustagglom";
}
