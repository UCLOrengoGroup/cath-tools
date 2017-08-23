/// \file
/// \brief The clustagglom_fixture class header

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

#ifndef _CATH_TOOLS_SOURCE_CLUSTAGGLOM_CLUSTAGGLOM_FIXTURE_H
#define _CATH_TOOLS_SOURCE_CLUSTAGGLOM_CLUSTAGGLOM_FIXTURE_H

#include <boost/filesystem/path.hpp>

namespace cath {
	namespace clust {

		/// \brief A test fixture for clustagglom tests that adds some extras to global_test_constants
		class clustagglom_fixture {
		protected:
			static boost::filesystem::path CLUSTAGGLOM_DIR();
		};

	} // namespace clust
} // namespace cath

#endif
