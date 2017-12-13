/// \file
/// \brief The resolve_hits_fixture header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_TEST_RESOLVE_HITS_FIXTURE_HPP
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_TEST_RESOLVE_HITS_FIXTURE_HPP

#include <string>

namespace cath {
	namespace rslv {

		/// \brief Store data that can be reused for multiple resolve_hits tests
		class resolve_hits_fixture {
		protected:
			static const std::string example_input_raw;
			static const std::string example_output;
		};

	} // namespace rslv
} // namespace cath
#endif
