/// \file
/// \brief The test_or_exe_run_mode literal header

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Tools project and then tweaked, eg namespaced in cath)
/// Copyright (C) 2007, Tony Lewis
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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_TEST_OR_EXE_RUN_MODE_H
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_TEST_OR_EXE_RUN_MODE_H

namespace cath {
	namespace common {

		/// \brief Represent whether running as part of a normal executable or a test
		enum class run_mode : bool {
			TEST, ///< Running as part of a test suite
			EXE   ///< Running in normal executable
		};

		/// \brief Store a static flag that acts as a global variable ( :o ) to record whether
		///        running as part of a normal binary or a test
		///
		/// The default value is run_mode::EXE; test executables are responsible for setting to run_mode::TEST
		struct run_mode_flag final {
			static run_mode value;
		};

	} // namespace common
} // namespace cath

#endif
