/// \file
/// \brief The options_block_tester class header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef OPTIONS_BLOCK_TESTER_H_INCLUDED
#define OPTIONS_BLOCK_TESTER_H_INCLUDED

#include "common/type_aliases.h"

#include <map>
#include <string>
#include <vector>

namespace cath { namespace opts { class options_block; } }

namespace cath {
	namespace opts {

		/// \brief Provide method of parsing vector of option strings into an options_block for testing purposes
		class options_block_tester {
		private:
			static str_vec prepend_dummy_program_name_copy(const str_vec &);

		protected:
			~options_block_tester() noexcept = default;

		public:
			static void parse_into_options_block(cath::opts::options_block &,
			                                     const str_vec &);

			static const std::string          UNKNOWN_OPT;
			static const std::string          TEST_OPTION_1;
			static const std::string          TEST_OPTION_2;
			static const std::string          TEST_HELP_1;
			static const std::string          TEST_HELP_2;
			static const str_str_str_pair_map TEST_DESC_AND_HELP_OF_OPTION_NAME;
		};
	}
}

#endif