/// \file
/// \brief The options_block_tester class header

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

#ifndef CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_OPTIONS_BLOCK_OPTIONS_BLOCK_TESTER_HPP
#define CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_OPTIONS_BLOCK_OPTIONS_BLOCK_TESTER_HPP

#include <map>
#include <string>
#include <string_view>
#include <vector>

#include "cath/common/type_aliases.hpp"

// clang-format off
namespace cath::opts { class options_block; }
// clang-format on

namespace cath::opts {

	/// \brief Provide method of parsing vector of option strings into an options_block for testing purposes
	class options_block_tester {
	private:
		static str_vec prepend_dummy_program_name_copy(const str_vec &);

	protected:
		options_block_tester() = default;
		~options_block_tester() noexcept = default;
		options_block_tester(const options_block_tester &) noexcept = default;
		options_block_tester(options_block_tester &&) noexcept = default;
		options_block_tester & operator=(const options_block_tester &) noexcept = default;
		options_block_tester & operator=(options_block_tester &&) noexcept = default;

	public:
		static void parse_into_options_block(cath::opts::options_block &,
		                                     const str_vec &);

		template <typename OB>
		static OB parse_into_options_block_copy(OB,
		                                        const str_vec &);

		static constexpr ::std::string_view UNKNOWN_OPT{ "--it-does-not-know-me" };
		static constexpr ::std::string_view TEST_OPTION_1{ "test_help_option_1" };
		static constexpr ::std::string_view TEST_OPTION_2{ "test_help_option_2" };
		static constexpr ::std::string_view TEST_HELP_1{ "This is the first piece of help" };
		static constexpr ::std::string_view TEST_HELP_2{ "This is the second piece of help" };

		static str_str_str_pair_map TEST_DESC_AND_HELP_OF_OPTION_NAME();
	};

	/// \brief For a concrete options_block, parse the specified options into a copy of the specified block
	template <typename OB>
	OB options_block_tester::parse_into_options_block_copy(OB             prm_options_block, ///< The options_block from which a copy should be taken that then has the options parsed into it
	                                                       const str_vec &prm_options        ///< A vector of options strings to parse into the options_block (without the program name at the start - a dummy program name will be added)
	                                                       ) {
		parse_into_options_block( prm_options_block, prm_options );
		return prm_options_block;
	}

} // namespace cath::opts

#endif // CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_OPTIONS_BLOCK_OPTIONS_BLOCK_TESTER_HPP
