/// \file
/// \brief The string_options_block class header

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

#ifndef _CATH_TOOLS_SOURCE_OPTIONS_OPTIONS_BLOCK_STRING_OPTIONS_BLOCK_H
#define _CATH_TOOLS_SOURCE_OPTIONS_OPTIONS_BLOCK_STRING_OPTIONS_BLOCK_H

#include "options/options_block/options_block.hpp"

#include <string>

namespace cath {
	namespace opts {

		/// \brief No-options options block that can be used to bung a string into the options usage text
		///
		/// The slight annoyance is that a ':' will be automatically appended to the string
		class string_options_block final : public options_block {
		private:
			/// \brief A useful type alias for the parent class
			using super = options_block;

			/// \brief The string to insert into the options usage
			std::string the_string;

			std::unique_ptr<options_block> do_clone() const final;
			std::string do_get_block_name() const final;
			void do_add_visible_options_to_description(boost::program_options::options_description &,
			                                           const size_t &) final;
			str_opt do_invalid_string(const boost::program_options::variables_map &) const final;
			str_vec do_get_all_options_names() const final;

		public:
			explicit string_options_block(std::string);

			const std::string & get_string() const;
		};

	} // namespace opts
} // namespace cath

#endif
