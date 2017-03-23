/// \file
/// \brief The detail_help_options_block class header

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

#ifndef _CATH_TOOLS_SOURCE_OPTIONS_OPTIONS_BLOCK_DETAIL_HELP_OPTIONS_BLOCK_H
#define _CATH_TOOLS_SOURCE_OPTIONS_OPTIONS_BLOCK_DETAIL_HELP_OPTIONS_BLOCK_H

#include "common/type_aliases.hpp"
#include "options/options_block/options_block.hpp"

namespace cath {
	namespace opts {

		/// \brief A reusable options_block for defining blocks of detailed help options
		///
		/// The client specifies the list of help options.
		///
		/// This is done by passing the ctor a map<string, pair<string, string> > in which
		///  * each key is a strings containing a help option name
		///  * each value's first entry is a string containing the option's description for the
		///    options_description output
		///  * each value's second entry is a string containing the help that should be output
		///    if that option is specified
		class detail_help_options_block final : public options_block {
		private:
			using super = options_block;

			using str_str_str_pair_map = std::map<std::string, str_str_pair>;
			using str_str_str_pair_pair = str_str_str_pair_map::value_type;

			/// \brief The map from option name (string) to pair of description (string) and help (string) that define
			///        the behaviour of this detail_help_options_block
			str_str_str_pair_map desc_and_help_of_option_name;

			using str_bool_map = std::map<std::string, bool>;
			using str_bool_pair = str_bool_map::value_type;

			/// \brief The bool values recording which help options are requested
			///
			/// Note that the ctor doesn't have to populate this map with the correct keys
			/// because that is done automatically by do_add_options_to_description() when
			/// adding the references to the options_description.
			str_bool_map values;

			std::unique_ptr<options_block> do_clone() const final;
			std::string do_get_block_name() const final;
			void do_add_visible_options_to_description(boost::program_options::options_description &) final;
			str_opt do_invalid_string(const boost::program_options::variables_map &) const final;

		public:
			explicit detail_help_options_block(str_str_str_pair_map);

			bool has_help_string() const;
			std::string help_string() const;
		};
	} // namespace opts
} // namespace cath

#endif
