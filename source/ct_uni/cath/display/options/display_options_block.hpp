/// \file
/// \brief The display_options_block class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_DISPLAY_OPTIONS_DISPLAY_OPTIONS_BLOCK_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_DISPLAY_OPTIONS_DISPLAY_OPTIONS_BLOCK_HPP

#include "cath/display/options/display_spec.hpp"
#include "cath/options/options_block/options_block.hpp"

namespace cath { namespace opts { class pdbs_acquirer; } }
namespace cath { class display_colour_list; }

namespace cath {
	namespace opts {

		/// \brief A block of program options that are common to all viewers (eg PyMOL, Jmol etc)
		///
		/// Unlike some other options_blocks, this one has been split so that the responsibility
		/// for holding/validating the options is passed to display_spec.
		class display_options_block final : public cath::opts::options_block {
		private:
			display_spec the_display_spec;

			[[nodiscard]] std::unique_ptr<options_block> do_clone() const final;
			[[nodiscard]] std::string                    do_get_block_name() const final;
			void do_add_visible_options_to_description(boost::program_options::options_description &,
			                                           const size_t &) final;
			[[nodiscard]] str_opt do_invalid_string( const boost::program_options::variables_map & ) const final;
			[[nodiscard]] str_view_vec do_get_all_options_names() const final;

		  public:
			[[nodiscard]] display_spec get_display_spec() const;
		};

	} // namespace opts
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_DISPLAY_OPTIONS_DISPLAY_OPTIONS_BLOCK_HPP
