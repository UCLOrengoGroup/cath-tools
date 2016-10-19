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

#ifndef DISPLAY_OPTIONS_BLOCK_H_INCLUDED
#define DISPLAY_OPTIONS_BLOCK_H_INCLUDED

#include "display/options/display_spec.h"
#include "options/options_block/options_block.h"

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
			static const std::string PO_VIEWER_COLOURS;
			static const std::string PO_GRADIENT_COLOUR_ALIGNMENT;
			static const std::string PO_SHOW_SCORES_IF_PRESENT;
			static const std::string PO_SCORES_TO_EQUIVS;
			static const std::string PO_NORMALISE_SCORES;

			static const str_vec ALL_BLOCK_POS;

			display_spec the_display_spec;

			virtual std::unique_ptr<options_block> do_clone() const override final;
			virtual std::string do_get_block_name() const override final;
			virtual void do_add_visible_options_to_description(boost::program_options::options_description &) override final;
			virtual str_opt do_invalid_string(const boost::program_options::variables_map &) const override final;

		public:
			virtual ~display_options_block() noexcept = default;

			display_spec get_display_spec() const;

			bool has_specified_options(const boost::program_options::variables_map &arg_vm) const;
		};

	}
}

#endif
