/// \file
/// \brief The cath_assign_domains_options class header

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

#ifndef _CATH_TOOLS_SOURCE_CATH_ASSIGN_DOMAINS_OPTIONS_CATH_ASSIGN_DOMAINS_OPTIONS_HPP
#define _CATH_TOOLS_SOURCE_CATH_ASSIGN_DOMAINS_OPTIONS_CATH_ASSIGN_DOMAINS_OPTIONS_HPP

#include "cath_assign_domains/options/cath_assign_domains_options_block.hpp"
#include "options/executable/executable_options.hpp"

namespace boost { namespace program_options { class options_description; } }

namespace cath {
	namespace opts {

		/// \brief TODOCUMENT
		///
		/// \todo Add options to allow identifiable PDBs to be loaded from a directory
		/// \todo Sort out the interaction between loading PDBs and loading an alignment
		/// \todo Add options to write/read a superposition
		class cath_assign_domains_options final : public executable_options {
		private:
			/// \brief Convenience type-alias for the parent class
			using super = executable_options;

			/// \brief The cath-assign-domains options block
			cath_assign_domains_options_block the_cath_assign_domains_options_block;

			std::string do_get_program_name() const final;
			str_opt do_get_error_or_help_string() const final;

			std::string do_get_help_prefix_string() const final;
			std::string do_get_help_suffix_string() const final;
			std::string do_get_overview_string() const final;

		public:
			cath_assign_domains_options();

			const boost::filesystem::path & get_rbf_svm_file() const;
			const boost::filesystem::path & get_data_data_file() const;
			const boost::filesystem::path & get_sf_of_dom_file() const;
			const str_vec & get_forbidden_nodes() const;

			static const std::string PROGRAM_NAME;
		};

	} // namespace opts
} // namespace cath

#endif
