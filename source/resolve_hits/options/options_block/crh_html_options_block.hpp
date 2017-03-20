/// \file
/// \brief The crh_html_options_block class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_OPTIONS_OPTIONS_BLOCK_CRH_HTML_OPTIONS_BLOCK_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_OPTIONS_OPTIONS_BLOCK_CRH_HTML_OPTIONS_BLOCK_H

#include "options/options_block/options_block.hpp"
#include "resolve_hits/options/spec/crh_html_spec.hpp"

namespace cath {
	namespace rslv {

		/// \brief Define an options_block for options specifying how cath-resolve-hits should read the input
		class crh_html_options_block final : public opts::options_block {
		private:
			using super = opts::options_block;

			/// \brief The spec this options_block configures
			crh_html_spec the_spec;

			std::unique_ptr<opts::options_block> do_clone() const final;
			std::string do_get_block_name() const final;
			void do_add_visible_options_to_description(boost::program_options::options_description &) final;
			str_opt do_invalid_string(const boost::program_options::variables_map &) const final;

		public:
			static const std::string PO_RESTRICT_HTML_WITHIN_BODY;
			static const std::string PO_MAX_NUM_NON_SOLN_HITS;
			static const std::string PO_EXCLUDE_REJECTED_HITS;

			static const str_vec ALL_BLOCK_POS;

			const crh_html_spec & get_crh_html_spec() const;
		};

		bool has_specified_crh_html_options(const boost::program_options::variables_map &arg_vm);

	} // namespace rslv
} // namespace cath

#endif
