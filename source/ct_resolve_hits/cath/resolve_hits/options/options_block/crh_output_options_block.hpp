/// \file
/// \brief The crh_output_options_block class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_OPTIONS_BLOCK_CRH_OUTPUT_OPTIONS_BLOCK_HPP
#define _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_OPTIONS_BLOCK_CRH_OUTPUT_OPTIONS_BLOCK_HPP

#include "cath/options/options_block/options_block.hpp"
#include "cath/resolve_hits/options/options_block/crh_single_output_options_block.hpp"
#include "cath/resolve_hits/options/spec/crh_output_spec.hpp"

namespace cath {
	namespace rslv {

		/// \brief Define an options_block for options specifying how cath-resolve-hits should write the output
		class crh_output_options_block final : public opts::options_block {
		private:
			using super = opts::options_block;

			// /// \brief The spec this options_block configures
			crh_output_spec the_spec;

			/// \brief The options_block that provides for the (soon-to-be) deprecated single-output options
			///
			/// \TODO Remove this after a reasonable time of it being deprecated
			crh_single_output_options_block deprecated_single_output_ob;

			std::unique_ptr<opts::options_block> do_clone() const final;
			std::string do_get_block_name() const final;
			void do_add_visible_options_to_description(boost::program_options::options_description &,
			                                           const size_t &) final;
			void do_add_hidden_options_to_description(boost::program_options::options_description &,
			                                          const size_t &) final;
			str_opt do_invalid_string(const boost::program_options::variables_map &) const final;
			str_vec do_get_all_options_names() const final;

			str_vec get_all_non_deprecated_option_names_that_clash_with_deprecated() const;
			str_vec get_all_non_deprecated_option_names_that_do_not_clash_with_deprecated() const;
			str_vec get_all_non_deprecated_option_names() const;

		public:
			static const std::string PO_HITS_TEXT_TO_FILE;
			static const std::string PO_QUIET;
			static const std::string PO_OUTPUT_TRIMMED_HITS;
			static const std::string PO_SUMMARISE_TO_FILE;
			static const std::string PO_HTML_OUTPUT_TO_FILE;
			static const std::string PO_JSON_OUTPUT_TO_FILE;
			static const std::string PO_EXPORT_CSS_FILE;
			static const std::string PO_OUTPUT_HMMER_ALN;

			const crh_output_spec & get_crh_output_spec() const;
			const crh_single_output_options_block & get_deprecated_single_output_options_block() const;
		};

		const crh_single_output_spec & get_deprecated_single_output_spec(const crh_output_options_block &);

		bool has_any_html_output(const crh_output_options_block &);

	} // namespace rslv
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_OPTIONS_BLOCK_CRH_OUTPUT_OPTIONS_BLOCK_HPP
