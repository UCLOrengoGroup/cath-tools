/// \file
/// \brief The crh_single_output_options_block class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_OPTIONS_BLOCK_CRH_SINGLE_OUTPUT_OPTIONS_BLOCK_HPP
#define _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_OPTIONS_BLOCK_CRH_SINGLE_OUTPUT_OPTIONS_BLOCK_HPP

#include <string_view>

#include "cath/options/options_block/options_block.hpp"
#include "cath/resolve_hits/options/spec/crh_single_output_spec.hpp"

namespace cath::rslv {

	/// \brief Define an options_block for options specifying how cath-resolve-hits should write the output
	class crh_single_output_options_block final : public opts::options_block {
	private:
		using super = opts::options_block;

		/// \brief The spec this options_block configures
		crh_single_output_spec the_spec;

		[[nodiscard]] std::unique_ptr<opts::options_block> do_clone() const final;
		[[nodiscard]] std::string                          do_get_block_name() const final;
		void do_add_visible_options_to_description(boost::program_options::options_description &,
		                                           const size_t &) final;
		[[nodiscard]] str_opt do_invalid_string( const boost::program_options::variables_map & ) const final;
		[[nodiscard]] str_view_vec do_get_all_options_names() const final;

	  public:
		/// \TODO Remove any options duplicated in crh_output_options_block once this
		///       is only used for deprecated options within that. Then:
		///        * propagate those changes to crh_single_output_spec
		///        * remove sort_copy_build in crh_output_options_block::do_get_all_options_names()
		[[nodiscard]] const crh_single_output_spec &get_crh_single_output_spec() const;

		/// \brief The option name for the output file to which data should be written (or unspecified for stdout)
		static constexpr ::std::string_view PO_OUTPUT_FILE{ "output-file" };

		/// \brief The option name for whether to output a summary of the input data
		static constexpr ::std::string_view PO_SUMMARISE{ "summarise" };

		/// \brief The option name for whether to output HTML describing the hits and the results
		static constexpr ::std::string_view PO_GENERATE_HTML_OUTPUT{ "html-output" };

		/// \brief The option name for whether to output the results in JSON format
		static constexpr ::std::string_view PO_JSON_OUTPUT{ "json-output" };
	};

	crh_out_format get_out_format(const crh_single_output_options_block &);

} // namespace cath::rslv

#endif // _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_OPTIONS_BLOCK_CRH_SINGLE_OUTPUT_OPTIONS_BLOCK_HPP
