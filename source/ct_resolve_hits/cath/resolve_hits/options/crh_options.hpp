/// \file
/// \brief The crh_options class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_CRH_OPTIONS_HPP
#define _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_CRH_OPTIONS_HPP

#include <iosfwd>
#include <string_view>

#include "cath/options/executable/executable_options.hpp"
#include "cath/options/options_block/detail_help_options_block.hpp"
#include "cath/resolve_hits/options/options_block/crh_filter_options_block.hpp"
#include "cath/resolve_hits/options/options_block/crh_html_options_block.hpp"
#include "cath/resolve_hits/options/options_block/crh_input_options_block.hpp"
#include "cath/resolve_hits/options/options_block/crh_output_options_block.hpp"
#include "cath/resolve_hits/options/options_block/crh_score_options_block.hpp"
#include "cath/resolve_hits/options/options_block/crh_segment_options_block.hpp"

// clang-format off
namespace cath::rslv { class crh_spec; }
// clang-format on

namespace cath::rslv {

	/// \brief Implement the executable_options for cath-resolve-hits
	class crh_options final : public opts::executable_options {
	private:
		using super = opts::executable_options;

		/// \brief The cath-resolve-hits input options_block
		crh_input_options_block         the_input_ob;

		/// \brief The cath-resolve-hits segment options_block
		crh_segment_options_block       the_segment_ob;

		/// \brief The cath-resolve-hits score options_block
		crh_score_options_block         the_score_ob;

		/// \brief The cath-resolve-hits filter options_block
		crh_filter_options_block        the_filter_ob;

		/// \brief The cath-resolve-hits output options_block
		crh_output_options_block        the_output_ob;

		/// \brief The cath-resolve-hits html options_block
		crh_html_options_block          the_html_ob;

		/// \brief The detailed help options_block
		opts::detail_help_options_block detail_help_ob;

		[[nodiscard]] std::string_view                         do_get_program_name() const final;
		boost::program_options::positional_options_description get_positional_options() final;
		[[nodiscard]] str_opt                                  do_get_error_or_help_string() const final;

		[[nodiscard]] std::string do_get_help_prefix_string() const final;
		[[nodiscard]] std::string do_get_help_suffix_string() const final;
		[[nodiscard]] std::string do_get_overview_string() const final;

		static str_str_str_pair_map detail_help_spec();

	public:
		crh_options();

		[[nodiscard]] crh_spec get_crh_spec() const;

		[[nodiscard]] const crh_input_spec &  get_crh_input_spec() const;
		[[nodiscard]] const crh_segment_spec &get_crh_segment_spec() const;
		[[nodiscard]] const crh_score_spec &  get_crh_score_spec() const;
		[[nodiscard]] const crh_filter_spec & get_crh_filter_spec() const;
		[[nodiscard]] const crh_output_spec & get_crh_output_spec() const;
		[[nodiscard]] const crh_html_spec &   get_crh_html_spec() const;

		/// The name of the program that uses this executable_options
		static constexpr ::std::string_view PROGRAM_NAME{ "cath-resolve-hits" };
	};

	std::string get_crh_raw_format_help_string();
	std::string get_crh_cath_rules_help_string();

} // namespace cath::rslv

#endif // _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_CRH_OPTIONS_HPP
