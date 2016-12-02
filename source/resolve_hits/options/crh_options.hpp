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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_OPTIONS_CRH_OPTIONS_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_OPTIONS_CRH_OPTIONS_H

#include "options/executable/executable_options.hpp"
#include "options/options_block/detail_help_options_block.hpp"
#include "resolve_hits/options/options_block/crh_filter_options_block.hpp"
#include "resolve_hits/options/options_block/crh_input_options_block.hpp"
#include "resolve_hits/options/options_block/crh_output_options_block.hpp"
#include "resolve_hits/options/options_block/crh_score_options_block.hpp"
#include "resolve_hits/options/options_block/crh_segment_options_block.hpp"

#include <iosfwd>

namespace cath { namespace rslv {class crh_spec; } }

namespace cath {
	namespace rslv {

		/// \brief Implement the executable_options for cath-resolve-hits
		class crh_options final : public opts::executable_options {
		private:
			using super = opts::executable_options;

			static const std::string STANDARD_USAGE_ERROR_STRING;

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

			/// \brief The detailed help options_block
			opts::detail_help_options_block detail_help_ob;

			virtual std::string do_get_program_name() const override final;
			virtual boost::program_options::positional_options_description get_positional_options() override final;
			virtual str_opt do_get_error_or_help_string() const override final;

			virtual std::string do_get_help_prefix_string() const override final;
			virtual std::string do_get_help_suffix_string() const override final;
			virtual std::string do_get_overview_string() const override final;

			static str_str_str_pair_map detail_help_spec();

		public:
			crh_options();

			crh_spec get_crh_spec() const;

			const crh_input_spec & get_crh_input_spec() const;
			const crh_segment_spec & get_crh_segment_spec() const;
			const crh_score_spec & get_crh_score_spec() const;
			const crh_filter_spec & get_crh_filter_spec() const;
			const crh_output_spec & get_crh_output_spec() const;

			static const std::string PROGRAM_NAME;
		};

		std::string get_crh_raw_format_help_string();
		std::string get_crh_cath_rules_help_string();
	} // namespace rslv
} // namespace cath

#endif
