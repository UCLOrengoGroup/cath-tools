/// \file
/// \brief The crh_input_options_block class header

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

#ifndef CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_OPTIONS_BLOCK_CRH_INPUT_OPTIONS_BLOCK_HPP
#define CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_OPTIONS_BLOCK_CRH_INPUT_OPTIONS_BLOCK_HPP

#include <string_view>

#include "cath/options/options_block/options_block.hpp"
#include "cath/resolve_hits/options/spec/crh_input_spec.hpp"

namespace cath::rslv {

	/// \brief Define an options_block for options specifying how cath-resolve-hits should read the input
	class crh_input_options_block final : public opts::options_block {
	private:
		using super = opts::options_block;

		/// \brief The spec this options_block configures
		crh_input_spec the_spec;

		[[nodiscard]] std::unique_ptr<opts::options_block> do_clone() const final;
		[[nodiscard]] std::string                          do_get_block_name() const final;
		void do_add_visible_options_to_description(boost::program_options::options_description &,
		                                           const size_t &) final;
		void do_add_hidden_options_to_description(boost::program_options::options_description &,
		                                          const size_t &) final;
		[[nodiscard]] str_opt do_invalid_string( const boost::program_options::variables_map & ) const final;
		[[nodiscard]] str_view_vec do_get_all_options_names() const final;

	  public:
		[[nodiscard]] const crh_input_spec &get_crh_input_spec() const;

		/// \brief The option name for the input file from which data should be read
		static constexpr ::std::string_view PO_INPUT_FILE_OR_STDIN{ "input-file-or-stdin" };

		/// \brief The option name for the format of the input data
		static constexpr ::std::string_view PO_INPUT_FORMAT{ "input-format" };

		/// \brief The option name for the minimum gap length to consider when parsing an alignment
		static constexpr ::std::string_view PO_MIN_GAP_LENGTH{ "min-gap-length" };

		/// \brief The option name for whether the code can assume that the input data is pre-grouped by query_id
		static constexpr ::std::string_view PO_INPUT_HITS_ARE_GROUPED{ "input-hits-are-grouped" };
	};

} // namespace cath::rslv

#endif // CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_OPTIONS_BLOCK_CRH_INPUT_OPTIONS_BLOCK_HPP
