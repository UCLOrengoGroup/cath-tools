/// \file
/// \brief The crh_input_spec class header

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

#ifndef CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC_CRH_INPUT_SPEC_HPP
#define CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC_CRH_INPUT_SPEC_HPP

#include <filesystem>

#include "cath/common/path_type_aliases.hpp"
#include "cath/resolve_hits/file/hits_input_format_tag.hpp"
#include "cath/resolve_hits/resolve_hits_type_aliases.hpp"

namespace cath::rslv {

	/// \brief Specify the inputs for cath-resolve-hits
	class crh_input_spec final {
	private:
		/// \brief The input file from which data should be read
		path_opt              input_file;

		/// \brief Whether to read the input data from stdin
		bool                  read_from_stdin        = DEFAULT_READ_FROM_STDIN;

		/// \brief The format of the input data
		hits_input_format_tag input_format           = DEFAULT_INPUT_FORMAT;

		/// \brief The minimum gap length to consider when parsing an alignment
		seq::residx_t         min_gap_length         = DEFAULT_MIN_GAP_LENGTH;

		/// \brief Whether the code can assume that the input data is pre-grouped by query_id
		bool                  input_hits_are_grouped = DEFAULT_INPUT_HITS_ARE_GROUPED;

	public:
		/// \brief The default value for whether to read the input data from stdin
		static constexpr bool                  DEFAULT_READ_FROM_STDIN        = false;

		/// \brief The default value for the format of the input data
		static constexpr hits_input_format_tag DEFAULT_INPUT_FORMAT           = hits_input_format_tag::RAW_WITH_SCORES;

		/// \brief The default value for the minimum gap length to consider when parsing an alignment
		static constexpr seq::residx_t         DEFAULT_MIN_GAP_LENGTH         = 30;

		/// \brief The default value for whether the code can assume that the input data is pre-grouped by query_id
		static constexpr bool                  DEFAULT_INPUT_HITS_ARE_GROUPED = false;

		[[nodiscard]] const path_opt &             get_input_file() const;
		[[nodiscard]] const bool &                 get_read_from_stdin() const;
		[[nodiscard]] const hits_input_format_tag &get_input_format() const;
		[[nodiscard]] const seq::residx_t &        get_min_gap_length() const;
		[[nodiscard]] const bool &                 get_input_hits_are_grouped() const;

		crh_input_spec & set_input_file(const ::std::filesystem::path &);
		crh_input_spec & set_read_from_stdin(const bool &);
		crh_input_spec & set_input_format(const hits_input_format_tag &);
		crh_input_spec & set_min_gap_length(const seq::residx_t &);
		crh_input_spec & set_input_hits_are_grouped(const bool &);
	};

	str_opt get_invalid_description(const crh_input_spec &);

} // namespace cath::rslv

#endif // CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC_CRH_INPUT_SPEC_HPP
