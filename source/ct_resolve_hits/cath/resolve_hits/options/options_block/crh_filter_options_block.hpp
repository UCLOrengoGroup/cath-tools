/// \file
/// \brief The crh_filter_options_block class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_OPTIONS_BLOCK_CRH_FILTER_OPTIONS_BLOCK_HPP
#define _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_OPTIONS_BLOCK_CRH_FILTER_OPTIONS_BLOCK_HPP

#include <string_view>

#include "cath/options/options_block/options_block.hpp"
#include "cath/resolve_hits/options/spec/crh_filter_spec.hpp"

namespace cath {
	namespace rslv {

		/// \brief Define an options_block for options specifying how cath-resolve-hits should filter input hits
		class crh_filter_options_block final : public opts::options_block {
		private:
			using super = opts::options_block;

			/// \brief The spec this options_block configures
			crh_filter_spec the_spec;

			[[nodiscard]] std::unique_ptr<opts::options_block> do_clone() const final;
			[[nodiscard]] std::string                          do_get_block_name() const final;
			void do_add_visible_options_to_description(boost::program_options::options_description &,
			                                           const size_t &) final;
			void do_add_hidden_options_to_description(boost::program_options::options_description &,
			                                          const size_t &) final;
			[[nodiscard]] str_opt do_invalid_string( const boost::program_options::variables_map & ) const final;
			[[nodiscard]] str_view_vec do_get_all_options_names() const final;

		  public:
			[[nodiscard]] const crh_filter_spec &get_crh_filter_spec() const;

			/// \brief The option name for the worst permissible evalue before a hit is ignored
			static constexpr ::std::string_view PO_WORST_PERMISSIBLE_EVALUE{ "worst-permissible-evalue" };

			/// \brief The option name for the worst permissible bitscore before a hit is ignored
			static constexpr ::std::string_view PO_WORST_PERMISSIBLE_BITSCORE{ "worst-permissible-bitscore" };

			/// \brief The option name for the worst permissible cath-resolve-hits score before a hit is ignored
			static constexpr ::std::string_view PO_WORST_PERMISSIBLE_SCORE{ "worst-permissible-score" };

			/// \brief The option name for the query IDs on which to filter the input, if any are present
			static constexpr ::std::string_view PO_FILTER_QUERY_ID{ "filter-query-id" };

			/// \brief The option name for the maximum number of query IDs to process
			static constexpr ::std::string_view PO_LIMIT_QUERIES{ "limit-queries" };

			/// \brief The option name for the (optional) minimum coverage fraction of an HMM for a hit to be considered
			static constexpr ::std::string_view PO_MIN_HMM_COVERAGE{ "min-hmm-coverage" };

			/// \brief The option name for the (optional) minimum coverage fraction of an HMM for a discontinuous hit (/^\dc_{32}$/) to be considered
			static constexpr ::std::string_view PO_MIN_DC_HMM_COVERAGE{ "min-dc-hmm-coverage" };
		};

	} // namespace rslv
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_OPTIONS_BLOCK_CRH_FILTER_OPTIONS_BLOCK_HPP
