/// \file
/// \brief The load_and_scan_metrics class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_SCAN_TOOLS_LOAD_AND_SCAN_METRICS_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_SCAN_TOOLS_LOAD_AND_SCAN_METRICS_HPP

#include <filesystem>

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/sliced.hpp>
#include <boost/range/join.hpp>

#include "cath/common/algorithm/copy_build.hpp"
#include "cath/common/algorithm/generate_n_build.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/cpp17/apply.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/scan/scan_tools/scan_metrics.hpp"

namespace cath::scan {

	namespace detail {

		/// \brief TODOCUMENT
		template <typename... STRs>
		std::string join_strings(const std::string &    prm_separator, ///< TODOCUMENT
		                         const STRs        &... prm_strings    ///< TODOCUMENT
		                         ) {
			const auto strings_init_list = { prm_strings... };
			return boost::algorithm::join( strings_init_list, prm_separator );
		}

		/// \brief TODOCUMENT
		class string_joiner final {
		private:
			/// \brief TODOCUMENT
			const std::string &separator;

		public:
			/// \brief TODOCUMENT
			explicit string_joiner(const std::string &prm_separator
			                       ) : separator ( prm_separator ) {
			}

			/// \brief TODOCUMENT
			template <typename... STRs>
			std::string operator()(const STRs &...prm_strings ///< TODOCUMENT
			                       ) {
				return join_strings( separator, prm_strings... );
			}
		};

		/// \brief TODOCUMENT
		template <typename... STRs>
		std::string join_string_tuple(const std::tuple<STRs...> &prm_strings,  ///< TODOCUMENT
		                              const std::string         &prm_separator ///< TODOCUMENT
		                              ) {
			return common::apply( string_joiner( prm_separator ), prm_strings );
		}

		/// \brief TODOCUMENT
		template <typename... STRs>
		std::string markdown_table(const std::vector<std::tuple<STRs...>> &prm_data ///< TODOCUMENT
		                           ) {
			const std::string bar        = "|";
			const std::string spaced_bar = " " + bar + " ";
			const std::string rule_septr = "-------:";

			// Grab the number of elements in the tuples
			constexpr size_t num_strings = sizeof...(STRs);
			if ( prm_data.empty() ) {
				return "";
			}

			// Build a string of headers separated by |s
			const std::string headers  = spaced_bar + join_string_tuple( prm_data.front(), spaced_bar ) + spaced_bar;

			// Build a string of rules separated by |s
			const auto num_fence_posts = ( num_strings > 0 ) ? ( num_strings + 1 ) : 0;
			const auto bar_str_closure = [&] { return bar; };
			const std::string rules    = boost::algorithm::join(
				common::generate_n_build<str_vec>( num_fence_posts, bar_str_closure ),
				rule_septr
			);

			//
			const auto headers_and_rules = { headers, rules };

			// Build a
			return boost::algorithm::join(
				boost::range::join(
					headers_and_rules,
					common::transform_build<str_vec>(
						prm_data | boost::adaptors::sliced( 1, prm_data.size() ),
						[&] (const std::tuple<STRs...> &x) {
							return spaced_bar + join_string_tuple( x, spaced_bar ) + spaced_bar;
						}
					)
				),
				"\n"
			);
		}
	} // namespace detail

	/// \brief TODOCUMENT
	class load_and_scan_metrics final {
	private:
		/// \brief TODOCUMENT
		hrc_duration load_files_durn;

		/// \brief TODOCUMENT
		scan_metrics the_scan_metrics;

	public:
		load_and_scan_metrics(const hrc_duration &,
		                      scan_metrics);

		[[nodiscard]] const hrc_duration &get_load_files_durn() const;
		[[nodiscard]] const scan_metrics &get_scan_metrics() const;
	};

	const durn_mem_pair & get_query_strucs_metrics(const load_and_scan_metrics &);
	const durn_mem_pair & get_query_index_metrics(const load_and_scan_metrics &);
	const durn_mem_pair & get_index_strucs_metrics(const load_and_scan_metrics &);
	const durn_mem_pair & get_index_index_metrics(const load_and_scan_metrics &);
	const hrc_duration & get_scan_durn(const load_and_scan_metrics &);

	std::string to_markdown_string(const load_and_scan_metrics &);
	void to_markdown_file(const load_and_scan_metrics &,
	                      const ::std::filesystem::path &);

} // namespace cath::scan

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_SCAN_TOOLS_LOAD_AND_SCAN_METRICS_HPP

