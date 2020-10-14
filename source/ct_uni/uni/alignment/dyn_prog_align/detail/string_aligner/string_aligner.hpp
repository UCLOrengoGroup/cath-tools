/// \file
/// \brief The string_aligner class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_STRING_ALIGNER_STRING_ALIGNER_HPP
#define _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_STRING_ALIGNER_STRING_ALIGNER_HPP

#include "common/type_aliases.hpp"

namespace cath { namespace align { class alignment; } }
namespace cath { namespace align { namespace gap { class gap_penalty; } } }

namespace cath {
	namespace align {

		namespace detail {

			/// \brief ABC-interface to define an interface for classes that align strings of upper-case letters
			///
			/// These are for simple string problems for which the identity substitution matrix is used.
			///
			/// These take pairs of strings like:
			///
			///     BFEDCBECB
			///     AAFCEDFECF
			///
			/// ...and a gap penalty and return a pair of aligned strings like:
			///
			/// This is in cath::align::detail because it is currently only used in tests
			class string_aligner {
			private:
				virtual str_str_pair do_align(const std::string &,
				                              const std::string &,
				                              const gap::gap_penalty &) const = 0;

			public:
				string_aligner() = default;
				virtual ~string_aligner() noexcept = default;

				string_aligner(const string_aligner &) = default;
				string_aligner(string_aligner &&) noexcept = default;
				string_aligner & operator=(const string_aligner &) = default;
				string_aligner & operator=(string_aligner &&) noexcept = default;

				cath::str_str_score_tpl align(const std::string &,
				                              const std::string &,
				                              const gap::gap_penalty &) const;
			};

			void check_aligned_string_matches_original(const std::string &,
			                                           const std::string &);

			void check_aligned_string_pair_is_valid(const std::string &,
			                                        const std::string &);

			void check_aligned_string_is_valid(const std::string &);

			size_size_pair get_num_gaps_and_extensions(const std::string &);

			score_type get_score_of_aligned_sequence_strings(const std::string &,
			                                                 const std::string &,
			                                                 const gap::gap_penalty &);

//			score_alignment_pair align_sequence_strings(const dyn_prog_aligner &,
//			                                            const std::string &,
//			                                            const std::string &,
//			                                            const score_type &);

			str_vec format_alignment_strings(const alignment &,
			                                 const str_vec &);

			str_str_pair format_alignment_strings(const alignment &,
			                                      const std::string &,
			                                      const std::string &);

//			str_str_pair align_and_format_sequence_strings(const std::string &,
//			                                               const std::string &,
//			                                               const score_type &,
//			                                               const size_t &);

//			str_str_pair align_and_format_sequence_strings(const std::string &,
//			                                               const std::string &,
//			                                               const score_type &);
		} // namespace detail
	} // namespace align
} // namespace cath

#endif
