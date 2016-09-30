/// \file
/// \brief The hits_input_format_tag header

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

#ifndef HITS_INPUT_FORMAT_TAG_H_INCLUDED
#define HITS_INPUT_FORMAT_TAG_H_INCLUDED

#include <boost/any.hpp>

#include "common/algorithm/constexpr_is_uniq.h"
#include "common/type_aliases.h"

#include <array>
#include <map>

namespace cath {
	namespace rslv {

		/// \brief Represent the different formats in which resolve-hits input data can be read
		enum class hits_input_format_tag {
			HMMER_DOMTMBLOUT, ///< HMMER domtblout fomat (must assume all hits are continuous)
			HMMSEARCH_OUT,    ///< HMMer hmmsearch output fomat (can be used to deduce discontinuous hits)
			RAW_WITH_SCORES,  ///< "raw" format with scores
			RAW_WITH_EVALUES  ///< "raw" format with evalues
		};

		/// \brief A constexpr list of all hits_input_format_tags
		static constexpr std::array<hits_input_format_tag, 4> all_hits_input_format_tags { {
			hits_input_format_tag::HMMER_DOMTMBLOUT,
			hits_input_format_tag::HMMSEARCH_OUT,
			hits_input_format_tag::RAW_WITH_SCORES,
			hits_input_format_tag::RAW_WITH_EVALUES
		} };

		// Compile-time check that there aren't any duplicates in all_hits_input_format_tags
		static_assert( common::constexpr_is_uniq( all_hits_input_format_tags ), "all_hits_input_format_tags shouldn't contain repeated values" );

		/// \brief Store a constexpr record of the number of hits_input_format_tags
		static constexpr size_t num_hits_input_format_tags = std::tuple_size< decltype( all_hits_input_format_tags ) >::value;
		// static constexpr size_t num_hits_input_format_tags = common::tuple_size_v< decltype( all_hits_input_format_tags ) >;

		namespace detail {

			/// \brief Class with static getter for a map from name to hits_input_format_tag
			struct hits_input_format_tag_by_name final {
				static std::map<std::string, hits_input_format_tag> get();
			};

		} // namespace detail

		/// \brief Class with static getter for a list of all the hits_input_format_tag names
		struct all_hits_input_format_tag_names final {
			static str_vec get();
		};

		std::string to_string(const hits_input_format_tag &);

		std::ostream & operator<<(std::ostream &,
		                          const hits_input_format_tag &);

		std::istream & operator>>(std::istream &,
		                          hits_input_format_tag &);

		std::string description_of_input_format(const hits_input_format_tag &);

		void validate(boost::any &,
		              const str_vec &,
		              hits_input_format_tag *,
		              int);

	}
}

#endif
