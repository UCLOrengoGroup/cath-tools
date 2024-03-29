/// \file
/// \brief The length_getter type lists header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_LENGTH_GETTER_LENGTH_GETTER_TYPES_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_LENGTH_GETTER_LENGTH_GETTER_TYPES_HPP

#include <boost/mpl/copy.hpp>
#include <boost/mpl/vector.hpp>

// clang-format off
namespace cath::score { class geometric_mean_length_getter; }
namespace cath::score { class length_of_first_getter; }
namespace cath::score { class length_of_longer_getter; }
namespace cath::score { class length_of_second_getter; }
namespace cath::score { class length_of_shorter_getter; }
namespace cath::score { class mean_length_getter; }
namespace cath::score { class num_aligned_length_getter; }
// clang-format on

namespace cath::score {

	/// \brief TODOCUMENT
	using sym_protein_only_length_getter_types = boost::mpl::vector<
		geometric_mean_length_getter,
		length_of_longer_getter,
		length_of_shorter_getter,
		mean_length_getter
	>;

	/// \brief TODOCUMENT
	using protein_only_length_getter_types = boost::mpl::copy<
		boost::mpl::vector<
			length_of_first_getter,
			length_of_second_getter
		>,
		boost::mpl::back_inserter< sym_protein_only_length_getter_types >
	>::type;

	/// \brief TODOCUMENT
	using length_getter_types = boost::mpl::copy<
		boost::mpl::vector<
			num_aligned_length_getter
		>,
		boost::mpl::back_inserter< protein_only_length_getter_types >
	>::type;

} // namespace cath::score

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_LENGTH_GETTER_LENGTH_GETTER_TYPES_HPP
