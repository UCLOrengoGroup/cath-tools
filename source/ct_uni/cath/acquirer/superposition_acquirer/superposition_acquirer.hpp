/// \file
/// \brief The superposition_acquirer class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_SUPERPOSITION_ACQUIRER_SUPERPOSITION_ACQUIRER_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_SUPERPOSITION_ACQUIRER_SUPERPOSITION_ACQUIRER_HPP

#include <boost/tuple/tuple.hpp>

#include "cath/common/type_aliases.hpp"

#include <iostream>

// clang-format off
namespace cath::opts { class cath_superpose_options; }
namespace cath::sup { class superposition_context; }
// clang-format on

namespace cath::opts {

	/// \brief TODOCUMENT
	class superposition_acquirer {
	private:
		virtual sup::superposition_context do_get_superposition(std::ostream &) const = 0;

	public:
		superposition_acquirer() = default;
		virtual ~superposition_acquirer() noexcept = default;

		superposition_acquirer(const superposition_acquirer &) = default;
		superposition_acquirer(superposition_acquirer &&) noexcept = default;
		superposition_acquirer & operator=(const superposition_acquirer &) = default;
		superposition_acquirer & operator=(superposition_acquirer &&) noexcept = default;

		sup::superposition_context get_superposition(std::ostream &) const;

		/// \brief TODOCUMENT
		static constexpr double TOLERANCE_FOR_EQUAL_RMSDS = 0.000001;
	};

} // namespace cath::opts

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_SUPERPOSITION_ACQUIRER_SUPERPOSITION_ACQUIRER_HPP
