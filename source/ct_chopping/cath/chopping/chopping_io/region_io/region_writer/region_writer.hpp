/// \file
/// \brief The region_writer class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_CHOPPING_IO_REGION_IO_REGION_WRITER_REGION_WRITER_HPP
#define _CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_CHOPPING_IO_REGION_IO_REGION_WRITER_REGION_WRITER_HPP

#include <string>

// clang-format off
namespace cath::chop { class region; }
// clang-format on

namespace cath::chop {

	/// \brief TODOCUMENT
	class region_writer {
	  public:
		[[nodiscard]] virtual std::string do_write_region( const region & ) const = 0;

	  protected:
		region_writer() = default;
		virtual ~region_writer() noexcept = default;

		region_writer(const region_writer &) = default;
		region_writer(region_writer &&) noexcept = default;
		region_writer & operator=(const region_writer &) = default;
		region_writer & operator=(region_writer &&) noexcept = default;

	  public:
		[[nodiscard]] std::string write_region( const region & ) const;
	};

} // namespace cath::chop

#endif // _CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_CHOPPING_IO_REGION_IO_REGION_WRITER_REGION_WRITER_HPP
