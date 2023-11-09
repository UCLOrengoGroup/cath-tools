/// \file
/// \brief The chopping_format class header

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

#ifndef CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_CHOPPING_FORMAT_CHOPPING_FORMAT_HPP
#define CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_CHOPPING_FORMAT_CHOPPING_FORMAT_HPP

#include "cath/chopping/domain/domain.hpp"
#include "cath/chopping/region/region.hpp"
#include "cath/common/clone/check_uptr_clone_against_this.hpp"

#include <memory>
#include <string>

// clang-format off
namespace cath::chop { class domain; }
// clang-format on

namespace cath::chop {

	/// \brief TODOCUMENT
	class chopping_format {
	  private:
		[[nodiscard]] virtual std::unique_ptr<chopping_format> do_clone() const = 0;

		/// \brief TODOCUMENT
		[[nodiscard]] virtual bool do_represents_fragments() const = 0;

		/// \brief TODOCUMENT
		[[nodiscard]] virtual domain do_parse_domain( const std::string & ) const = 0;

		/// \brief Pure virtual method with which each concrete format must define how to write a region to the specified string
		[[nodiscard]] virtual std::string do_write_region( const region & ) const = 0;

		/// \brief Pure virtual method with which each concrete format must define how to write a domain to the specified string
		[[nodiscard]] virtual std::string do_write_domain( const domain & ) const = 0;

	  public:
		chopping_format() = default;
		virtual ~chopping_format() noexcept = default;

		chopping_format(const chopping_format &) = default;
		chopping_format(chopping_format &&) noexcept = default;
		chopping_format & operator=(const chopping_format &) = default;
		chopping_format & operator=(chopping_format &&) noexcept = default;

		[[nodiscard]] bool represents_fragments() const;

		[[nodiscard]] domain parse_domain( const std::string & ) const;

		[[nodiscard]] std::string write_region( const region & ) const;
		[[nodiscard]] std::string write_domain( const domain & ) const;

		[[nodiscard]] std::unique_ptr<chopping_format> clone() const;
	};

	/// \brief TODOCUMENT
	inline bool chopping_format::represents_fragments() const {
		return do_represents_fragments();
	}

	/// \brief TODOCUMENT
	inline domain chopping_format::parse_domain(const std::string &prm_domain_chopping_string ///< TODOCUMENT
	                                            ) const {
		return do_parse_domain( prm_domain_chopping_string );
	}

	/// \brief Write the specified region to a string
	inline std::string chopping_format::write_region(const region &prm_region ///< The region to write to a string
	                                                 ) const {
		return do_write_region( prm_region );
	}

	/// \brief Write the specified domain to a string
	inline std::string chopping_format::write_domain(const domain &prm_domain ///< The domain to write to a string
	                                                 ) const {
		return do_write_domain( prm_domain );
	}

	/// \brief Standard approach to achieving a virtual copy-ctor
	inline std::unique_ptr<chopping_format> chopping_format::clone() const {
	 return common::check_uptr_clone_against_this( do_clone(), *this );
	}

	domain parse_domain(const chopping_format &,
	                    const std::string &,
	                    const std::string &);

} // namespace cath::chop

#endif // CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_CHOPPING_FORMAT_CHOPPING_FORMAT_HPP
