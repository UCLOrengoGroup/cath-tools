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

#ifndef _CATH_TOOLS_SOURCE_CHOPPING_CHOPPING_FORMAT_CHOPPING_FORMAT_H
#define _CATH_TOOLS_SOURCE_CHOPPING_CHOPPING_FORMAT_CHOPPING_FORMAT_H

#include "chopping/domain/domain.hpp"
#include "common/clone/check_uptr_clone_against_this.hpp"

#include <memory>
#include <string>

namespace cath { namespace chop { class domain; } }

namespace cath {
	namespace chop {

		/// \brief TODOCUMENT
		class chopping_format {
		private:
			virtual std::unique_ptr<chopping_format> do_clone() const = 0;

			/// \brief TODOCUMENT
			virtual bool do_represents_fragments() const = 0;

			/// \brief TODOCUMENT
			virtual domain do_parse_domain(const std::string &) const = 0;

			/// \brief Pure virtual method with which each concrete format must define how to write a region to the specified string
			virtual std::string do_write_region(const region &) const = 0;

			/// \brief Pure virtual method with which each concrete format must define how to write a domain to the specified string
			virtual std::string do_write_domain(const domain &) const = 0;

		public:
			chopping_format() = default;
			virtual ~chopping_format() noexcept = default;

			chopping_format(const chopping_format &) = default;
			chopping_format(chopping_format &&) noexcept = default;
			chopping_format & operator=(const chopping_format &) = default;
			chopping_format & operator=(chopping_format &&) noexcept = default;

			bool represents_fragments() const;

			domain parse_domain(const std::string &) const;

			std::string write_region(const region &) const;
			std::string write_domain(const domain &) const;

			std::unique_ptr<chopping_format> clone() const;
		};

		/// \brief TODOCUMENT
		inline bool chopping_format::represents_fragments() const {
			return do_represents_fragments();
		}

		/// \brief TODOCUMENT
		inline domain chopping_format::parse_domain(const std::string &arg_domain_chopping_string ///< TODOCUMENT
		                                            ) const {
			return do_parse_domain( arg_domain_chopping_string );
		}

		/// \brief Write the specified region to a string
		inline std::string chopping_format::write_region(const region &arg_region ///< The region to write to a string
		                                                 ) const {
			return do_write_region( arg_region );
		}

		/// \brief Write the specified domain to a string
		inline std::string chopping_format::write_domain(const domain &arg_domain ///< The domain to write to a string
		                                                 ) const {
			return do_write_domain( arg_domain );
		}

		/// \brief Standard approach to achieving a virtual copy-ctor
		inline std::unique_ptr<chopping_format> chopping_format::clone() const {
		 return common::check_uptr_clone_against_this( do_clone(), *this );
		}

		domain parse_domain(const chopping_format &,
		                    const std::string &,
		                    const std::string &);

	} // namespace chop
} // namespace cath

#endif
