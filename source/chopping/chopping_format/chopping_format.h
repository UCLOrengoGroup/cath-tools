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

#ifndef CHOPPING_FORMAT_H_INCLUDED
#define CHOPPING_FORMAT_H_INCLUDED

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

		public:
			virtual ~chopping_format() noexcept = default;

			bool represents_fragments() const;

			domain parse_domain(const std::string &) const;

			std::unique_ptr<chopping_format> clone() const;
		};

		domain parse_domain(const chopping_format &,
		                    const std::string &,
		                    const std::string &);
	}
}

#endif
