/// \file
/// \brief The cath_hit_resolver class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_CATH_HIT_RESOLVER_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_CATH_HIT_RESOLVER_H

#include "common/type_aliases.hpp"

#include <iostream>

namespace cath { namespace rslv { class crh_options; } }
namespace cath { namespace rslv { class crh_spec; } }

namespace cath {
	namespace rslv {

		void perform_resolve_hits(const str_vec &,
		                          std::istream & = std::cin,
		                          std::ostream & = std::cout,
		                          std::ostream & = std::cerr);

		void perform_resolve_hits(const crh_options &,
		                          std::istream & = std::cin,
		                          std::ostream & = std::cout,
		                          std::ostream & = std::cerr);

		void perform_resolve_hits(const crh_spec &,
		                          std::istream & = std::cin,
		                          std::ostream & = std::cout);

	} // namespace rslv
} // namespace cath

#endif
