/// \file
/// \brief The parse_hmmer_out class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_FILE_PARSE_HMMER_OUT_HPP
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_FILE_PARSE_HMMER_OUT_HPP

#include <boost/filesystem/path.hpp>

#include "common/type_aliases.hpp"
#include "resolve_hits/file/hmmer_format.hpp"
#include "resolve_hits/resolve_hits_type_aliases.hpp"

namespace cath { namespace rslv { class calc_hit_list; } }
namespace cath { namespace rslv { class read_and_process_mgr; } }

namespace cath {
	namespace rslv {

		void parse_hmmer_out_file(read_and_process_mgr &,
		                          const boost::filesystem::path &,
		                          const hmmer_format &,
		                          const bool &,
		                          const seq::residx_t &,
		                          const bool &);

		void parse_hmmer_out(read_and_process_mgr &,
		                     std::istream &,
		                     const hmmer_format &,
		                     const bool &,
		                     const seq::residx_t &,
		                     const bool &);

		str_calc_hit_list_pair_vec parse_hmmer_out_file(const boost::filesystem::path &,
		                                                const hmmer_format &,
		                                                const bool &,
		                                                const seq::residx_t &,
		                                                const bool &);

	} // namespace rslv
} // namespace cath

#endif
