/// \file
/// \brief The hmmer_hmmsearch_out class header

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

#ifndef HMMER_HMMSEARCH_OUT_H_INCLUDED
#define HMMER_HMMSEARCH_OUT_H_INCLUDED

#include <boost/filesystem/path.hpp>

#include "common/type_aliases.h"
#include "resolve_hits/resolve_hits_type_aliases.h"

namespace cath { namespace rslv { class read_and_process_mgr; } }

namespace cath {
	namespace rslv {

		void parse_hmmsearch_out_file(read_and_process_mgr &,
		                              const boost::filesystem::path &,
		                              const bool &,
		                              const residx_t &);

		void parse_hmmsearch_out(read_and_process_mgr &,
		                         std::istream &,
		                         const bool &,
		                         const residx_t &);

	}
}

#endif
