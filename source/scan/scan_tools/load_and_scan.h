/// \file
/// \brief The load_and_scan class header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef LOAD_AND_SCAN_H_INCLUDED
#define LOAD_AND_SCAN_H_INCLUDED

#include "common/clone/clone_ptr.h"
#include "scan/detail/scan_type_aliases.h"
#include "scan/scan_tools/scan_metrics.h"
#include "scan/scan_tools/scan_type.h"
#include "structure/protein/protein_loader/protein_list_loader.h"

namespace cath { namespace scan { class load_and_scan_metrics; } }

namespace cath {
	namespace scan {

		/// \brief TODOCUMENT
		class load_and_scan final {
		private:
			/// \brief TODOCUMENT
			protein_list_loader query_protein_loader;

			/// \brief TODOCUMENT
			protein_list_loader match_protein_loader;

			/// \brief TODOCUMENT
			common::clone_ptr<scan_type> scan_ptr;

			/// \brief TODOCUMENT
			boost::optional<protein_list> query_proteins;

			/// \brief TODOCUMENT
			boost::optional<protein_list> match_proteins;

			/// \brief TODOCUMENT
			hrc_duration_opt load_files_duration;

			/// \brief TODOCUMENT
			boost::optional<scan_metrics> the_scan_metrics;

			void perform_load();
			void perform_scan();

		public:
			load_and_scan(const protein_list_loader &,
			              const protein_list_loader &,
			              const scan_type &);

			const protein_list & get_query_proteins() const;
			const protein_list & get_match_proteins() const;
			load_and_scan_metrics get_load_and_scan_metrics() const;
		};

	}
}

#endif
