/// \file
/// \brief The bioplib_pdb class header

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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_BIOPLIB_FACADE_BIOPLIB_PDB_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_BIOPLIB_FACADE_BIOPLIB_PDB_H

#include "file/pdb/pdb_base.hpp"

extern "C" {
	/// \todo This hits a problem when trying to compile against bioplib >= v3.0
	///       The included files (pdb.hpp -> general.hpp -> deprecated.hpp -> macros.hpp) define:
	/// ~~~~~
	/// extern int isatty(int);
	/// ~~~~~
	///       but this conflicts with the one with a `throw ()` specification
	///       that gets included into C++ code from unistd.h.
	///       This problem's rather tricky. The include of the `throw ()` version
	///       needn't happen in within the include tree of pdb.hpp for the problem to occur.
	///       For now, just use bioplib <= v2.1.2
	#include <unistd.h>
	#include "bioplib/pdb.h"
}

#include <vector>

namespace cath { namespace geom { class coord; } }
namespace cath { namespace geom { class rotation; } }
namespace cath { namespace test { struct bioplib_pdb_test_suite_fixture; } }

namespace cath {
	namespace file {

		/// \brief A simple wrapper around a bioplib PDB to hold it with the RAII idiom
		///        and to provide the abilities that are used in the CATH algorithms.
		///
		/// This is an implementation of the ABC pdb_base
		class bioplib_pdb final : public pdb_base {
		private:
			/// \brief TODOCUMENT
			PDB * pdb_ptr = nullptr;

			/// \brief TODOCUMENT
			int natoms = 0;

			friend struct cath::test::bioplib_pdb_test_suite_fixture;

			void free_pdb_internals();
			void check_ptr() const;
			PDB * get_ptr();
			const PDB * get_ptr() const;
			size_t get_natoms() const;

			void do_read_file(const boost::filesystem::path &) final;
			void do_append_to_file(const boost::filesystem::path &) const final;
			void do_set_chain_label(const chain_label &) final;
			residue_id_vec do_get_residue_ids_of_first_chain__backbone_unchecked() const final;
			geom::coord do_get_residue_ca_coord_of_index__backbone_unchecked(const size_t &) const final;
			size_t do_get_num_atoms() const final;

			void do_rotate(const geom::rotation &) final;
			void do_add(const geom::coord &) final;
			void do_subtract(const geom::coord &) final;

		public:
			bioplib_pdb() = default;
			bioplib_pdb(const bioplib_pdb &);
			bioplib_pdb(bioplib_pdb &&) = default;
			bioplib_pdb & operator=(const bioplib_pdb &);
			bioplib_pdb & operator=(bioplib_pdb &&) = default;
			~bioplib_pdb() noexcept override;
		};

	} // namespace file
} // namespace cath

#endif
