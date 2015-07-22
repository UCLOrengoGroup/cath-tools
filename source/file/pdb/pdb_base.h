/// \file
/// \brief The pdb_base class header

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

#ifndef PDB_BASE_H_INCLUDED
#define PDB_BASE_H_INCLUDED

#include <boost/filesystem.hpp>

#include "common/type_aliases.h"
#include "structure/chain_label.h"
#include "structure/structure_type_aliases.h"

#include <iosfwd>
#include <vector>

namespace cath { namespace geom { class coord; } }
namespace cath { namespace geom { class rotation; } }

namespace cath {
	namespace file {

		/// \brief TODOCUMENT
		///
		/// A pdb_base represents a parsed PDB file more than the protein itself.
		///
		/// \todo Possibly change this to always read on construction and only read on construction?
		///
		/// \todo This interface is currently the natural one that fits bioplib_pdb but it should
		///       be changed to be the interface that would be designed from scratch. Specifics follow:
		///
		/// \todo Specifically: pdb_residue and pdb_atom should also be made concrete instances of ABCs.
		///
		/// \todo Specifically: via, get_residue_names_of_first_chain(), get_residue_cref_of_index() etc, this interface
		///       should do provide the access that you would expect. bioplib_pdb should do the hard (possibly inefficient)
		///       work to meet this interface rather than vice versa.
		///
		/// \todo Specifically things like reading, writing, setting the chain label, getting a list of residue names,
		///       getting the CA coord of a residue of specified index etc should all be done via this interface.
		///
		/// \todo bioplib_pdb should keep its own methods that wrap ReadPDB() and WritePDB() and tests should
		///       test that these give identical output.
		///
		/// \todo get_residue_names_of_first_chain() should be get_residue_names_of_chain(const char &) and there should be an
		///       easy way to get the chain_code of the first chain.
		class pdb_base {
		private:
			virtual void do_read_file(const boost::filesystem::path &) = 0;
			virtual void do_append_to_file(const boost::filesystem::path &) const = 0;
			virtual void do_set_chain_label(const chain_label &) = 0;
			virtual residue_name_vec do_get_residue_names_of_first_chain__backbone_unchecked() const = 0;
			virtual geom::coord do_get_residue_ca_coord_of_index__backbone_unchecked(const size_t &) const = 0;
			virtual size_t do_get_num_atoms() const = 0;

			virtual void do_rotate(const geom::rotation &) = 0;
			virtual void do_add(const geom::coord &) = 0;
			virtual void do_subtract(const geom::coord &) = 0;

		public:
			virtual ~pdb_base() noexcept = default;

			void read_file(const boost::filesystem::path &);
			void append_to_file(const boost::filesystem::path &) const;
			void set_chain_label(const chain_label &);
			residue_name_vec get_residue_names_of_first_chain__backbone_unchecked() const;
			geom::coord get_residue_ca_coord_of_index__backbone_unchecked(const size_t &) const;
			size_t get_num_atoms() const;

			void rotate(const geom::rotation &);
			void operator+=(const geom::coord &);
			void operator-=(const geom::coord &);

			static const std::string PDB_RECORD_STRING_TER;
			static const size_t      MAX_NUM_PDB_COLS;
		};
	}
}

#endif