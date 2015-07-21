/// \file
/// \brief The pdb_atom class header

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

#ifndef PDB_ATOM_H_INCLUDED
#define PDB_ATOM_H_INCLUDED

#include <boost/tuple/tuple.hpp>

#include "file/pdb/pdb_record.h"
#include "structure/chain_label.h"
#include "structure/geometry/coord.h"
#include "structure/protein/amino_acid.h"
#include "structure/residue_name.h"

#include <iosfwd>

namespace cath { namespace geom { class rotation; } }

namespace cath {
	namespace file {

		/// \brief TODOCUMENT
		class pdb_atom final {
		private:
			/// \brief TODOCUMENT
			pdb_record  record_type;

			/// \brief TODOCUMENT
			size_t      atom_serial;

			/// \brief The untrimmed string describing the element type of this atom.
			///
			/// This is whitespace-untrimmed so that the correct whitespace can be returned whilst writing
			std::string element_type_untrimmed;

			/// \brief The trimmed string describing the element type of this atom.
			std::string element_type;

			/// \brief TODOCUMENT
			char        alt_locn;

			/// \brief TODOCUMENT
			amino_acid  the_amino_acid;

			/// \brief TODOCUMENT
			geom::coord atom_coord;

			/// \brief TODOCUMENT
			double      occupancy;

			/// \brief TODOCUMENT
			double      temp_factor;

		public:
			pdb_atom(const pdb_record &,
			         const size_t &,
			         const std::string &,
			         const char &,
			         const amino_acid &,
			         const geom::coord &,
			         const double &,
			         const double &);

			const pdb_record & get_record_type() const;
			const size_t & get_atom_serial() const;
			const std::string & get_element_type_untrimmed() const;
			const std::string & get_element_type() const;
			const char & get_alt_locn() const;
			const amino_acid & get_amino_acid() const;
			const geom::coord & get_coord() const;
			const double & get_occupancy() const;
			const double & get_temp_factor() const;

			void rotate(const geom::rotation &);
			void operator+=(const geom::coord &);
			void operator-=(const geom::coord &);

			static const std::string PDB_ID_NITROGEN;
			static const std::string PDB_ID_CARBON_ALPHA;
			static const std::string PDB_ID_CARBON;

			static const std::string PDB_ID_CARBON_BETA;
			static const std::string PDB_ID_OXYGEN;
		};

		char get_amino_acid_letter(const pdb_atom &);
		std::string get_amino_acid_code(const pdb_atom &);
		std::string get_amino_acid_name(const pdb_atom &);
		bool is_pdb_record_of_type(const std::string &,
		                           const pdb_record &);
		std::string pdb_record_parse_problem(const std::string &);

		using chain_resname_atom_tuple = std::tuple<chain_label, residue_name, pdb_atom>;

		chain_resname_atom_tuple parse_pdb_atom_record(std::string);

		std::ostream & write_pdb_file_entry(std::ostream &,
		                                    const chain_label &,
		                                    const residue_name &,
		                                    const pdb_atom &);
		std::ostream & operator<<(std::ostream &,
		                          const pdb_atom &);
	}
}

#endif
