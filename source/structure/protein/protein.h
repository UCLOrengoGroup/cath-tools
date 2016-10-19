/// \file
/// \brief The protein class header

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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_PROTEIN_PROTEIN_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_PROTEIN_PROTEIN_H

#include <boost/lexical_cast.hpp>
#include <boost/range/sub_range.hpp>

#include "common/temp_check_offset_1.h"
#include "exception/invalid_argument_exception.h"
#include "structure/geometry/coord.h"
#include "structure/protein/amino_acid.h"
#include "structure/protein/residue.h"
#include "structure/structure_type_aliases.h"

#include <iosfwd>
#include <string>
#include <cstddef>

namespace cath { class residue; }
namespace cath { class residue_name; }
namespace cath { class sec_struc; }

namespace cath {

	/// \brief The data on a protein as grabbed from the WOLF file and sec file
	class protein final {
	private:
		/// \brief TODOCUMENT
		std::string title;

		/// \brief TODOCUMENT
		residue_vec residues;

		/// \brief TODOCUMENT
		sec_struc_vec sec_strucs;

	private:
		void check_residue_index_is_valid(const size_t &) const;
		void check_sec_struc_is_valid(const size_t &) const;

	public:
		using iterator         = residue_vec::iterator;
		using const_iterator   = residue_vec::const_iterator;
		using sec_struc_crange = boost::sub_range<const sec_struc_vec>;

		protein() = default;
		protein(const std::string &, const residue_vec &);
		void set_title(const std::string &);
		void set_residues(const residue_vec &);
		void set_sec_strucs(const sec_struc_vec &);

		std::string get_title() const;

		inline residue & get_residue_ref_of_index(const size_t &);
		inline const residue & get_residue_ref_of_index(const size_t &) const;

		sec_struc & get_sec_struc_ref_of_index(const size_t &);
		const sec_struc & get_sec_struc_ref_of_index(const size_t &) const;

		inline size_t get_length() const;
		size_t get_num_sec_strucs() const;

		iterator begin();
		iterator end();

		const_iterator begin() const;
		const_iterator end() const;

		sec_struc_crange get_sec_strucs() const;
	};

	/// \brief TODOCUMENT
	inline void protein::check_residue_index_is_valid(const size_t &arg_index ///< TODOCUMENT
	                                                  ) const {
#ifndef NDEBUG
		if ( arg_index >= get_length() ) {
			BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception(
				"Residue index "
				+ boost::lexical_cast<std::string>( arg_index    )
				+ " is out of range in protein of length "
				+ boost::lexical_cast<std::string>( get_length() )
			));
		}
#else
		boost::ignore_unused( arg_index );
#endif
	}

	/// \brief TODOCUMENT
	inline void protein::check_sec_struc_is_valid(const size_t &arg_index ///< TODOCUMENT
	                                              ) const {
#ifndef NDEBUG
		if (arg_index >= get_num_sec_strucs()) {
			BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception(
				"Secondary structure index "
				+ boost::lexical_cast<std::string>( arg_index            )
				+ " is out of range in protein with "
				+ boost::lexical_cast<std::string>( get_num_sec_strucs() )
				+ " secondary structures"
			));
		}
#else
		boost::ignore_unused( arg_index );
#endif
	}

	/// \brief TODOCUMENT
	inline residue & protein::get_residue_ref_of_index(const size_t &arg_index ///< TODOCUMENT
	                                                   ) {
		check_residue_index_is_valid(arg_index);
		return residues[arg_index];
	}

	/// \brief TODOCUMENT
	inline const residue & protein::get_residue_ref_of_index(const size_t &arg_index ///< TODOCUMENT
	                                                         ) const {
		check_residue_index_is_valid(arg_index);
		return residues[arg_index];
	}

	/// \brief TODOCUMENT
	inline size_t protein::get_length() const {
		return residues.size();
	}

	residue & get_residue_ref_of_index__offset_1(protein &,
	                                             const size_t &);

	const residue & get_residue_ref_of_index__offset_1(const protein &,
	                                                   const size_t &);

	protein build_protein(const residue_vec &);

	protein build_protein(const residue_vec &,
	                      const sec_struc_vec &);

	residue_name get_pdb_residue_name_of_index(const protein &,
	                                           const size_t &);

	std::string get_pdb_residue_name_string_of_index(const protein &,
	                                                 const size_t &);

	size_t get_index_of_pdb_residue_name(const protein &,
	                                     const residue_name &);

	amino_acid_vec get_amino_acid_list(const protein &);

	amino_acid get_amino_acid_of_index(const protein &,
	                                   const size_t &);

	char get_amino_acid_letter_of_index(const protein &,
	                                    const size_t &);

	residue_name_vec get_residue_names(const protein &);

	void wipe_sec_strucs_of_residues(protein &);

	void label_residues_with_sec_strucs(protein &,
	                                    std::ostream &);

	geom::coord calculate_inter_sec_struc_vector(const protein &,
	                                             const size_t &,
	                                             const size_t &);

	geom::coord view_vector(const protein &,
	                        const size_t &,
	                        const size_t &);


	/// \brief TODOCUMENT
	///
	/// This is a temporary function during the transition away from using offsets of 1
	///
	/// \relates protein
	inline residue & get_residue_ref_of_index__offset_1(protein      &arg_protein,       ///< The protein to query
	                                                    const size_t &arg_index_offset_1 ///< The index (offset 1)
	                                                    ) {
		check_offset_1( arg_index_offset_1 );
		return arg_protein.get_residue_ref_of_index(arg_index_offset_1 - 1);
	}

	/// \brief TODOCUMENT
	///
	/// This is a temporary function during the transition away from using offsets of 1
	///
	/// \relates protein
	inline const residue & get_residue_ref_of_index__offset_1(const protein &arg_protein,       ///< The protein to query
	                                                          const size_t  &arg_index_offset_1 ///< The index (offset 1)
	                                                          ) {
		check_offset_1( arg_index_offset_1 );
		return arg_protein.get_residue_ref_of_index(arg_index_offset_1 - 1);
	}

} // namespace cath
#endif
