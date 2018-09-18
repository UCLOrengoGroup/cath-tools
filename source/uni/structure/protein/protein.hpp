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

#ifndef _CATH_TOOLS_SOURCE_UNI_STRUCTURE_PROTEIN_PROTEIN_H
#define _CATH_TOOLS_SOURCE_UNI_STRUCTURE_PROTEIN_PROTEIN_H

#include <boost/lexical_cast.hpp>
#include <boost/range/sub_range.hpp>

#include "chopping/chopping_type_aliases.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/temp_check_offset_1.hpp"
#include "file/name_set/name_set.hpp"
#include "structure/geometry/coord.hpp"
#include "structure/protein/amino_acid.hpp"
#include "structure/protein/residue.hpp"
#include "structure/structure_type_aliases.hpp"

#include <cstddef>
#include <iosfwd>
#include <string>

namespace cath { class residue; }
namespace cath { class residue_id; }
namespace cath { class sec_struc; }

namespace cath {

	/// \brief The data on a protein as grabbed from the WOLF file and sec file
	class protein final {
	private:
		/// \brief TODOCUMENT
		file::name_set the_name_set;

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
		protein(file::name_set,
		        residue_vec);
		protein & set_name_set(file::name_set);
		protein & set_residues(residue_vec);
		protein & set_sec_strucs(sec_struc_vec);

		file::name_set & get_name_set();
		const file::name_set & get_name_set() const;

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
	inline void protein::check_residue_index_is_valid(const size_t &prm_index ///< TODOCUMENT
	                                                  ) const {
#ifndef NDEBUG
		if ( prm_index >= get_length() ) {
			BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception(
				"Residue index "
				+ boost::lexical_cast<std::string>( prm_index    )
				+ " is out of range in protein of length "
				+ boost::lexical_cast<std::string>( get_length() )
			));
		}
#else
		boost::ignore_unused( prm_index );
#endif
	}

	/// \brief TODOCUMENT
	inline void protein::check_sec_struc_is_valid(const size_t &prm_index ///< TODOCUMENT
	                                              ) const {
#ifndef NDEBUG
		if (prm_index >= get_num_sec_strucs()) {
			BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception(
				"Secondary structure index "
				+ boost::lexical_cast<std::string>( prm_index            )
				+ " is out of range in protein with "
				+ boost::lexical_cast<std::string>( get_num_sec_strucs() )
				+ " secondary structures"
			));
		}
#else
		boost::ignore_unused( prm_index );
#endif
	}

	/// \brief TODOCUMENT
	inline residue & protein::get_residue_ref_of_index(const size_t &prm_index ///< TODOCUMENT
	                                                   ) {
		check_residue_index_is_valid(prm_index);
		return residues[prm_index];
	}

	/// \brief TODOCUMENT
	inline const residue & protein::get_residue_ref_of_index(const size_t &prm_index ///< TODOCUMENT
	                                                         ) const {
		check_residue_index_is_valid(prm_index);
		return residues[prm_index];
	}

	/// \brief TODOCUMENT
	inline size_t protein::get_length() const {
		return residues.size();
	}

	residue & get_residue_ref_of_index__offset_1(protein &,
	                                             const size_t &);

	const residue & get_residue_ref_of_index__offset_1(const protein &,
	                                                   const size_t &);

	protein build_protein(residue_vec);

	protein build_protein(residue_vec,
	                      sec_struc_vec);

	residue_id get_pdb_residue_id_of_index(const protein &,
	                                       const size_t &);

	std::string get_pdb_residue_id_string_of_index(const protein &,
	                                               const size_t &);

	size_t get_index_of_pdb_residue_id(const protein &,
	                                   const residue_id &);

	amino_acid_vec get_amino_acid_list(const protein &);

	amino_acid get_amino_acid_of_index(const protein &,
	                                   const size_t &);

	char get_amino_acid_letter_of_index_tolerantly(const protein &,
	                                               const size_t &);

	residue_id_vec get_residue_ids(const protein &);

	void set_accessibilities(protein &,
	                         const size_vec &);

	void set_accessibilities(protein &,
	                         const doub_vec &);

	void set_sec_struc_types(protein &,
	                         const sec_struc_type_vec &);

	void wipe_sec_strucs_of_residues(protein &);

	void label_residues_with_sec_strucs(protein &,
	                                    const ostream_ref_opt &);

	geom::coord calculate_inter_sec_struc_vector(const protein &,
	                                             const size_t &,
	                                             const size_t &);

	geom::coord view_vector(const protein &,
	                        const size_t &,
	                        const size_t &);

	size_vec get_indices_of_residues_within_regions(const protein &,
	                                                const chop::region_vec_opt &);

	void restrict_to_regions(protein &,
	                         const chop::region_vec_opt &);

	protein restrict_to_regions_copy(protein,
	                                 const chop::region_vec_opt &);

	std::string get_domain_or_specified_or_name_from_acq(const protein &);

	/// \brief TODOCUMENT
	///
	/// This is a temporary function during the transition away from using offsets of 1
	///
	/// \relates protein
	inline residue & get_residue_ref_of_index__offset_1(protein      &prm_protein,       ///< The protein to query
	                                                    const size_t &prm_index_offset_1 ///< The index (offset 1)
	                                                    ) {
		check_offset_1( prm_index_offset_1 );
		return prm_protein.get_residue_ref_of_index(prm_index_offset_1 - 1);
	}

	/// \brief TODOCUMENT
	///
	/// This is a temporary function during the transition away from using offsets of 1
	///
	/// \relates protein
	inline const residue & get_residue_ref_of_index__offset_1(const protein &prm_protein,       ///< The protein to query
	                                                          const size_t  &prm_index_offset_1 ///< The index (offset 1)
	                                                          ) {
		check_offset_1( prm_index_offset_1 );
		return prm_protein.get_residue_ref_of_index(prm_index_offset_1 - 1);
	}

} // namespace cath
#endif
