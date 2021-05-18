/// \file
/// \brief The pdb_residue class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_PDB_RESIDUE_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_PDB_RESIDUE_HPP

#include <optional>
#include <vector>

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/operators.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/algorithm/find_if.hpp>

#include "cath/biocore/chain_label.hpp"
#include "cath/biocore/residue_id.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/file/file_type_aliases.hpp"
#include "cath/file/pdb/pdb_atom.hpp"
#include "cath/file/pdb/residue_makeup.hpp"
#include "cath/structure/geometry/coord_linkage.hpp"
#include "cath/structure/structure_type_aliases.hpp"

using namespace ::cath::common::literals;

namespace cath { class amino_acid; }
namespace cath { namespace file { class pdb_atom; } }
namespace cath { namespace geom { class coord; } }
namespace cath { namespace geom { class rotation; } }
namespace cath { class residue; }

namespace cath {
	namespace file {

		/// \brief Represent a batch of pdb_atoms all belonging to the same residue
		///
		/// The amino acid cannot be stored here, rather than in the pdb_atom entries because
		/// some residues (eg residues 22 and 25 in 1cbnA) have different amino acid values in
		/// the same residue. Storing that the residue level would mean that information could
		/// not be reproduced when outputting the pdb.
		class pdb_residue final : boost::additive<pdb_residue, geom::coord> {
		private:
			/// \brief The name of the residue
			residue_id   the_residue_id;

			/// \brief The pdb_atoms making up this residue
			pdb_atom_vec atoms;

			/// \brief An array of the indices of the core element atoms (N/CA/C/CB/O) within atoms, above,
			///        or nullopt for any core element atoms not present in this residue
			std::array<size_opt, 5> core_atom_indices;

			static constexpr size_t get_core_atom_index_ref_index(const coarse_element_type &);
			const ::std::optional<size_t> & get_core_atom_index_ref(const coarse_element_type &) const;
			static std::array<::std::optional<size_t>, 5> make_core_atom_indices(const pdb_atom_vec &);

		public:
			/// \brief Type alias for the const_iterator type for the range of atoms
			using const_iterator = pdb_atom_vec::const_iterator;

			pdb_residue(residue_id,
			            pdb_atom_vec);

			const residue_id & get_residue_id() const;
			bool empty() const;
			size_t get_num_atoms() const;
			const pdb_atom & get_atom_cref_of_index(const size_t &) const;

			bool has_nitrogen() const;
			bool has_carbon_alpha() const;
			bool has_carbon() const;
			bool has_carbon_beta() const;
			bool has_oxygen() const;

			const pdb_atom & get_nitrogen() const;
			const pdb_atom & get_carbon_alpha() const;
			const pdb_atom & get_carbon() const;
			const pdb_atom & get_carbon_beta() const;
			const pdb_atom & get_oxygen() const;

			const amino_acid & get_amino_acid() const;

			pdb_residue & set_chain_label(const chain_label &);

			pdb_residue & rotate(const geom::rotation &);
			pdb_residue & operator+=(const geom::coord &);
			pdb_residue & operator-=(const geom::coord &);

			const_iterator begin() const;
			const_iterator end() const;
		};

		const chain_label & get_chain_label(const pdb_residue &);
		const residue_name & get_residue_name(const pdb_residue &);

		const geom::coord & get_nitrogen_coord(const pdb_residue &);
		const geom::coord & get_carbon_alpha_coord(const pdb_residue &);
		const geom::coord & get_carbon_coord(const pdb_residue &);
		const geom::coord & get_carbon_beta_coord(const pdb_residue &);
		const geom::coord & get_oxygen_coord(const pdb_residue &);

		char_opt get_letter_if_amino_acid(const pdb_residue &);
		char get_amino_acid_letter_tolerantly(const pdb_residue &);
		char_3_arr get_amino_acid_code(const pdb_residue &);
		std::string get_amino_acid_code_string(const pdb_residue &);
		std::string_view get_amino_acid_name( const pdb_residue & );

		residue_makeup contains_non_proper_amino_acids(const pdb_residue &);

		bool is_backbone_complete(const pdb_residue &);

		geom::coord fake_carbon_beta_coord_of_residue(const pdb_residue &);
		geom::coord get_or_predict_carbon_beta_coord_of_residue(const pdb_residue &);
		geom::rotation get_ssap_frame_of_residue(const pdb_residue &);

		double max_dist_from_coord(const pdb_residue &,
		                           const geom::coord &);

		std::ostream & write_pdb_file_entry(std::ostream &,
		                                    const pdb_residue &);

		geom::angle_angle_pair<double> get_psi_of_this_and_phi_of_next(const pdb_residue &,
		                                                               const pdb_residue &);

		residue build_residue_of_pdb_residue(const pdb_residue &,
		                                     const geom::doub_angle &,
		                                     const geom::doub_angle &,
		                                     const size_t &);

		bool dssp_will_skip_residue(const pdb_residue &);

		geom::coord_vec get_all_coords(const pdb_residue_vec &);

		geom::coord_coord_linkage_pair_vec get_all_coords_with_linkage(const pdb_residue_vec &);

		std::ostream & operator<<(std::ostream &,
		                          const pdb_residue_vec &);

		std::ostream & operator<<(std::ostream &,
		                          const pdb_residue &);

		//// \brief Get the core_atom_indices index corresponding to the specified coarse_element_type
		inline constexpr size_t pdb_residue::get_core_atom_index_ref_index(const coarse_element_type &prm_element_type ///< The coarse_element_type to be accessed in core_atom_indices
		                                                                   ) {
			return
				( prm_element_type == coarse_element_type::NITROGEN     ) ? 0_z :
				( prm_element_type == coarse_element_type::CARBON_ALPHA ) ? 1_z :
				( prm_element_type == coarse_element_type::CARBON       ) ? 2_z :
				( prm_element_type == coarse_element_type::CARBON_BETA  ) ? 3_z :
				( prm_element_type == coarse_element_type::OXYGEN       ) ? 4_z :
					throw std::invalid_argument("Cannot get_core_atom_index_ref() for non-core element type");
		}

		/// \brief Get a const reference to the core atom index corresponding to the specified coarse_element_type
		///
		/// The returned value is the index within atoms of the relevant atom or nullopt if no such is present
		inline const ::std::optional<size_t> & pdb_residue::get_core_atom_index_ref(const coarse_element_type &prm_element_type ///< The coarse_element_type to be accessed in core_atom_indices
		                                                                            ) const {
			return core_atom_indices.at( get_core_atom_index_ref_index( prm_element_type ) );
		}

		/// \brief Make an array of the core element atoms' indices within the specified pdb_atoms
		///        (or nullopt for elements that have no corresponding atoms)
		///
		/// The positions of the elements in the array are defined by get_core_atom_index_ref_index()
		///
		/// This is arranged so that, in the cases where there are multiple atoms of the same type,
		/// it will replicate DSSP's implicit policy :
		///
		///  * Prefer an atom that meets alt_locn_is_dssp_accepted() (ie altlocn of ' ' or 'A'), then
		///  * Prefer an atom nearer the back
		inline std::array<::std::optional<size_t>, 5> pdb_residue::make_core_atom_indices(const pdb_atom_vec &prm_pdb_atoms ///< The list of atoms to index
		                                                                                  ) {
			std::array<::std::optional<size_t>, 5> results;
			for (const size_t &atom_ctr : common::indices( prm_pdb_atoms.size() ) ) {

				const pdb_atom            &the_atom = prm_pdb_atoms[ atom_ctr ];
				const coarse_element_type &element  = get_coarse_element_type( the_atom );
				if ( element != coarse_element_type::NON_CORE ) {
					::std::optional<size_t>   &result   = results.at( get_core_atom_index_ref_index( element ) );
					if ( ! result || alt_locn_is_dssp_accepted( the_atom ) || ! alt_locn_is_dssp_accepted( prm_pdb_atoms[ *result ] ) ) {
						result = atom_ctr;
					}
				}
			}
			return results;
		}

		/// \brief Ctor for pdb_residue
		///
		inline pdb_residue::pdb_residue(residue_id   prm_residue_id, ///< The name of the residue
		                                pdb_atom_vec prm_atoms       ///< The pdb_atoms making up this residue
		                                ) : the_residue_id   { std::move( prm_residue_id     ) },
		                                    atoms            { std::move( prm_atoms          ) },
		                                    core_atom_indices( make_core_atom_indices( atoms ) ) //< Don't change these brackets to braces - it breaks the build on the older Clang on Travis-CI
		                                    {
		}

		/// \brief Getter for the residue name
		inline const residue_id & pdb_residue::get_residue_id() const {
			return the_residue_id;
		}

		/// \brief Return whether the pdb_residue is empty (ie contains no pdb_atom entries)
		inline bool pdb_residue::empty() const {
			return atoms.empty();
		}

		/// \brief Get the number of atoms held in this residue
		inline size_t pdb_residue::get_num_atoms() const {
			return atoms.size();
		}

		/// \brief Get the pdb_atom for the specified index
		inline const pdb_atom & pdb_residue::get_atom_cref_of_index(const size_t &prm_index ///< The index of the pbd_atom to return
		                                                            ) const {
			return atoms[ prm_index ];
		}

		/// \brief Get whether this has a nitrogen atom
		inline bool pdb_residue::has_nitrogen() const {
			return static_cast<bool>( get_core_atom_index_ref( coarse_element_type::NITROGEN     ) );
		}

		/// \brief Get whether this has a carbon_alpha atom
		inline bool pdb_residue::has_carbon_alpha() const {
			return static_cast<bool>( get_core_atom_index_ref( coarse_element_type::CARBON_ALPHA ) );
		}

		/// \brief Get whether this has a carbon atom
		inline bool pdb_residue::has_carbon() const {
			return static_cast<bool>( get_core_atom_index_ref( coarse_element_type::CARBON       ) );
		}

		/// \brief Get whether this has a carbon_beta atom
		inline bool pdb_residue::has_carbon_beta() const {
			return static_cast<bool>( get_core_atom_index_ref( coarse_element_type::CARBON_BETA  ) );
		}

		/// \brief Get whether this has a oxygen atom
		inline bool pdb_residue::has_oxygen() const {
			return static_cast<bool>( get_core_atom_index_ref( coarse_element_type::OXYGEN       ) );
		}

		/// \brief Get the nitrogen atom
		///
		/// \pre `has_nitrogen()`
		inline const pdb_atom & pdb_residue::get_nitrogen() const {
			return get_atom_cref_of_index( *get_core_atom_index_ref( coarse_element_type::NITROGEN     ) );
		}

		/// \brief Get the carbon_alpha atom
		///
		/// \pre `has_carbon_alpha()`
		inline const pdb_atom & pdb_residue::get_carbon_alpha() const {
			return get_atom_cref_of_index( *get_core_atom_index_ref( coarse_element_type::CARBON_ALPHA ) );
		}

		/// \brief Get the carbon atom
		///
		/// \pre `has_carbon()`
		inline const pdb_atom & pdb_residue::get_carbon() const {
			return get_atom_cref_of_index( *get_core_atom_index_ref( coarse_element_type::CARBON       ) );
		}

		/// \brief Get the carbon_beta atom
		///
		/// \pre `has_carbon_beta()`
		inline const pdb_atom & pdb_residue::get_carbon_beta() const {
			return get_atom_cref_of_index( *get_core_atom_index_ref( coarse_element_type::CARBON_BETA  ) );
		}

		/// \brief Get the oxygen atom
		///
		/// \pre `has_oxygen()`
		inline const pdb_atom & pdb_residue::get_oxygen() const {
			return get_atom_cref_of_index( *get_core_atom_index_ref( coarse_element_type::OXYGEN       ) );
		}

		/// \brief Return the amino acid for this PDB residue
		///
		/// \pre `! empty()` else this throws an invalid_argument_exception
		///
		/// This returns the amino acid for the first present core atom
		/// or the amino acid of the first atom if there are no core atoms
		inline const amino_acid & pdb_residue::get_amino_acid() const {
			// Find the first present core atom
			const auto first_present_core_itr = boost::range::find_if(
				core_atom_indices,
				[] (const size_opt &x) { return static_cast<bool>( x ); }
			);

			// If there was one, return its amino acid
			if ( first_present_core_itr != ::std::cend( core_atom_indices ) ) {
				return get_atom_cref_of_index( first_present_core_itr->value() ).get_amino_acid();
			}

			// If there are no atoms, throw an exception
			if ( empty() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Can't get amino acid for empty pdb_residue"));
			}

			// Otherwise, just return the first atom's amino acid
			return get_atom_cref_of_index( 0 ).get_amino_acid();
		}

		/// \brief Setter for the chain label
		inline pdb_residue & pdb_residue::set_chain_label(const chain_label &prm_chain_label ///< The new chain label to set
		                                                  ) {
			the_residue_id = residue_id{ prm_chain_label, get_residue_name( *this ) };
			return *this;
		}

		/// \brief Rotate each of the atoms in this residue by the specified rotation
		inline pdb_residue & pdb_residue::rotate(const geom::rotation &prm_rotation ///< The rotation to apply to each of the atoms in this residue
		                                         ) {
			for (pdb_atom &my_pdb_atom : atoms) {
				my_pdb_atom.rotate( prm_rotation );
			}
			return *this;
		}

		/// \brief Add the specified coord from each of the atoms in this residue
		inline pdb_residue & pdb_residue::operator+=(const geom::coord &prm_coord ///< The coord to add to each of the atoms in this residue
		                                             ) {
			for (pdb_atom &my_pdb_atom : atoms) {
				my_pdb_atom += prm_coord;
			}
			return *this;
		}

		/// \brief Subtract the specified coord from each of the atoms in this residue
		inline pdb_residue & pdb_residue::operator-=(const geom::coord &prm_coord ///< The coord to subtract from each of the atoms in this residue
		                                             ) {
			for (pdb_atom &my_pdb_atom : atoms) {
				my_pdb_atom -= prm_coord;
			}
			return *this;
		}

		/// \brief Standard const begin() method for the range of pdb_atoms
		inline auto pdb_residue::begin() const -> const_iterator {
			return ::std::cbegin( atoms );
		}

		/// \brief Standard const end() method for the range of pdb_atoms
		inline auto pdb_residue::end() const -> const_iterator {
			return ::std::cend  ( atoms );
		}

		/// \brief Get the chain_label from the specified pdb_residue
		inline const chain_label & get_chain_label(const pdb_residue &prm_pdb_residue ///< The residue to query
		                                           ) {
			return prm_pdb_residue.get_residue_id().get_chain_label();
		}

		/// \brief Get the residue_name from the specified pdb_residue
		inline const residue_name & get_residue_name(const pdb_residue &prm_pdb_residue ///< The residue to query
		                                             ) {
			return prm_pdb_residue.get_residue_id().get_residue_name();
		}

		/// \brief Get the coord for the specified residue's nitrogen atom
		inline const geom::coord & get_nitrogen_coord(const pdb_residue &prm_pdb_residue ///< The residue to query
		                                              ) {
			return prm_pdb_residue.get_nitrogen().get_coord();
		}

		/// \brief Get the coord for the specified residue's carbon_alpha atom
		inline const geom::coord & get_carbon_alpha_coord(const pdb_residue &prm_pdb_residue ///< The residue to query
		                                                  ) {
			return prm_pdb_residue.get_carbon_alpha().get_coord();
		}

		/// \brief Get the coord for the specified residue's carbon atom
		inline const geom::coord & get_carbon_coord(const pdb_residue &prm_pdb_residue ///< The residue to query
		                                            ) {
			return prm_pdb_residue.get_carbon().get_coord();
		}

		/// \brief Get the coord for the specified residue's carbon_beta atom
		inline const geom::coord & get_carbon_beta_coord(const pdb_residue &prm_pdb_residue ///< The residue to query
		                                                 ) {
			return prm_pdb_residue.get_carbon_beta().get_coord();
		}

		/// \brief Get the coord for the specified residue's oxygen atom
		inline const geom::coord & get_oxygen_coord(const pdb_residue &prm_pdb_residue ///< The residue to query
		                                            ) {
			return prm_pdb_residue.get_oxygen().get_coord();
		}

		// /// \brief Get an atom of the specified type from the specified residue
		// ///
		// /// This is organised so that, in the cases where there are multiple atoms of the same type,
		// /// it will replicate DSSP's implicit policy :
		// ///
		// ///  * Prefer an atom that meets alt_locn_is_dssp_accepted() (ie altlocn of ' ' or 'A'), then
		// ///  * Prefer an atom nearer the back
		// ///
		// /// \relates pdb_residue
		// inline const geom::coord & get_atom_of_id_of_residue(const pdb_residue         &prm_pdb_residue, ///< The residue from which to extract the atom coordinates
		//                                                      const coarse_element_type &prm_element      ///< The type of atom to retrieve
		//                                                      ) {
		// 	size_opt best_seen_atom_ctr;
		// 	for (const size_t &atom_ctr : common::indices( prm_pdb_residue.get_num_atoms() ) | boost::adaptors::reversed) {
		// 		const pdb_atom &the_atom = prm_pdb_residue.get_atom_cref_of_index( atom_ctr );
		// 		if ( get_coarse_element_type( the_atom ) == prm_element ) {
		// 			if ( alt_locn_is_dssp_accepted( the_atom ) ) {
		// 				return the_atom.get_coord();
		// 			}
		// 			if ( ! best_seen_atom_ctr ) {
		// 				best_seen_atom_ctr = atom_ctr;
		// 			}
		// 		}
		// 	}
		// 	if ( best_seen_atom_ctr ) {
		// 		return prm_pdb_residue.get_atom_cref_of_index( *best_seen_atom_ctr ).get_coord();
		// 	}
		// 	BOOST_THROW_EXCEPTION(common::invalid_argument_exception(
		// 		"Cannot find atom of type "
		// 		+ to_string( prm_element )
		// 		+ " within the "
		// 		+ boost::lexical_cast<std::string>( prm_pdb_residue.get_num_atoms() )
		// 		+ " atom(s) of residue "
		// 		+ boost::lexical_cast<std::string>( prm_pdb_residue.get_residue_name() )
		// 	));
		// 	return geom::ORIGIN_COORD; // To appease Eclipse's syntax parser
		// }

	} // namespace file
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_PDB_RESIDUE_HPP
