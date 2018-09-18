/// \file
/// \brief The bifur_hbond_list class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_STRUCTURE_SEC_STRUC_CALC_DSSP_BIFUR_HBOND_LIST_HPP
#define _CATH_TOOLS_SOURCE_UNI_STRUCTURE_SEC_STRUC_CALC_DSSP_BIFUR_HBOND_LIST_HPP

#include <boost/optional.hpp>

#include "common/cpp14/cbegin_cend.hpp"

#include <utility>
#include <vector>

namespace cath {
	namespace sec {

		/// \brief Type to use for storing the index of an hbond partner residue
		using hbond_partner_t = unsigned int;

		/// \brief Type to use for storing the energy of an hbond
		using hbond_energy_t  = double;

		/// \brief Represent the details of an hbond from a given residue
		///        (called "_half" because the from-residue and type of bond (NH->CO or CO->NH) is
		///        unspecified, left to be implied by context)
		struct hbond_half final {
			/// \brief The index of the hbond partner residue
			hbond_partner_t index;

			/// \brief The energy of the hbond
			hbond_energy_t energy;

			/// The energy cutoff below which DSSP considers a candidate to be a true hbond
			static constexpr hbond_energy_t HBOND_ENERGY_CUTOFF = -0.5;
		};

		/// \brief Equality operator for hbond_half
		///
		/// \relates hbond_half
		inline constexpr bool operator==(const hbond_half &prm_lhs, ///< The first  hbond_half to compare
		                                 const hbond_half &prm_rhs  ///< The second hbond_half to compare
		                                 ) {
			return (
				prm_lhs.index  == prm_rhs.index
				&&
				prm_lhs.energy == prm_rhs.energy
			);
		}

		/// \brief Return whether the first hbond_half represents a stronger bond than the second
		///        or has the same strength but an earlier index
		///
		/// \relates hbond_half
		inline constexpr bool is_bondier_than(const hbond_half &prm_lhs, ///< The first  hbond_half to compare
		                                      const hbond_half &prm_rhs  ///< The second hbond_half to compare
		                                      ) {
			return (
				( prm_lhs.energy < prm_rhs.energy )
				||
				(
					( prm_lhs.energy == prm_rhs.energy )
					&&
					( prm_lhs.index  <  prm_rhs.index  )
				)
			);
		}

		/// \brief Type alias for an optional hbond_half
		using hbond_half_opt = boost::optional<hbond_half>;

		std::string to_string(const hbond_half_opt &);

		/// \brief Whether the specified hbond_half has a strong enough energy to be considered a true hbond
		inline bool is_bondy_enough(const hbond_half &prm_hbond_half ///< The hbond_half to consider
		                            ) {
			return prm_hbond_half.energy < hbond_half::HBOND_ENERGY_CUTOFF ;
		}

		/// \brief Whether the specified hbond_half_opt is present and has a strong enough energy to be considered a true hbond
		inline bool is_bondy_enough(const hbond_half_opt &prm_hbond_half_opt ///< The hbond_half_opt to consider
		                            ) {
			return ( prm_hbond_half_opt && is_bondy_enough( *prm_hbond_half_opt ) );
		}

		/// \brief Wipe the specified hbond_half_opt to none if 
		///         * (a) it isn't already and
		///         * (b) it isn't strong enough to be considered a true hbond
		inline void remove_not_bondy_enough(hbond_half_opt &prm_hbond_half_opt ///< The hbond_half_opt to examine (and potentially wipe)
		                                    ) {
			if ( prm_hbond_half_opt && ! is_bondy_enough( *prm_hbond_half_opt ) ) {
				prm_hbond_half_opt = boost::none;
			}
		}

		/// \brief Type alias for a pair of hbond_half_opts
		///
		/// This is used to represent a pair of hbond_halfs of the same type from the same residue
		/// and are kept in descending order of strength, ie the following invariants:
		///  * `! ( ! first && second )`
		///  * `! second || ! is_bondier_than( second, first )`
		using hbond_half_opt_pair = std::pair<hbond_half_opt, hbond_half_opt>;

		/// \brief Update an half_bond_opt_pair with a new hbond_half
		///
		/// If new hbond_half is active and is better then either of the entries in the hbond_half_opt_pair
		/// then insert it in the correct position
		inline void update_half_bond_pair(hbond_half_opt_pair &prm_hbond_pair, ///< The hbond_half_opt_pair to update
		                                  const hbond_half    &prm_hbond       ///< The new hbond_half with which to update the hbond_half_opt_pair
		                                  ) {
			// If the new one's better than the first, copy the first to the second and set
			// first to the new one
			if ( ! prm_hbond_pair.first || is_bondier_than( prm_hbond, *prm_hbond_pair.first ) ) {
				prm_hbond_pair.second = prm_hbond_pair.first;
				prm_hbond_pair.first = prm_hbond;
			}
			// Else if the new one's better than the second, overwrite the second with the new one
			else if ( ! prm_hbond_pair.second || is_bondier_than( prm_hbond, *prm_hbond_pair.second ) ) {
				prm_hbond_pair.second = prm_hbond;
			}
		}

		/// \brief Copy the specified hbond_half_opt_pair and update_half_bond_pair() and return the copy
		inline hbond_half_opt_pair update_half_bond_pair_copy(hbond_half_opt_pair  prm_hbond_pair, ///< The hbond_half_opt_pair from which to take a copy to update_half_bond_pair() and return
		                                                      const hbond_half    &prm_hbond       ///< The hbond_half with which the copy of the hbond_half_opt_pair should be updated
		                                                      ) {
			update_half_bond_pair( prm_hbond_pair, prm_hbond );
			return prm_hbond_pair;
		}

		/// \brief Wipe each half of the specified hbond_half_opt_pair to none if 
		///         * (a) it isn't already and
		///         * (b) it isn't strong enough to be considered a true hbond
		///
		/// \pre If only one bond is present, it must be the first
		///
		/// \pre If both bonds are present, the first mustn't be weaker
		inline void remove_not_bondy_enough(hbond_half_opt_pair &prm_hbond_pair ///< The hbond_half_opt_pair to examine (and potentially wipe)
		                                    ) {
			remove_not_bondy_enough( prm_hbond_pair.first  );
			remove_not_bondy_enough( prm_hbond_pair.second );
		}

		/// \brief Represent a possibly-bifurcating set of bonds from a given residue
		class bifur_hbond final {
		private:
			/// \brief The pair of best h-bonds from the NH atoms of this residue
			hbond_half_opt_pair for_this_nh;

			/// \brief The pair of best h-bonds from the CO atoms of this residue
			hbond_half_opt_pair for_this_co;

		public:
			bifur_hbond & update_for_this_nh(const hbond_half &);
			bifur_hbond & update_for_this_co(const hbond_half &);

			bifur_hbond & remove_not_bondy_enough();

			const hbond_half_opt_pair & get_bound_pair_for_this_nh() const;
			const hbond_half_opt_pair & get_bound_pair_for_this_co() const;
		};

		/// \brief Update the best h-bonds from the NH atoms of this residue with the new hbond_half
		inline bifur_hbond & bifur_hbond::update_for_this_nh(const hbond_half &prm_hbond ///< The new hbond_half with which to update
		                                                     ) {
			update_half_bond_pair( for_this_nh, prm_hbond );
			return *this;
		}

		/// \brief Update the best h-bonds from the CO atoms of this residue with the new hbond_half
		inline bifur_hbond & bifur_hbond::update_for_this_co(const hbond_half &prm_hbond ///< The new hbond_half with which to update
		                                                     ) {
			update_half_bond_pair( for_this_co, prm_hbond );
			return *this;
		}

		/// \brief Remove any parts of the bifur_hbond that aren't strong enough to be considered a true hbond
		inline bifur_hbond & bifur_hbond::remove_not_bondy_enough() {
			cath::sec::remove_not_bondy_enough( for_this_nh );
			cath::sec::remove_not_bondy_enough( for_this_co );
			return *this;
		}

		/// \brief Get the best seen pair of bonds from the NH atoms of this residue
		inline const hbond_half_opt_pair & bifur_hbond::get_bound_pair_for_this_nh() const {
			return for_this_nh;
		}

		/// \brief Get the best seen pair of bonds from the CO atoms of this residue
		inline const hbond_half_opt_pair & bifur_hbond::get_bound_pair_for_this_co() const {
			return for_this_co;
		}

		/// \brief Type alias for a vector of bifur_hbonds
		using bifur_hbond_vec = std::vector<bifur_hbond>;

		std::string to_string(const bifur_hbond &);

		std::ostream & operator<<(std::ostream &,
		                          const bifur_hbond &);

		/// \brief Represent the bifur_hbonds for the residues in a structure
		class bifur_hbond_list final {
		private:
			/// \brief The vector of bifur_hbonds
			bifur_hbond_vec bifur_hbonds;

		public:
			/// \brief A type alias for the const_iterator type to iterate over the bifur_hbonds
			using const_iterator = bifur_hbond_vec::const_iterator;

			explicit bifur_hbond_list(const size_t &);

			bool empty() const;
			size_t size() const;

			const bifur_hbond & operator[](const size_t &) const;

			bifur_hbond_list & remove_not_bondy_enough();

			bifur_hbond_list & update_with_nh_idx_co_idx_energy(const hbond_partner_t &,
			                                                    const hbond_partner_t &,
			                                                    const hbond_energy_t &);

			const_iterator begin() const;
			const_iterator end() const;
		};

		/// \brief Ctor
		inline bifur_hbond_list::bifur_hbond_list(const size_t &prm_size ///< The number of residues to represent 
		                                          ) : bifur_hbonds( prm_size ) {
		}

		/// \brief Whether this bifur_hbond_list is empty
		inline bool bifur_hbond_list::empty() const {
			return bifur_hbonds.empty();
		}

		/// \brief The number of residues this bifur_hbond_list represents
		inline size_t bifur_hbond_list::size() const {
			return bifur_hbonds.size();
		}

		/// \brief Return (a const-reference to) the bifur_hbond corresponding to the residue at the specified index
		inline const bifur_hbond & bifur_hbond_list::operator[](const size_t &prm_index ///< The index of the residue for which the bifur_hbond should be returned
		                                                        ) const {
			return bifur_hbonds[ prm_index ];
		}

		/// \brief Remove any candidate bonds that aren't strong enough to be considered true hbonds
		inline bifur_hbond_list & bifur_hbond_list::remove_not_bondy_enough() {
			for (bifur_hbond &the_bifur_hbond : bifur_hbonds) {
				the_bifur_hbond.remove_not_bondy_enough();
			}
			return *this;
		}

		/// \brief Update the bifur_hbond_list with the specified h-bond energy between
		///        the residues of specified indices
		inline bifur_hbond_list & bifur_hbond_list::update_with_nh_idx_co_idx_energy(const hbond_partner_t &prm_nh_index, ///< The index of the residue at the NH side of the h-bond
		                                                                             const hbond_partner_t &prm_co_index, ///< The index of the residue at the CO side of the h-bond
		                                                                             const hbond_energy_t  &prm_energy    ///< The energy of the h-bond
		                                                                             ) {
			bifur_hbonds[ prm_nh_index ].update_for_this_nh( { prm_co_index, prm_energy } );
			bifur_hbonds[ prm_co_index ].update_for_this_co( { prm_nh_index, prm_energy } );
			return *this;
		}

		/// \brief Standard begin() method to allow iteration over the bifur_hbonds
		inline auto bifur_hbond_list::begin() const -> const_iterator {
			return common::cbegin( bifur_hbonds );
		}

		/// \brief Standard end() method to allow iteration over the bifur_hbonds
		inline auto bifur_hbond_list::end() const -> const_iterator {
			return common::cend  ( bifur_hbonds );
		}

		/// \brief Remove any candidate bonds that aren't strong enough to be considered true hbonds
		///        in a copy of the specified bifur_hbond_list and then return the copy
		inline bifur_hbond_list remove_not_bondy_enough_copy(bifur_hbond_list prm_bifur_hbond_list ///< The bifur_hbond_list from which a copy should be taken, stripped and returned
		                                                     ) {
			prm_bifur_hbond_list.remove_not_bondy_enough();
			return prm_bifur_hbond_list;
		}

		std::string to_string(const bifur_hbond_list &);

		std::ostream & operator<<(std::ostream &,
		                          const bifur_hbond_list &);

	} // namespace sec
} // namespace cath

#endif
