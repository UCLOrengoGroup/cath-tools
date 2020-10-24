/// \file
/// \brief The clust_mapping_spec class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_OPTIONS_SPEC_CLUST_MAPPING_SPEC_HPP
#define _CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_OPTIONS_SPEC_CLUST_MAPPING_SPEC_HPP

#include <boost/operators.hpp>

#include <stdexcept>

namespace cath {
	namespace clust {

		/// \brief Specify the details of how to clusters should be mapped
		///
		/// Note that this stores as absolute fractions, not as percentages
		///
		/// Possible extra option:
		///
		///      --renumber-new-start-num             For --map-from-member-file, renumber unmapped clusters starting from <num>
		///                                           (rather than from one above the highest number in the <file>))
		class clust_mapping_spec final : boost::equality_comparable<clust_mapping_spec> {
		private:
			/// \brief The fraction that the overlap over the longest of two domains must exceed for them to be considered equivalent
			double min_equiv_dom_ol   = DEFAULT_MIN_EQUIV_DOM_OL;

			/// \brief The fraction of the old cluster's entries that must map to a new cluster for them to be considered equivalent
			double min_equiv_clust_ol = DEFAULT_MIN_EQUIV_CLUST_OL;

			static constexpr double check_frac_against_strict_min_and_return(const double &,
			                                                                 const double &);

		public:
			/// \brief Default value for the fraction that the overlap over the longest of two domains must exceed for them to be considered equivalent
			static constexpr double DEFAULT_MIN_EQUIV_DOM_OL           = 0.6;

			/// \brief Default value for the fraction of the old cluster's entries that must map to a new cluster for them to be considered equivalent
			static constexpr double DEFAULT_MIN_EQUIV_CLUST_OL         = 0.6;

			/// \brief Strict minimum value for the fraction that the overlap over the longest of two domains must exceed for them to be considered equivalent
			static constexpr double MIN_MIN_EQUIV_DOM_OL               = 0.5;

			/// \brief Strict minimum value for the fraction of the old cluster's entries that must map to a new cluster for them to be considered equivalent
			static constexpr double MIN_MIN_EQUIV_CLUST_OL             = 0.5;

			/// \brief Strict minimum value for the fraction of the new cluster's entries that must map to a map-from cluster for them to be considered equivalent
			static constexpr double MIN_EQUIV_FRAC_OF_NEW_CLUST        = 0.2;

			/// \brief Strict minimum value for the fraction of the new cluster's entries that have mapped anywhere that must map to a map-from cluster for them to be considered equivalent
			static constexpr double MIN_EQUIV_FRAC_OF_NEW_CLUST_EQUIVS = 0.5;

			/// \brief Default ctor
			constexpr clust_mapping_spec() = default;

			explicit constexpr clust_mapping_spec(const double &,
			                                      const double & = MIN_MIN_EQUIV_CLUST_OL);

			constexpr const double & get_min_equiv_dom_ol() const;
			constexpr const double & get_min_equiv_clust_ol() const;

			clust_mapping_spec & set_min_equiv_dom_ol(const double &);
			clust_mapping_spec & set_min_equiv_clust_ol(const double &);
		};

		/// \brief Check a fraction is within [0, 1] and is strictly greater than the specified minimum.
		///        Throw if not or return the value otherwise
		inline constexpr double clust_mapping_spec::check_frac_against_strict_min_and_return(const double &prm_frac, ///< The fraction to check
		                                                                                     const double &prm_min   ///< The strict minimum to compare the fraction against
		                                                                                     ) {
			return ( prm_frac >= 0.0 && prm_frac <= 1.0 && prm_frac >= prm_min )
				? prm_frac
				: throw std::invalid_argument("The specified cluster mapping fraction is out of range");
		}

		/// \brief Ctor from the domain and cluster overlap fractions
		inline constexpr clust_mapping_spec::clust_mapping_spec(const double &prm_min_equiv_dom_ol,  ///< The fraction that the overlap over the longest of two domains must exceed for them to be considered equivalent
		                                                        const double &prm_min_equiv_clust_ol ///< The fraction of the old cluster's entries that must map to a new cluster for them to be considered equivalent
		                                                        ) : min_equiv_dom_ol   { check_frac_against_strict_min_and_return( prm_min_equiv_dom_ol,   MIN_MIN_EQUIV_DOM_OL   ) },
		                                                            min_equiv_clust_ol { check_frac_against_strict_min_and_return( prm_min_equiv_clust_ol, MIN_MIN_EQUIV_CLUST_OL ) } {
		}

		/// \brief Getter for the fraction that the overlap over the longest of two domains must exceed for them to be considered equivalent
		inline constexpr const double & clust_mapping_spec::get_min_equiv_dom_ol() const {
			return min_equiv_dom_ol;
		}

		/// \brief Getter for the fraction of the old cluster's entries that must map to a new cluster for them to be considered equivalent
		inline constexpr const double & clust_mapping_spec::get_min_equiv_clust_ol() const {
			return min_equiv_clust_ol;
		}

		/// \brief Setter for the fraction that the overlap over the longest of two domains must exceed for them to be considered equivalent
		inline clust_mapping_spec & clust_mapping_spec::set_min_equiv_dom_ol(const double &prm_min_equiv_dom_ol ///< The fraction that the overlap over the longest of two domains must exceed for them to be considered equivalent
		                                                                     ) {
			min_equiv_dom_ol = check_frac_against_strict_min_and_return( prm_min_equiv_dom_ol, MIN_MIN_EQUIV_DOM_OL );
			return *this;
		}

		/// \brief Setter for the fraction of the old cluster's entries that must map to a new cluster for them to be considered equivalent
		inline clust_mapping_spec & clust_mapping_spec::set_min_equiv_clust_ol(const double &prm_min_equiv_clust_ol ///< The fraction of the old cluster's entries that must map to a new cluster for them to be considered equivalent
		                                                                       ) {
			min_equiv_clust_ol = check_frac_against_strict_min_and_return( prm_min_equiv_clust_ol, MIN_MIN_EQUIV_CLUST_OL );
			return *this;
		}

		/// \brief Return whether the two specified clust_mapping_specs are identical
		///
		/// \relates clust_mapping_spec
		inline constexpr bool operator==(const clust_mapping_spec &prm_lhs, ///< The first  clust_mapping_spec to compare
		                                 const clust_mapping_spec &prm_rhs  ///< The second clust_mapping_spec to compare
		                                 ) {
			return (
				prm_lhs.get_min_equiv_dom_ol()   == prm_rhs.get_min_equiv_dom_ol()
				&&
				prm_lhs.get_min_equiv_clust_ol() == prm_rhs.get_min_equiv_clust_ol()
			);
		}

	} // namespace clust
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_OPTIONS_SPEC_CLUST_MAPPING_SPEC_HPP
