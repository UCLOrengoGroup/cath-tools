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

#ifndef _CATH_TOOLS_SOURCE_CLUSTER_OPTIONS_SPEC_CLUST_MAPPING_SPEC_H
#define _CATH_TOOLS_SOURCE_CLUSTER_OPTIONS_SPEC_CLUST_MAPPING_SPEC_H

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
			/// \brief The minimum fraction overlap for domains TODOCUMENT FURTHER
			double min_equiv_dom_ol   = DEFAULT_MIN_EQUIV_DOM_OL;

			/// \brief The minimum fraction overlap for clusters TODOCUMENT FURTHER
			double min_equiv_clust_ol = DEFAULT_MIN_EQUIV_CLUST_OL;

			static constexpr double check_frac_against_strict_min_and_return(const double &,
			                                                                 const double &);

		public:
			/// \brief Default value for the minimum fraction overlap for domains TODOCUMENT FURTHER
			static constexpr double DEFAULT_MIN_EQUIV_DOM_OL   = 0.6;

			/// \brief Default value for the minimum fraction overlap for clusters TODOCUMENT FURTHER
			static constexpr double DEFAULT_MIN_EQUIV_CLUST_OL = 0.6;

			/// \brief Strict minimum value for the minimum fraction overlap for domains TODOCUMENT FURTHER
			static constexpr double MIN_MIN_EQUIV_DOM_OL       = 0.5;

			/// \brief Strict minimum value for the minimum fraction overlap for clusters TODOCUMENT FURTHER
			static constexpr double MIN_MIN_EQUIV_CLUST_OL     = 0.5;

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
		inline constexpr double clust_mapping_spec::check_frac_against_strict_min_and_return(const double &arg_frac, ///< The fraction to check
		                                                                                     const double &arg_min   ///< The strict minimum to compare the fraction against
		                                                                                     ) {
			return ( arg_frac >= 0.0 && arg_frac <= 1.0 && arg_frac > arg_min )
				? arg_frac
				: throw std::invalid_argument("The specified cluster mapping fraction is out of range");
		}

		/// \brief Ctor from the domain and cluster overlap fractions
		inline constexpr clust_mapping_spec::clust_mapping_spec(const double &arg_min_equiv_dom_ol,  ///< The minimum fraction overlap for domains TODOCUMENT FURTHER
		                                                        const double &arg_min_equiv_clust_ol ///< The minimum fraction overlap for clusters TODOCUMENT FURTHER
		                                                        ) : min_equiv_dom_ol   { check_frac_against_strict_min_and_return( arg_min_equiv_dom_ol,   MIN_MIN_EQUIV_DOM_OL   ) },
		                                                            min_equiv_clust_ol { check_frac_against_strict_min_and_return( arg_min_equiv_clust_ol, MIN_MIN_EQUIV_CLUST_OL ) } {
		}

		/// \brief Getter for the minimum fraction overlap for domains TODOCUMENT FURTHER
		inline constexpr const double & clust_mapping_spec::get_min_equiv_dom_ol() const {
			return min_equiv_dom_ol;
		}

		/// \brief Getter for the minimum fraction overlap for clusters TODOCUMENT FURTHER
		inline constexpr const double & clust_mapping_spec::get_min_equiv_clust_ol() const {
			return min_equiv_clust_ol;
		}

		/// \brief Setter for the minimum fraction overlap for domains TODOCUMENT FURTHER
		inline clust_mapping_spec & clust_mapping_spec::set_min_equiv_dom_ol(const double &arg_min_equiv_dom_ol ///< The minimum fraction overlap for domains TODOCUMENT FURTHER
		                                                                     ) {
			min_equiv_dom_ol = check_frac_against_strict_min_and_return( arg_min_equiv_dom_ol, MIN_MIN_EQUIV_DOM_OL );
			return *this;
		}

		/// \brief Setter for the minimum fraction overlap for clusters TODOCUMENT FURTHER
		inline clust_mapping_spec & clust_mapping_spec::set_min_equiv_clust_ol(const double &arg_min_equiv_clust_ol ///< The minimum fraction overlap for clusters TODOCUMENT FURTHER
		                                                                       ) {
			min_equiv_clust_ol = check_frac_against_strict_min_and_return( arg_min_equiv_clust_ol, MIN_MIN_EQUIV_CLUST_OL );
			return *this;
		}

		/// \brief Return whether the two specified clust_mapping_specs are identical
		///
		/// \relates clust_mapping_spec
		inline constexpr bool operator==(const clust_mapping_spec &arg_lhs, ///< The first  clust_mapping_spec to compare
		                                 const clust_mapping_spec &arg_rhs  ///< The second clust_mapping_spec to compare
		                                 ) {
			return (
				arg_lhs.get_min_equiv_dom_ol()   == arg_rhs.get_min_equiv_dom_ol()
				&&
				arg_lhs.get_min_equiv_clust_ol() == arg_rhs.get_min_equiv_clust_ol()
			);
		}

	} // namespace clust
} // namespace cath

#endif
