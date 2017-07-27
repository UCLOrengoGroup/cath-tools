/// \map
/// \brief The overlap_frac_distn class header

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

#ifndef _CATH_TOOLS_SOURCE_CLUSTER_MAP_OVERLAP_FRAC_DISTN_H
#define _CATH_TOOLS_SOURCE_CLUSTER_MAP_OVERLAP_FRAC_DISTN_H

#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/operators.hpp>
#include <boost/range/combine.hpp>

#include "common/cpp14/cbegin_cend.hpp"
#include "common/size_t_literal.hpp"
#include "common/type_aliases.hpp"
#include "exception/invalid_argument_exception.hpp"

using namespace cath::common::literals;

namespace cath {
	namespace clust {

		namespace detail {

			/// \brief Compute ten to the specified power in a constexpr way
			constexpr size_t power_of_ten(const size_t &arg_expt,        ///< The exponent to which ten should be raised
			                              const size_t &arg_result = 1_z ///< Implementation value used in recursive calls
			                              ) {
				return (
					( arg_expt == 0 )
					? arg_result
					: power_of_ten( arg_expt - 1, 10_z * arg_result )
				);
			}
		}

		/// \brief Store the distribution of overlap fractions using a bunch of bin counts
		class overlap_frac_distn final : private boost::addable<overlap_frac_distn> {
		private:
			/// \brief The number of decimal places at which fractions should be distinguished
			static constexpr size_t num_dec_places = 4;

			/// \brief The number of gaps between the bins
			///
			/// This must be 10 ^ num_dec_places
			static constexpr size_t num_gaps = detail::power_of_ten( num_dec_places );

			static_assert( num_gaps == 10000, "" );
			
			/// \brief The total number of bins
			///
			/// This must be 1 + 10 ^ num_dec_places
			static constexpr size_t num_posns = 1_z + num_gaps;

			static_assert( num_posns == 10001, "" );

			/// \brief The bins: counts of the number seen of each fraction between 0 and 1 to num_dec_places decimal places
			size_vec fraction_counts{ size_vec( num_gaps + 1, 0 ) };

			/// \brief The number of overlap fractions stored so farm
			size_t num_counts = 0;

			size_t find_index_of_nth(const size_t &) const;

			static constexpr double get_fraction_of_index(const size_t &);
			static constexpr double get_percentage_of_index(const size_t &);
			static constexpr double get_index_of_fraction(const double &);
			static constexpr double get_index_of_percentage(const double &);

		public:

			/// \brief What to do about zero overlap factions when calculating statistics
			enum class zeroes_policy {
				INCLUDE, ///< Include the zeroes as normal
				EXCLUDE  ///< Exclude the zeroes from calculations (ie calculate as if they weren't there)
			};

			/// \brief A const_iterator type alias as part of making this a range over the bins
			///
			/// NOTE: be careful because the range is not over the values, it's over a bunch of bin-counts
			using const_iterator = size_vec::const_iterator;

			/// \brief An iterator type alias that just duplicates const_iterator to appease some Boost code
			///
			/// NOTE: be careful because the range is not over the values, it's over a bunch of bin-counts
			using iterator = const_iterator;

			overlap_frac_distn() = default;

			size_t size() const;
			bool empty() const;
			const_iterator begin() const;
			const_iterator end() const;

			overlap_frac_distn & add_overlap_fraction(const double &,
			                                          const size_t & = 1);
			overlap_frac_distn & operator+=(const overlap_frac_distn &);
			size_t get_num_in_range(const double &,
			                        const double &) const;
			double get_frac_at_percentile(const double &,
			                              const zeroes_policy & = zeroes_policy::INCLUDE) const;
		};

		overlap_frac_distn build_overlap_frac_distn_from_overlap_fractions(const doub_vec &);
		doub_doub_pair_vec percentile_data(const overlap_frac_distn &,
		                                   const doub_vec &,
		                                   const overlap_frac_distn::zeroes_policy & = overlap_frac_distn::zeroes_policy::INCLUDE);
		std::string percentile_markdown_table(const overlap_frac_distn &,
		                                      const doub_vec &,
		                                      const std::string &,
		                                      const std::string &,
		                                      const overlap_frac_distn::zeroes_policy & = overlap_frac_distn::zeroes_policy::INCLUDE);

		/// \brief Get the overlap fraction corresponding to the bin of the specified index
		inline constexpr double overlap_frac_distn::get_fraction_of_index(const size_t &arg_index ///< The index of the bin of interest
		                                                                  ) {
			return static_cast<double>( arg_index ) / static_cast<double>( num_gaps );
		}

		/// \brief Get the overlap percentage corresponding to the bin of the specified index
		inline constexpr double overlap_frac_distn::get_percentage_of_index(const size_t &arg_index ///< The index of the bin of interest
		                                                                    ) {
			return 100.0 * get_fraction_of_index( arg_index );
		}

		/// \brief Get the bin index associated with the specified overlap fraction
		inline constexpr double overlap_frac_distn::get_index_of_fraction(const double &arg_fraction ///< The overlap fraction to look for
		                                                                  ) {
			return static_cast<double>( arg_fraction ) * static_cast<double>( num_gaps );
		}

		/// \brief Get the bin index associated with the specified overlap percentage
		inline constexpr double overlap_frac_distn::get_index_of_percentage(const double &arg_percentage ///< The overlap percentage to look for
		                                                                    ) {
			return get_index_of_fraction( arg_percentage / 100.0 );
		}

		/// \brief The number of overlap fractions
		inline size_t overlap_frac_distn::size() const {
			return num_counts;
		}

		/// \brief Whether no overlap fractions have been stored
		inline bool overlap_frac_distn::empty() const {
			return ( size() == 0 );
		}

		/// \brief Standard const begin() method, as part of making this into a range over the bins
		///
		/// NOTE: be careful because the range is not over the values, it's over a bunch of bin-counts
		inline auto overlap_frac_distn::begin() const -> const_iterator {
			return common::cbegin( fraction_counts );
		}

		/// \brief Standard const begin() method, as part of making this into a range over the bins
		///
		/// NOTE: be careful because the range is not over the values, it's over a bunch of bin-counts
		inline auto overlap_frac_distn::end() const -> const_iterator {
			return common::cend  ( fraction_counts );
		}

		/// \brief Add an overlap fraction to the overlap_frac_distn
		inline overlap_frac_distn & overlap_frac_distn::add_overlap_fraction(const double &arg_overlap_fraction, ///< The overlap fraction to add
		                                                                     const size_t &arg_count             ///< The number of instances of this fraction to add (default 1)
		                                                                     ) {
			using std::to_string;
			if ( ! ::boost::math::isfinite( arg_overlap_fraction ) || arg_overlap_fraction < 0.0 || arg_overlap_fraction > 1.0 ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception(
					"Unable to add invalid overlap fraction "
					+ to_string( arg_overlap_fraction )
					+ " to overlap_frac_distn"
				));
			}

			const auto index = static_cast<size_t>( round( static_cast<double>( num_gaps ) * arg_overlap_fraction ) );
			( fraction_counts[ index ] ) += arg_count;
			num_counts += arg_count;

			return *this;
		}

		/// \brief Add all the overlap fractions from another overlap_frac_distn into this one
		inline overlap_frac_distn & overlap_frac_distn::operator+=(const overlap_frac_distn &arg_distn ///< The overlap_frac_distn to add
		                                                           ) {
			for (boost::tuple<size_t &, const size_t &> &&x : boost::range::combine( fraction_counts, arg_distn ) ) {
				x.get<0>() += x.get<1>();
			}
			num_counts += arg_distn.size();
			return *this;
		}

	} // namespace clust
} // namespace cath

#endif
