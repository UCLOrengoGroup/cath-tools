/// \file
/// \brief The view_cache_index_entry test suite

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

#include <boost/numeric/conversion/cast.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/numeric.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "alignment/pair_alignment.h"
#include "common/algorithm/generate_n_build.h"
#include "common/algorithm/transform_build.h"
#include "common/chrono/chrono_type_aliases.h"
#include "common/random/pick_random_pair.h"
#include "common/size_t_literal.h"
#include "exception/invalid_argument_exception.h"
#include "structure/protein/protein_source_file_set/protein_source_from_pdb.h"
#include "structure/protein/sec_struc.h"
#include "structure/protein/sec_struc_planar_angles.h"
#include "structure/view_cache/index/detail/vcie_match_criteria.h"
#include "structure/view_cache/index/view_cache_index_entry.h"
#include "test/global_test_constants.h"

#include <chrono>
#include <random>

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::index;
using namespace cath::index::detail;
using namespace std;
using namespace std::chrono;

using boost::accumulate;
using boost::ignore_unused;
using boost::irange;
using boost::numeric_cast;
using boost::range::combine;

namespace cath {
	namespace test {

		/// \brief The view_cache_index_entry_test_suite_fixture to assist in testing
		struct view_cache_index_entry_test_suite_fixture : protected global_test_constants {
		private:
			static float_score_type compare(const view_cache_index_entry_vec &,
			                                const view_cache_index_entry_vec &,
											const vcie_match_criteria &);
			static hrc_duration time_comparison(const vcie_vcie_vec_pair &);
			static view_cache_index_entry_vec build_random_vcies(const size_t  &,
			                                                     const protein &,
			                                                     mt19937 &);
		protected:
			~view_cache_index_entry_test_suite_fixture() noexcept = default;

			static vcie_vcie_vec_pair build_random_vcies_pair(const size_t &,
			                                                  const protein &,
			                                                  const protein &,
			                                                  mt19937 &);

			static vcie_vcie_vec_pair build_alignment_vcies_pair(const size_t &,
			                                                     const protein &,
			                                                     const protein &,
			                                                     const alignment &,
			                                                     mt19937 &);

			hrc_duration_vec time_comparisons(const vcie_vcie_vec_pair &,
		                                      const size_t &);

			static constexpr size_t NUM_ENTRIES = 1000000;
//			static constexpr size_t NUM_ENTRIES =  500000;
			static constexpr size_t NUM_REPEATS =      50;
		};

		/// \brief TODOCUMENT
		float_score_type view_cache_index_entry_test_suite_fixture::compare(const view_cache_index_entry_vec &arg_vcies_a,  ///< TODOCUMENT
		                                                                    const view_cache_index_entry_vec &arg_vcies_b,  ///< TODOCUMENT
								                                            const vcie_match_criteria        &arg_criteria ///< TODOCUMENT
								                                            ) {
			float_score_type total_score = 0.0;
			auto vcies_a_itr = cath::common::cbegin( arg_vcies_a );
			auto vcies_b_itr = cath::common::cbegin( arg_vcies_b );
			const auto vcies_end_itr = cath::common::cend( arg_vcies_a );
			while ( vcies_a_itr != vcies_end_itr ) {
				if ( arg_criteria( *vcies_a_itr, *vcies_b_itr ) ) {
					const float_score_type distance = sqrt( squared_distance( *vcies_a_itr, *vcies_b_itr ) );
					const float_score_type score    = static_cast<float_score_type>( 1.0 ) - (distance / static_cast<float_score_type>( 7.0 ) );
					if ( score > 0 ) {
						total_score += score;
					}
				}
				++vcies_a_itr;
				++vcies_b_itr;
			}
			return total_score;

//			return accumulate(
//				combine( arg_vcies_a, arg_vcies_b ),
//				static_cast<float_score_type>( 0.0 ),
//				[&] (const float_score_type &x, const boost::tuple<const view_cache_index_entry &, const view_cache_index_entry &> &y) {
//					if ( ! arg_criteria( y.get<0>(), y.get<1>() ) ) {
//						return x;
//					}
//					const float_score_type distance = sqrt( squared_distance( y.get<0>(), y.get<1>() ) );
//					return x + std::max(
//						static_cast<float_score_type>( 0.0 ),
//						static_cast<float_score_type>( 1.0 ) - (distance / static_cast<float_score_type>( 7.0 ) )
//					);
//				}
//			);
		}

		/// \brief TODOCUMENT
		hrc_duration view_cache_index_entry_test_suite_fixture::time_comparison(const vcie_vcie_vec_pair &arg_vcies_pair ///< TODOCUMENT
		                                                                        ) {
			const auto &vcies_1 = arg_vcies_pair.first;
			const auto &vcies_2 = arg_vcies_pair.second;

			if ( vcies_1.size() != vcies_2.size() ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception(""));
			}

			const hrc_time_point before_time_point = high_resolution_clock::now();
			compare( vcies_1, vcies_2, make_default_vcie_match_criteria() );
			return high_resolution_clock::now() - before_time_point;
		}

		/// \brief TODOCUMENT
		view_cache_index_entry_vec view_cache_index_entry_test_suite_fixture::build_random_vcies(const size_t  &arg_num_entries, ///< TODOCUMENT
		                                                                                         const protein &arg_protein,     ///< TODOCUMENT
		                                                                                         mt19937       &arg_rng          ///< TODOCUMENT
		                                                                                         ) {
			if ( arg_protein.get_length() < 2 ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot build random vcies if the protein has fewer than two residues"));
			}
			const auto &length = arg_protein.get_length();
			return generate_n_build<view_cache_index_entry_vec>(
				arg_num_entries,
				[&] () {
					const auto index_pair = pick_random_pair( 0_z, length - 1, arg_rng );
					return make_view_cache_index_entry(
						arg_protein,
						index_pair.first,
						index_pair.second
					);
				}
			);
		}

		/// \brief TODOCUMENT
		vcie_vcie_vec_pair view_cache_index_entry_test_suite_fixture::build_random_vcies_pair(const size_t  &arg_num_entries, ///< TODOCUMENT
		                                                                                      const protein &arg_protein_a,   ///< TODOCUMENT
		                                                                                      const protein &arg_protein_b,   ///< TODOCUMENT
		                                                                                      mt19937       &arg_rng          ///< TODOCUMENT
		                                                                                      ) {
			return make_pair(
				build_random_vcies( arg_num_entries, arg_protein_a, arg_rng ),
				build_random_vcies( arg_num_entries, arg_protein_b, arg_rng )
			);
		}

		/// \brief TODOCUMENT
		vcie_vcie_vec_pair view_cache_index_entry_test_suite_fixture::build_alignment_vcies_pair(const size_t    &arg_num_entries, ///< TODOCUMENT
		                                                                                         const protein   &arg_protein_a,   ///< TODOCUMENT
		                                                                                         const protein   &arg_protein_b,   ///< TODOCUMENT
		                                                                                         const alignment &arg_alignment,   ///< TODOCUMENT
		                                                                                         mt19937         &arg_rng          ///< TODOCUMENT
		                                                                                         ) {
			const auto &present_posn_indices     = indices_of_present_positions_of_both_entries( arg_alignment );
			const auto  num_present_posn_indices = present_posn_indices.size();
			if ( num_present_posn_indices < 2 ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot build random_vcie pairs if the alignment has fewer than two residues"));
			}
			vcie_vcie_vec_pair results;
			for (const auto &entry_ctr : irange( 0_z, arg_num_entries ) ) {
				ignore_unused( entry_ctr );
				const auto  index_pair = pick_random_pair( 0_z, num_present_posn_indices - 1, arg_rng );
				const auto &index_1    = present_posn_indices[ index_pair.first  ];
				const auto &index_2    = present_posn_indices[ index_pair.second ];
				results.first.push_back( make_view_cache_index_entry(
					arg_protein_a,
					get_a_position_of_index( arg_alignment, index_1 ),
					get_a_position_of_index( arg_alignment, index_2 )
				) );
				results.second.push_back( make_view_cache_index_entry(
					arg_protein_b,
					get_b_position_of_index( arg_alignment, index_1 ),
					get_b_position_of_index( arg_alignment, index_2 )
				) );
			}
			return results;
		}

		/// \brief TODOCUMENT
		hrc_duration_vec view_cache_index_entry_test_suite_fixture::time_comparisons(const vcie_vcie_vec_pair &arg_vcies_pair, ///< TODOCUMENT
		                                                                             const size_t             &arg_num_repeats ///< TODOCUMENT
		                                                                             ) {
			const auto vcies_1 = arg_vcies_pair.first;
			const auto vcies_2 = arg_vcies_pair.second;

			if ( vcies_1.size() != vcies_1.size() ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception(""));
			}
			return generate_n_build<hrc_duration_vec>(
				arg_num_repeats,
				[&] () { return time_comparison( arg_vcies_pair ); }
			);
		}
	}
}

constexpr size_t cath::test::view_cache_index_entry_test_suite_fixture::NUM_ENTRIES;
constexpr size_t cath::test::view_cache_index_entry_test_suite_fixture::NUM_REPEATS;

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(view_cache_index_entry_test_suite, cath::test::view_cache_index_entry_test_suite_fixture)

///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(basic) {
//	mt19937 rng( random_device{}() );
//	ostringstream parse_ss;
//	const protein protein_a      = read_protein_from_files( protein_source_from_pdb(), TEST_SOURCE_DATA_DIR(), "1n3lA01", parse_ss );
//	const protein protein_b      = read_protein_from_files( protein_source_from_pdb(), TEST_SOURCE_DATA_DIR(), "1r6xA02", parse_ss );
//	const auto random_vcies_pair = build_random_vcies_pair( NUM_ENTRIES, protein_a, protein_b, rng );
//	const auto durns             = time_comparisons( random_vcies_pair, NUM_REPEATS );
//	const auto rates             = transform_build<doub_vec>(
//		durns,
//		[] (const hrc_duration &x) {
//			const double durn_in_seconds = numeric_cast<double>( duration_cast<nanoseconds>( x ).count() ) / 1000000000.0;
//			return numeric_cast<double>( NUM_ENTRIES ) / durn_in_seconds;
//		}
//	);
//	const auto rate_sum          = accumulate( rates, 0.0, plus<double>() );
//	const auto mean_rate         = rate_sum / numeric_cast<double>( rates.size() );
//	for (const auto &x : rates) {
//		cerr << "Rate per second for " << NUM_ENTRIES << " is " << fixed << x << endl;
//	}
//	cerr << "Mean rate per second is " << mean_rate << endl;
//	BOOST_CHECK_EQUAL(0, 0);
//}

BOOST_AUTO_TEST_SUITE_END()
