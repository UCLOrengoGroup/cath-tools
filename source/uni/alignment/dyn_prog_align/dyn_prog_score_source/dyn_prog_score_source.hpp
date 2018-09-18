/// \file
/// \brief The dyn_prog_score_source class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_SCORE_SOURCE_DYN_PROG_SCORE_SOURCE_HPP
#define _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_SCORE_SOURCE_DYN_PROG_SCORE_SOURCE_HPP

#include "common/exception/invalid_argument_exception.hpp"
#include "common/temp_check_offset_1.hpp"
#include "common/type_aliases.hpp"

namespace cath {
	namespace align {

		/// \brief Provide ABC interface for classes that provide scores to be used for aligning with dynamic-programming
		///
		/// Note: Consider whether this dynamic-polymorphism approach may be scrapped in favour of
		///       static-polymorphism. The main motivation for switching would probably be if the run-time approach
		///       turned out to be too slow.
		class dyn_prog_score_source {
		private:
			/// \brief Return the number of elements in the first entry to by aligned with dynamic-programming
			virtual size_t do_get_length_a() const = 0;

			/// \brief Return the number of elements in the second entry to by aligned with dynamic-programming
			virtual size_t do_get_length_b() const = 0;

			/// \brief Return the score the number of elements in the first entry to by aligned with dynamic-programming
			virtual score_type do_get_score(const size_t &,
			                                const size_t &) const = 0;

		public:
			dyn_prog_score_source() = default;
			virtual ~dyn_prog_score_source() noexcept = default;

			dyn_prog_score_source(const dyn_prog_score_source &) = default;
			dyn_prog_score_source(dyn_prog_score_source &&) noexcept = default;
			dyn_prog_score_source & operator=(const dyn_prog_score_source &) = default;
			dyn_prog_score_source & operator=(dyn_prog_score_source &&) noexcept = default;

			size_t get_length_a() const;
			size_t get_length_b() const;
			score_type get_score(const size_t &,
			                     const size_t &) const;
		};

		score_type get_score__offset_1(const dyn_prog_score_source &,
		                               const size_t &,
		                               const size_t &);

		/// \brief An NVI pass-through method to get the score for the specified indices to be used for aligning with dynamic-programming
		inline score_type dyn_prog_score_source::get_score(const size_t &prm_index_a, ///< The index of the element of interest in the first sequence
		                                                   const size_t &prm_index_b  ///< The index of the element of interest in the second sequence
		                                                   ) const {
#ifndef NDEBUG
			// Check that the two indices are valid
			if ( prm_index_a >= get_length_a() ) {
				BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("First index is out of range when getting score for aligning with dynamic-programming"));
			}
			if ( prm_index_b >= get_length_b() ) {
				BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Second index is out of range when getting score for aligning with dynamic-programming"));
			}
#endif

			// Pass-through to the concrete do_get_score() to do the real work
			return do_get_score( prm_index_a, prm_index_b );
		}

		/// \brief A non-member, non-friend helper function that gets a score from dyn_prog_score_source objects with offset_1 indices
		inline score_type get_score__offset_1(const dyn_prog_score_source &prm_dyn_prog_score_source, ///< The dyn_prog_score_source from which to get the score
		                                      const size_t                &prm_index_a__offset_1,     ///< The index of the element of interest in the first  sequence (using offset 1)
		                                      const size_t                &prm_index_b__offset_1      ///< The index of the element of interest in the second sequence (using offset 1)
		                                      ) {
			// Check the offsets are valid
			check_offset_1( prm_index_a__offset_1 );
			check_offset_1( prm_index_b__offset_1 );

			// Call the member get_score() function with adjusted indices
			return prm_dyn_prog_score_source.get_score( prm_index_a__offset_1 - 1,
			                                            prm_index_b__offset_1 - 1 );
		}
	} // namespace align

} // namespace cath
#endif
