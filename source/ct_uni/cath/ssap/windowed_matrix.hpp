/// \file
/// \brief The windowed_matrix class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SSAP_WINDOWED_MATRIX_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SSAP_WINDOWED_MATRIX_HPP

#include <boost/algorithm/clamp.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "cath/common/config.hpp"
#include "cath/common/debug_numeric_cast.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/common/temp_check_offset_1.hpp"
#include "cath/common/type_aliases.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace cath {

	using cath::common::literals::operator""_z;

	/// \brief A memory-efficient 2d matrix that's restricted to a window around the leading diagonal
	///
	/// Here's an example of a 7 x 9 matrix with a requested window size of 5:
	///
	///     *  *  *  +  .  .  .  .  .
	///     +  *  *  *  +  .  .  .  .
	///     .  +  *  *  *  +  .  .  .
	///     .  .  +  *  *  *  +  .  .
	///     .  .  .  +  *  *  *  +  .
	///     .  .  .  .  +  *  *  *  +
	///     .  .  .  .  .  +  *  *  *
	///
	/// (a full stop indicates a cell that's inactive because it's outside the window)
	///
	/// The asterisks indicate a central window of size 3 (within the window of size 5)
	/// that's defined by the two diagonals through the top-left and bottom right.
	/// This central window's size is always the difference in lengths plus one.
	/// Any requested window size must be at least this large.
	///
	/// Where the requested window size minus the central window size is an even number
	/// (as in the example above) then the excess can be evenly split between the two sides of
	/// the central window.
	///
	/// Where the requested window size minus the central window size is an odd number,
	/// the excess cannot be split evenly so the spare one is always placed on the side nearest
	/// the start of the longer edge, eg:
	///
	///     *  *  .  .  .
	///     +  *  *  .  .
	///     .  +  *  *  .
	///     .  .  +  *  *
	///
	/// Where the size of the excess is odd and the two lengths are even, there is no way to
	/// decide which side should get the spare one. In this case, the code just adds one to
	/// both sides, thus creating a window one larger than was requested.
	///
	/// This way the code can guarantee that swapping the dimensions of the matrix will always
	/// just transpose the same matrix.
	///
	/// Note that the \c a and \c b used in the names throughout this code are arbitrary;
	/// they happen to correspond to the \c a and \c b used throughout the SSAP routines
	/// but they needn't necessarily.
	class windowed_matrix final {
	public:
		using size_type = doub_vec::size_type;
		using difference_type = doub_vec::difference_type;

	private:
		/// \brief The length of the first entry (ie the number of rows in the matrix)
		size_type length_a;

		/// \brief The length of the second entry (ie the number of columns in the matrix)
		size_type length_b;

		/// \brief The requested width (or equivalently, height) of the window
		///
		/// Note the actual window size may be one larger than requested if the lengths
		/// are equal and the requested window size is even. See windowed_matrix notes
		/// for more information.
		size_type requested_window_size;

		/// \brief The container that is used to store the data.
		///
		/// This is a single dimensional vector; get_internal_index() is used
		/// to translate the two dimensional coordinates to access this vector.
		///
		/// This allows the windowed to be memory efficient.
		doub_vec data;

		[[nodiscard]] size_t get_internal_index( const size_type &, const size_type & ) const;

	  public:
		windowed_matrix(const size_type &,
		                const size_type &,
		                const size_type &);

		[[nodiscard]] size_type get_length_a() const;
		[[nodiscard]] size_type get_length_b() const;
		[[nodiscard]] size_type get_window_size() const;

		void set_value(const size_type &,
		               const size_type &,
		               const double &);
		[[nodiscard]] double get_value( const size_type &, const size_type & ) const;
	};

	void check_indices_are_within_window(const windowed_matrix::size_type &,
	                                     const windowed_matrix::size_type &,
	                                     const windowed_matrix::size_type &,
	                                     const windowed_matrix::size_type &,
	                                     const windowed_matrix::size_type &);

	void check_lengths_and_window_size_are_valid(const windowed_matrix::size_type &,
	                                             const windowed_matrix::size_type &,
	                                             const windowed_matrix::size_type &);

	size_size_pair get_window_upper_and_lower_part_widths(const windowed_matrix::size_type &,
	                                                      const windowed_matrix::size_type &,
	                                                      const windowed_matrix::size_type &);

	int get_window_indexing_offset_for_b__offset_1(const windowed_matrix::size_type &,
	                                               const windowed_matrix::size_type &,
	                                               const windowed_matrix::size_type &,
	                                               const windowed_matrix::size_type &);

	int get_window_matrix_a_index(const windowed_matrix::size_type &,
	                              const windowed_matrix::size_type &,
	                              const windowed_matrix::size_type &,
	                              const windowed_matrix::size_type &,
	                              const windowed_matrix::size_type &);

	int get_window_matrix_a_index__offset_1(const windowed_matrix::size_type &,
	                                        const windowed_matrix::size_type &,
	                                        const windowed_matrix::size_type &,
	                                        const windowed_matrix::size_type &,
	                                        const windowed_matrix::size_type &);

	size_t get_window_start_a_for_b__offset_1(const windowed_matrix::size_type &,
	                                          const windowed_matrix::size_type &,
	                                          const windowed_matrix::size_type &,
	                                          const windowed_matrix::size_type &);

	size_t get_window_stop_a_for_b__offset_1(const windowed_matrix::size_type &,
	                                         const windowed_matrix::size_type &,
	                                         const windowed_matrix::size_type &,
	                                         const windowed_matrix::size_type &);

	size_t get_window_width_for_full_matrix(const windowed_matrix::size_type &,
	                                        const windowed_matrix::size_type &);


	/// \brief TODOCUMENT
	///
	/// \relates windowed_matrix
	inline void check_indices_are_within_window(const windowed_matrix::size_type &prm_length_a,         ///< The length of the first  entry (or the number of rows    in     matrix)
	                                            const windowed_matrix::size_type &prm_length_b,         ///< The length of the second entry (or the number of columns in the matrix)
	                                            const windowed_matrix::size_type &prm_requested_window, ///< The requested width (or equivalently, height) of the window
	                                            const windowed_matrix::size_type &prm_index_a,          ///< TODOCUMENT
	                                            const windowed_matrix::size_type &prm_index_b           ///< TODOCUMENT
	                                            ) {
		if constexpr ( common::IS_IN_DEBUG_MODE ) {
			check_lengths_and_window_size_are_valid(
				prm_length_a,
				prm_length_b,
				prm_requested_window
			);
			const size_t window_start_a = get_window_start_a_for_b__offset_1(
				prm_length_a,
				prm_length_b,
				prm_requested_window,
				prm_index_b + 1
			) - 1;
			const size_t window_stop_a  = get_window_stop_a_for_b__offset_1(
				prm_length_a,
				prm_length_b,
				prm_requested_window,
				prm_index_b + 1
			) - 1;
			if (prm_index_a < window_start_a || prm_index_a > window_stop_a) {
				BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception(
					"Indices "
					+ boost::lexical_cast<std::string>(prm_index_a)
					+ ", "
					+ boost::lexical_cast<std::string>(prm_index_b)
					+ " are not valid in a matrix of size "
					+ boost::lexical_cast<std::string>(prm_length_a)
					+ ", "
					+ boost::lexical_cast<std::string>(prm_length_b)
					+ " with a window width of "
					+ boost::lexical_cast<std::string>(prm_requested_window)
					+ " (window for a at b-index of "
					+ boost::lexical_cast<std::string>(prm_index_b)
					+ " is [ "
					+ boost::lexical_cast<std::string>(window_start_a)
					+ ", "
					+ boost::lexical_cast<std::string>(window_stop_a)
					+ "])"
				));
			}
		}
	}

	/// \brief TODOCUMENT
	///
	/// \relates windowed_matrix
	inline void check_lengths_and_window_size_are_valid(const windowed_matrix::size_type &prm_length_a,        ///< The length of the first  entry (or the number of rows    in     matrix)
	                                                    const windowed_matrix::size_type &prm_length_b,        ///< The length of the second entry (or the number of columns in the matrix)
	                                                    const windowed_matrix::size_type &prm_requested_window ///< The requested width (or equivalently, height) of the window
	                                                    ) {
		// Calculate the length difference and check that the window size is large enough to handle the length difference
		const auto diff_t_length_a   = debug_numeric_cast<windowed_matrix::difference_type>( prm_length_a );
		const auto diff_t_length_b   = debug_numeric_cast<windowed_matrix::difference_type>( prm_length_b );
		const windowed_matrix::difference_type length_difference = diff_t_length_a - diff_t_length_b;
		if ( prm_requested_window < 1 + debug_numeric_cast<windowed_matrix::size_type>( labs( length_difference ) ) ) {
			BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception(
				"The window size is "
				+ boost::lexical_cast<std::string>( prm_requested_window )
				+ " but it should be at least as large as the magnitude of the difference in lengths (the lengths are "
				+ boost::lexical_cast<std::string>( prm_length_a )
				+ " and "
				+ boost::lexical_cast<std::string>( prm_length_b )
				+ " so their difference is "
				+ boost::lexical_cast<std::string>( labs( length_difference ) )
				+ ")"));
		}
	}

	/// \brief Implementation function to calculate window to be added above and below the leading diagonal
	///
	/// \relates windowed_matrix
	///
	/// Where prm_length_a is less than prm_length_b, the lower part will be less than the upper part:
	///
	///     *  +  +  +  .  .
	///     +  *  +  +  +  .
	///     .  +  *  +  +  +
	///     .  .  +  *  +  +
	///
	/// Where prm_length_a is greater than prm_length_b, the lower part will be greater than the upper part:
	///
	///     *  +  .  .
	///     +  *  +  .
	///     +  +  *  +
	///     +  +  +  *
	///     .  +  +  +
	///     .  .  +  +
	inline size_size_pair get_window_upper_and_lower_part_widths(const windowed_matrix::size_type &prm_length_a,        ///< The length of the first  entry (or the number of rows    in     matrix)
	                                                             const windowed_matrix::size_type &prm_length_b,        ///< The length of the second entry (or the number of columns in the matrix)
	                                                             const windowed_matrix::size_type &prm_requested_window ///< The requested width (or equivalently, height) of the window
	                                                             ) {
		check_lengths_and_window_size_are_valid(prm_length_a, prm_length_b, prm_requested_window);

		// Calculate the length difference
		const auto diff_t_length_a   = debug_numeric_cast<windowed_matrix::difference_type>(prm_length_a);
		const auto diff_t_length_b   = debug_numeric_cast<windowed_matrix::difference_type>(prm_length_b);
		const windowed_matrix::difference_type length_difference = diff_t_length_a - diff_t_length_b;

		// Divide the window size less one (for the leading diagonal itself) between the two parts,
		// ensuring that the correct part is length_difference greater than the other
		const auto upper_width = debug_numeric_cast<size_t>( ( debug_numeric_cast<int>(prm_requested_window - 1 ) - length_difference) / 2 );
		const auto lower_width = debug_numeric_cast<size_t>( ( debug_numeric_cast<int>(prm_requested_window - 1 ) + length_difference) / 2 );

		// Depending on rounding, the extra value may have meant that the resulting window is too small, so calculate any necessary increments
		const bool   is_too_narrow   = upper_width + lower_width + 1 < prm_requested_window;
		const size_t upper_increment = is_too_narrow && (upper_width <= lower_width) ? 1 : 0;
		const size_t lower_increment = is_too_narrow && (upper_width >= lower_width) ? 1 : 0;

		// Return the upper and lower values, minus any necessary increments
		return std::make_pair(
			upper_width + upper_increment,
			lower_width + lower_increment
		);
	}

	/// \brief [Temporary] Calculate the offset that is added to the a index to calculate the a index to be used to access a given cell
	///
	/// \relates windowed_matrix
	///
	/// At present this is offset_1: the prm_index_b is in [1, prm_length_b],
	/// the maximum return value is prm_length_a and, for a square matrix, the minimum return value is 1
	///
	/// \todo Eradicate this function and move its functionality into a private method of windowed_matrix
	///       In particular, used naively, this would cause indexing trouble around the very ends of the diagonal.
	inline int get_window_indexing_offset_for_b__offset_1(const windowed_matrix::size_type &prm_length_a,         ///< The length of the first  entry (or the number of rows    in     matrix)
	                                                      const windowed_matrix::size_type &prm_length_b,         ///< The length of the second entry (or the number of columns in the matrix)
	                                                      const windowed_matrix::size_type &prm_requested_window, ///< The requested width (or equivalently, height) of the window
	                                                      const windowed_matrix::size_type &prm_index_b           ///< TODOCUMENT
	                                                      ) {
		const size_size_pair window_upper_and_lower_part_widths = get_window_upper_and_lower_part_widths(
			prm_length_a,
			prm_length_b,
			prm_requested_window
		);
		const int &leading_diagonal_index  = debug_numeric_cast<int>( prm_index_b );
		const int &window_upper_part_width = debug_numeric_cast<int>( window_upper_and_lower_part_widths.first );
		return window_upper_part_width - leading_diagonal_index;
	}

	/// \brief [Temporary] Calculate the a index to be used to access a particular cell in the matrix
	///
	/// \relates windowed_matrix
	///
	/// \todo Eradicate this function and move its functionality into a private method of windowed_matrix
	inline int get_window_matrix_a_index(const windowed_matrix::size_type &prm_length_a,         ///< The length of the first  entry (or the number of rows    in     matrix)
	                                     const windowed_matrix::size_type &prm_length_b,         ///< The length of the second entry (or the number of columns in the matrix)
	                                     const windowed_matrix::size_type &prm_requested_window, ///< The requested width (or equivalently, height) of the window
	                                     const windowed_matrix::size_type &prm_index_a,          ///< TODOCUMENT
	                                     const windowed_matrix::size_type &prm_index_b           ///< TODOCUMENT
	                                     ) {
		const int offset = get_window_indexing_offset_for_b__offset_1(
			prm_length_a,
			prm_length_b,
			prm_requested_window,
			prm_index_b + 1
		);
		const int index_a = debug_numeric_cast<int>(prm_index_a);
		return index_a + offset;
	}

	/// \brief [Temporary] offset_1 wrapper to get_window_matrix_a_index()
	///
	/// Currently just replicating existing behaviour but this involves sometimes producing answers
	/// outside the range [1, prm_length_a] (including negative numbers)
	inline int get_window_matrix_a_index__offset_1(const windowed_matrix::size_type &prm_length_a,          ///< The length of the first  entry (or the number of rows    in     matrix)
	                                               const windowed_matrix::size_type &prm_length_b,          ///< The length of the second entry (or the number of columns in the matrix)
	                                               const windowed_matrix::size_type &prm_requested_window,  ///< The requested width (or equivalently, height) of the window
	                                               const windowed_matrix::size_type &prm_index_a__offset_1, ///< TODOCUMENT
	                                               const windowed_matrix::size_type &prm_index_b__offset_1  ///< TODOCUMENT
	                                               ) {
		check_offset_1( prm_index_a__offset_1 );
		check_offset_1( prm_index_b__offset_1 );
		return get_window_matrix_a_index(
			prm_length_a,
			prm_length_b,
			prm_requested_window,
			prm_index_a__offset_1 - 1,
			prm_index_b__offset_1 - 1
		) + 1;
	}

	/// \brief Calculate the index of the start of the window in a particular column
	///
	/// At present this is offset_1: the prm_index_b  is in [1, prm_length_b],
	///                              the return value is in [1, prm_length_a]
	///
	/// The code guarantees that swapping the lengths should give the transpose so this
	/// can be used to get the index of the start of the window in a particular row.
	///
	/// \relates windowed_matrix
	inline size_t get_window_start_a_for_b__offset_1(const windowed_matrix::size_type &prm_length_a,         ///< The length of the first  entry (or the number of rows    in     matrix)
	                                                 const windowed_matrix::size_type &prm_length_b,         ///< The length of the second entry (or the number of columns in the matrix)
	                                                 const windowed_matrix::size_type &prm_requested_window, ///< The requested width (or equivalently, height) of the window
	                                                 const windowed_matrix::size_type &prm_index_b           ///< TODOCUMENT
	                                                 ) {
		if constexpr ( common::IS_IN_DEBUG_MODE ) {
			if (prm_index_b < 1 || prm_index_b > prm_length_b) {
				BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Index for B is out of range"));
			}
		}
		const size_size_pair window_upper_and_lower_part_widths = get_window_upper_and_lower_part_widths(
			prm_length_a,
			prm_length_b,
			prm_requested_window
		);
		const int    window_upper_part_width = debug_numeric_cast<int>( window_upper_and_lower_part_widths.first );
		const int    leading_diagonal_index  = debug_numeric_cast<int>( prm_index_b );
		const int    window_start_unclamped  = leading_diagonal_index - window_upper_part_width;
		const auto window_start            = debug_numeric_cast<size_t>( boost::algorithm::clamp( window_start_unclamped, 1, debug_numeric_cast<int>( prm_length_a ) ) );
	//	cerr << "prm_length_a            : " << prm_length_a            << endl;
	//	cerr << "prm_length_b            : " << prm_length_b            << endl;
	//	cerr << "prm_requested_window    : " << prm_requested_window    << endl;
	//	cerr << "prm_index_b             : " << prm_index_b             << endl;
	//	cerr << "window_upper_part_width : " << window_upper_part_width << endl;
	//	cerr << "leading_diagonal_index  : " << leading_diagonal_index  << endl;
	//	cerr << "window_start            : " << window_start            << endl;
		return window_start;
	}

	/// \brief Calculate the index of the stop of the window in a particular column
	///
	/// At present this is offset_1: the prm_index_b  is in [1, prm_length_b],
	///                              the return value is in [1, prm_length_a]
	///
	/// The code guarantees that swapping the lengths should give the transpose so this
	/// can be used to get the index of the stop of the window in a particular row.
	///
	/// Note that the stop is the index of the last cell in the window (not the STL-style
	/// one-past-the-end).
	///
	/// \relates windowed_matrix
	inline size_t get_window_stop_a_for_b__offset_1(const windowed_matrix::size_type &prm_length_a,         ///< The length of the first  entry (or the number of rows    in     matrix)
	                                                const windowed_matrix::size_type &prm_length_b,         ///< The length of the second entry (or the number of columns in the matrix)
	                                                const windowed_matrix::size_type &prm_requested_window, ///< The requested width (or equivalently, height) of the window
	                                                const windowed_matrix::size_type &prm_index_b           ///< TODOCUMENT
	                                                ) {
		if constexpr ( common::IS_IN_DEBUG_MODE ) {
			if (prm_index_b < 1 || prm_index_b > prm_length_b) {
				BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Index for B is out of range"));
			}
		}

		const size_size_pair window_upper_and_lower_part_widths = get_window_upper_and_lower_part_widths(
			prm_length_a,
			prm_length_b,
			prm_requested_window
		);
		const size_t &window_lower_part_width = window_upper_and_lower_part_widths.second;
		const size_t &leading_diagonal_index  = prm_index_b;
		const size_t  window_stop_unclamped   = leading_diagonal_index + window_lower_part_width;
		const size_t  window_stop             = boost::algorithm::clamp(window_stop_unclamped, 1_z, prm_length_a);
	//	cerr << "prm_length_a            : " << prm_length_a            << endl;
	//	cerr << "prm_length_b            : " << prm_length_b            << endl;
	//	cerr << "prm_requested_window    : " << prm_requested_window    << endl;
	//	cerr << "prm_index_b             : " << prm_index_b             << endl;
	//	cerr << "window_lower_part_width : " << window_lower_part_width << endl;
	//	cerr << "leading_diagonal_index  : " << leading_diagonal_index  << endl;
	//	cerr << "window_stop             : " << window_stop             << endl;
		return window_stop;
	}

	/// \brief Return a window size that will ensure a full matrix is used for the specified matrix dimensions
	///
	/// \relates windowed_matrix
	inline size_t get_window_width_for_full_matrix(const windowed_matrix::size_type &prm_length_a, ///< The length of the first  entry (or the number of rows    in     matrix)
	                                               const windowed_matrix::size_type &prm_length_b  ///< The length of the second entry (or the number of columns in the matrix)
	                                               ) {
		return prm_length_a + prm_length_b;
	}

} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SSAP_WINDOWED_MATRIX_HPP
