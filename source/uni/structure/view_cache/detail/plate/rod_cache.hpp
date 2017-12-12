/// \file
/// \brief The rod_cache class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_STRUCTURE_VIEW_CACHE_DETAIL_PLATE_ROD_CACHE_H
#define _CATH_TOOLS_SOURCE_UNI_STRUCTURE_VIEW_CACHE_DETAIL_PLATE_ROD_CACHE_H

#include "common/exception/invalid_argument_exception.hpp"
#include "common/type_aliases.hpp"
#include "structure/view_cache/index/view_cache_index_entry.hpp"

namespace cath {
	namespace index {
		namespace detail {

			/// \brief TODOCUMENT
			///
			/// \todo Eliminate redundancy between this and view_cache
			///
			/// \todo If profiles reveal that this is critical in performance:
			///        * ensure that access the entries is inlined
			///        * consider trying giving access to a particular rod (ie sub-vector)
			///          so that the using code can directly index into that vector
			///          And then profile because it might well all have been optimised away anyway
			class rod_cache final {
			public:
				const view_cache_index_entry & get_entry(const size_t &,
				                                         const size_t &);
			};

			namespace detail {

				/// \brief Check that an element is in-range in a conversion between indices and rod/notch
				///
				/// This is used conversions in both directions
				inline void check_indices(const size_t &arg_index_a, ///< The index to check in the first dimension
				                          const size_t &arg_index_b, ///< The index to check in the second dimension
				                          const size_t &arg_size_a,  ///< The size of the matrix in the first dimension
				                          const size_t &arg_size_b   ///< The size of the matrix in the second dimension
				                          ) {
					// Check that both indices are in range
					if ( arg_index_a >= arg_size_a) {
						BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Unable to convert between rod/notch and indices because the input specify a location that is out of range with respect to size_a"));
					}
					if ( arg_index_b >= arg_size_b) {
						BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Unable to convert between rod/notch and indices because the input specify a location that is out of range with respect to size_a"));
					}
				}

				/// \brief Perform the rod/notch -> index_a conversion without any out-of-range checks
				inline size_t get_index_a_of_rod_and_notch_unchecked(const size_t &arg_rod_index,    ///< The rod of the element
				                                                     const size_t &arg_notch_index,  ///< The notch of the element
				                                                     const size_t &arg_size_a        ///< The size of the matrix in the first dimension
				                                                     ) {
					return std::max( arg_size_a - 1, arg_rod_index ) + arg_notch_index - arg_rod_index;
				}

				/// \brief Perform the rod/notch -> index_b conversion without any out-of-range checks
				inline size_t get_index_b_of_rod_and_notch_unchecked(const size_t &arg_rod_index,    ///< The rod of the element
				                                                     const size_t &arg_notch_index,  ///< The notch of the element
				                                                     const size_t &arg_size_a        ///< The size of the matrix in the first dimension
				                                                     ) {
					return std::max( arg_size_a - 1, arg_rod_index ) + arg_notch_index - ( arg_size_a - 1 );
				}

			} // namespace detail

			/// \brief Get the rod of the entry in the matrix of specified dimensions with the specified indices
			///
			/// The rod is the number of the upper-left to lower-right diagonal, numbered starting from 0 for the bottom left element.
			/// The notch is the number along that rod, starting from 0 at the top left.
			///
			/// Examples rod/cache numbers for a 3x3 matrix
			///
			///     (2,0)  (3,0)  (4,0)
			///     (1,0)  (2,1)  (3,1)
			///     (0,0)  (1,1)  (2,2)
			///
			/// \pre If arg_check is true then the values must specify a valid position in the matrix,
			///      else an invalid_argument_exception will be thrown
			inline size_t get_rod_of_indices(const size_t &arg_index_a,      ///< The index in the first dimension
			                                 const size_t &arg_index_b,      ///< The index in the second dimension
			                                 const size_t &arg_size_a,       ///< The size of the matrix in the first dimension
			                                 const size_t &arg_size_b,       ///< The size of the matrix in the second dimension
#ifndef NDEBUG
			                                 const bool   &arg_check = true  ///< Whether to check for out-of-range errors (default: true if NDEBUG not defined)
#else
			                                 const bool   &arg_check = false ///< Whether to check for out-of-range errors (default: false if NDEBUG defined)
#endif
			                                 ) {
				if ( arg_check ) {
					detail::check_indices( arg_index_a, arg_index_b, arg_size_a, arg_size_b );
				}
				return ( arg_size_a - 1 - arg_index_a ) + arg_index_b;
			}

			/// \brief Get the notch of the entry in the matrix of specified dimensions with the specified indices
			///
			/// The rod is the number of the upper-left to lower-right diagonal, numbered starting from 0 for the bottom left element.
			/// The notch is the number along that rod, starting from 0 at the top left.
			///
			/// Examples rod/cache numbers for a 3x3 matrix
			///
			///     (2,0)  (3,0)  (4,0)
			///     (1,0)  (2,1)  (3,1)
			///     (0,0)  (1,1)  (2,2)
			///
			/// \pre If arg_check is true then the values must specify a valid position in the matrix,
			///      else an invalid_argument_exception will be thrown
			inline size_t get_notch_of_indices(const size_t &arg_index_a,      ///< The index in the first dimension
			                                   const size_t &arg_index_b,      ///< The index in the second dimension
			                                   const size_t &arg_size_a,       ///< The size of the matrix in the first dimension
			                                   const size_t &arg_size_b,       ///< The size of the matrix in the second dimension
#ifndef NDEBUG
			                                   const bool   &arg_check = true  ///< Whether to check for out-of-range errors (default: true if NDEBUG not defined)
#else
			                                   const bool   &arg_check = false ///< Whether to check for out-of-range errors (default: false if NDEBUG defined)
#endif
			                                   ) {
				if ( arg_check ) {
					detail::check_indices( arg_index_a, arg_index_b, arg_size_a, arg_size_b );
				}
				return std::min( arg_index_a, arg_index_b );
			}

			/// \brief Get the rod and notch of the entry in the matrix of specified dimensions with the specified indices
			///
			/// The rod is the number of the upper-left to lower-right diagonal, numbered starting from 0 for the bottom left element.
			/// The notch is the number along that rod, starting from 0 at the top left.
			///
			/// Examples rod/cache numbers for a 3x3 matrix
			///
			///     (2,0)  (3,0)  (4,0)
			///     (1,0)  (2,1)  (3,1)
			///     (0,0)  (1,1)  (2,2)
			///
			/// \pre If arg_check is true then the values must specify a valid position in the matrix,
			///      else an invalid_argument_exception will be thrown
			inline size_size_pair get_rod_and_notch_of_indices(const size_t &arg_index_a,      ///< The index in the first dimension
			                                                   const size_t &arg_index_b,      ///< The index in the second dimension
			                                                   const size_t &arg_size_a,       ///< The size of the matrix in the first dimension
			                                                   const size_t &arg_size_b,       ///< The size of the matrix in the second dimension
#ifndef NDEBUG
			                                                   const bool   &arg_check = true  ///< Whether to check for out-of-range errors (default: true if NDEBUG not defined)
#else
			                                                   const bool   &arg_check = false ///< Whether to check for out-of-range errors (default: false if NDEBUG defined)
#endif
			                                                   ) {
				if ( arg_check ) {
					detail::check_indices( arg_index_a, arg_index_b, arg_size_a, arg_size_b );
				}
				return std::make_pair(
					get_rod_of_indices  ( arg_index_a, arg_index_b, arg_size_a, arg_size_b, false ),
					get_notch_of_indices( arg_index_a, arg_index_b, arg_size_a, arg_size_b, false )
				);
			}

			/// \brief Get the index_a of the entry in the matrix of specified dimensions with the specified rod and notch
			///
			/// The rod is the number of the upper-left to lower-right diagonal, numbered starting from 0 for the bottom left element.
			/// The notch is the number along that rod, starting from 0 at the top left.
			///
			/// Examples rod/cache numbers for a 3x3 matrix
			///
			///     (2,0)  (3,0)  (4,0)
			///     (1,0)  (2,1)  (3,1)
			///     (0,0)  (1,1)  (2,2)
			///
			/// \pre If arg_check is true then the values must specify a valid position in the matrix,
			///      else an invalid_argument_exception will be thrown
			inline size_t get_index_a_of_rod_and_notch(const size_t &arg_rod_index,    ///< The rod of the element
			                                           const size_t &arg_notch_index,  ///< The notch of the element
			                                           const size_t &arg_size_a,       ///< The size of the matrix in the first dimension
			                                           const size_t &arg_size_b,       ///< The size of the matrix in the second dimension
#ifndef NDEBUG
			                                           const bool   &arg_check = true  ///< Whether to check for out-of-range errors (default: true if NDEBUG not defined)
#else
			                                           const bool   &arg_check = false ///< Whether to check for out-of-range errors (default: false if NDEBUG defined)
#endif
			                                           ) {
				const size_t index_a = detail::get_index_a_of_rod_and_notch_unchecked( arg_rod_index, arg_notch_index, arg_size_a );
				if ( arg_check ) {
					const size_t index_b = detail::get_index_b_of_rod_and_notch_unchecked( arg_rod_index, arg_notch_index, arg_size_a );
					detail::check_indices( index_a, index_b, arg_size_a, arg_size_b );
				}
				return index_a;
			}

			/// \brief Get the index_b of the entry in the matrix of specified dimensions with the specified rod and notch
			///
			/// The rod is the number of the upper-left to lower-right diagonal, numbered starting from 0 for the bottom left element.
			/// The notch is the number along that rod, starting from 0 at the top left.
			///
			/// Examples rod/cache numbers for a 3x3 matrix
			///
			///     (2,0)  (3,0)  (4,0)
			///     (1,0)  (2,1)  (3,1)
			///     (0,0)  (1,1)  (2,2)
			///
			/// \pre If arg_check is true then the values must specify a valid position in the matrix,
			///      else an invalid_argument_exception will be thrown
			inline size_t get_index_b_of_rod_and_notch(const size_t &arg_rod_index,    ///< The rod of the element
			                                           const size_t &arg_notch_index,  ///< The notch of the element
			                                           const size_t &arg_size_a,       ///< The size of the matrix in the first dimension
			                                           const size_t &arg_size_b,       ///< The size of the matrix in the second dimension
#ifndef NDEBUG
			                                           const bool   &arg_check = true  ///< Whether to check for out-of-range errors (default: true if NDEBUG not defined)
#else
			                                           const bool   &arg_check = false ///< Whether to check for out-of-range errors (default: false if NDEBUG defined)
#endif
			                                           ) {
				const size_t index_b = detail::get_index_b_of_rod_and_notch_unchecked( arg_rod_index, arg_notch_index, arg_size_a );
				if ( arg_check ) {
					const size_t index_a = detail::get_index_a_of_rod_and_notch_unchecked( arg_rod_index, arg_notch_index, arg_size_a );
					detail::check_indices( index_a, index_b, arg_size_a, arg_size_b );
				}
				return index_b;
			}

			/// \brief Get the indices of the entry in the matrix of specified dimensions with the specified rod and notch
			///
			/// The rod is the number of the upper-left to lower-right diagonal, numbered starting from 0 for the bottom left element.
			/// The notch is the number along that rod, starting from 0 at the top left.
			///
			/// Examples rod/cache numbers for a 3x3 matrix
			///
			///     (2,0)  (3,0)  (4,0)
			///     (1,0)  (2,1)  (3,1)
			///     (0,0)  (1,1)  (2,2)
			///
			/// \pre If arg_check is true then the values must specify a valid position in the matrix,
			///      else an invalid_argument_exception will be thrown
			inline size_size_pair get_indices_of_rod_and_notch(const size_t &arg_rod_index,    ///< The rod of the element
			                                                   const size_t &arg_notch_index,  ///< The notch of the element
			                                                   const size_t &arg_size_a,       ///< The size of the matrix in the first dimension
			                                                   const size_t &arg_size_b,       ///< The size of the matrix in the second dimension
#ifndef NDEBUG
			                                                   const bool   &arg_check = true  ///< Whether to check for out-of-range errors (default: true if NDEBUG not defined)
#else
			                                                   const bool   &arg_check = false ///< Whether to check for out-of-range errors (default: false if NDEBUG defined)
#endif
			                                                   ) {
				const size_t index_a = detail::get_index_a_of_rod_and_notch_unchecked( arg_rod_index, arg_notch_index, arg_size_a );
				const size_t index_b = detail::get_index_b_of_rod_and_notch_unchecked( arg_rod_index, arg_notch_index, arg_size_a );
				if ( arg_check ) {
					detail::check_indices( index_a, index_b, arg_size_a, arg_size_b );
				}
				return std::make_pair( index_a, index_b );
			}

		} // namespace detail
	} // namespace index
} // namespace cath

#endif
