/// \file
/// \brief The eigen class header

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

#ifndef _CATH_TOOLS_SOURCE_EIGEN_EIGEN_H
#define _CATH_TOOLS_SOURCE_EIGEN_EIGEN_H

#include "exception/runtime_error_exception.hpp"

#include <iostream> // ***** TEMPORARY *****

#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>

namespace cath {
	namespace chop {

		/// \brief TODOCUMENT
		class gsl_matrix_wrp final {
		private:
			/// \brief TODOCUMENT
			gsl_matrix * ptr;

		public:
			/// \brief TODOCUMENT
			gsl_matrix_wrp() = delete;

			/// \brief TODOCUMENT
			gsl_matrix_wrp(const size_t &arg_n1, ///< TODOCUMENT
			               const size_t &arg_n2  ///< TODOCUMENT
			               ) : ptr{ gsl_matrix_alloc( arg_n1, arg_n2 ) } {
				if ( ptr == nullptr ) {
					BOOST_THROW_EXCEPTION(common::runtime_error_exception("Was unable to allocate gsl_matrix"));
				}
			}

			/// \brief TODOCUMENT
			~gsl_matrix_wrp() noexcept {
				try {
					gsl_matrix_free( ptr );
				}
				catch (...) {
				}
			}

			/// \brief TODOCUMENT
			const gsl_matrix & get_ref() const {
				return *ptr;
			}

			/// \brief TODOCUMENT
			gsl_matrix & get_ref() {
				return *ptr;
			}

			/// \brief TODOCUMENT
			const gsl_matrix * get_ptr() const {
				return ptr;
			}

			/// \brief TODOCUMENT
			gsl_matrix * get_ptr() {
				return ptr;
			}
		};

		/// \brief TODOCUMENT
		class gsl_vector_wrp final {
		private:
			/// \brief TODOCUMENT
			gsl_vector * ptr;

		public:
			/// \brief TODOCUMENT
			gsl_vector_wrp() = delete;

			/// \brief TODOCUMENT
			gsl_vector_wrp(const size_t &arg_n ///< TODOCUMENT
			               ) : ptr{ gsl_vector_alloc( arg_n ) } {
				if ( ptr == nullptr ) {
					BOOST_THROW_EXCEPTION(common::runtime_error_exception("Was unable to allocate gsl_vector"));
				}
			}

			/// \brief TODOCUMENT
			~gsl_vector_wrp() noexcept {
				try {
					gsl_vector_free( ptr );
				}
				catch (...) {
				}
			}

			/// \brief TODOCUMENT
			const gsl_vector & get_ref() const {
				return *ptr;
			}

			/// \brief TODOCUMENT
			gsl_vector & get_ref() {
				return *ptr;
			}

			/// \brief TODOCUMENT
			const gsl_vector * get_ptr() const {
				return ptr;
			}

			/// \brief TODOCUMENT
			gsl_vector * get_ptr() {
				return ptr;
			}
		};

			/// \brief TODOCUMENT
		void pretty_print(const gsl_matrix_wrp &arg_matrix ///< TODOCUMENT
		                  ) {
			// Get the dimension of the matrix.
			size_t rows = arg_matrix.get_ref().size1;
			size_t cols = arg_matrix.get_ref().size2;
			// Now print out the data in a square format.
			for (size_t i = 0; i < rows; ++i) {
				for (size_t j = 0; j < cols; ++j) {
					std::cout << gsl_matrix_get( arg_matrix.get_ptr(), i, j ) << " ";
				}
				std::cout << "\n";
			}
		}


		int my_function() {
			double a[] = { 1.0, 1.0, 0.0,
			               1.0, 2.0, 1.0,
			               0.0, 1.0, 4.0 };
			gsl_matrix_view A = gsl_matrix_view_array( a, 3, 3 );

			gsl_matrix_wrp  V    { 3, 3 };
			gsl_vector_wrp  S    { 3    };
			 gsl_vector_wrp work { 3    };

			// From the gsl documentation: The gsl_linalg_SV_decomp function
			// factorizes the M-by-N matrix A into the singular value
			// decomposition A = U S V^T for M >= N. On output the matrix A is
			// replaced by U. The diagonal elements of the singular value matrix
			// S are stored in the vector S. The singular values are non-negative
			// and form a non-increasing sequence from S_1 to S_N. The matrix V
			// contains the elements of V in untransposed form. To form the
			// product U S V^T it is necessary to take the transpose of V. A
			// workspace of length N is required in work.
			gsl_linalg_SV_decomp(
				&A.matrix,
				V.get_ptr(),
				S.get_ptr(),
				work.get_ptr()
			);

			//  gsl_matrix_fprintf (stdout, &A.matrix, "%g"); cout<<"\n";
			// pretty_print( &A.matrix );
			std::cout<<"\n";
			//  gsl_matrix_fprintf (stdout, V, "%g"); cout<<"\n";
			pretty_print( V );
			std::cout<<"\n";
			gsl_vector_fprintf( stdout, S.get_ptr(), "%g" );

			double data[] = { 1.0  , 1/2.0, 1/3.0, 1/4.0,
			                  1/2.0, 1/3.0, 1/4.0, 1/5.0,
			                  1/3.0, 1/4.0, 1/5.0, 1/6.0,
			                  1/4.0, 1/5.0, 1/6.0, 1/7.0 };

			gsl_matrix_view m  = gsl_matrix_view_array (data, 4, 4);

			gsl_vector *eval = gsl_vector_alloc (4);
			gsl_matrix *evec = gsl_matrix_alloc (4, 4);

			gsl_eigen_symmv_workspace * w =  gsl_eigen_symmv_alloc( 4 );

			gsl_eigen_symmv (&m.matrix, eval, evec, w);

			gsl_eigen_symmv_free (w);

			gsl_eigen_symmv_sort (eval, evec,  GSL_EIGEN_SORT_ABS_ASC);

			for (size_t i = 0; i < 4; ++i) {
				double          eval_i  = gsl_vector_get   ( eval, i );
				gsl_vector_view evec_i  = gsl_matrix_column( evec, i );

				printf( "eigenvalue = %g\n", eval_i );
				printf( "eigenvector = \n" );
				gsl_vector_fprintf( stdout,  &evec_i.vector, "%g" );
			}

			gsl_vector_free( eval );
			gsl_matrix_free( evec );

			return 0;
		}

	} // namespace chop
} // namespace cath

#endif
