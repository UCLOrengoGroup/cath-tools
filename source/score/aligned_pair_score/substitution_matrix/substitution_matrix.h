/// \file
/// \brief The substitution_matrix class header

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

#ifndef SUBSTITUTION_MATRIX_H_INCLUDED
#define SUBSTITUTION_MATRIX_H_INCLUDED

//#include <boost/serialization/nvp.hpp>
//#include <boost/serialization/vector.hpp>

#include "score/score_type_aliases.h"
#include "structure/structure_type_aliases.h"

namespace cath { namespace score { class substitution_matrix; } }
namespace cath { namespace score { bool operator<(const substitution_matrix &, const substitution_matrix &); } }

namespace cath {
	namespace score {

		/// \brief Represent a substitution matrix of scores for a substitution of one amino acid for another
		///
		/// Invariants:
		///  * amino_acids are sorted
		///  * scores' indices correspond to those of amino_acids (for both dimensions)
		///  * scores is symmetric
		class substitution_matrix final {
		private:
			friend bool operator<(const substitution_matrix &,
			                      const substitution_matrix &);

//			friend class boost::serialization::access;
//
//			template<class archive> void serialize(archive &ar,
//			                                       const size_t /*version*/
//			                                       ) {
//				ar & BOOST_SERIALIZATION_NVP( amino_acids );
//				ar & BOOST_SERIALIZATION_NVP( scores );
//				ar & BOOST_SERIALIZATION_NVP( score_for_one_unknown_aa );
//				ar & BOOST_SERIALIZATION_NVP( score_for_two_unknown_aas );
//				ar & BOOST_SERIALIZATION_NVP( name );
//			}

			/// \brief The sorted list of amino acids for which the substitution scores are provided
			amino_acid_vec amino_acids;

			/// \brief The all-against-all scores, with both dimensions sorted in the same order as the sorted amino_acids list
			diff_vec_vec scores;

			/// \brief The score that's returned if one of the query amino acids is unrecognised
			ptrdiff_t score_for_one_unknown_aa;

			/// \brief The score that's returned if both query amino acids are unrecognised
			ptrdiff_t score_for_two_unknown_aas;

			/// \brief A name for the substitution matrix
			std::string name;

			const amino_acid_vec & get_amino_acids() const;
			const diff_vec_vec & get_scores() const;
			const ptrdiff_t & get_score_for_one_unknown_aa() const;
			const ptrdiff_t & get_score_for_two_unknown_aas() const;

			static size_vec order_permutation(const amino_acid_vec &,
			                                  const amino_acid_vec &);

			static diff_vec_vec reorder_scores(const amino_acid_vec &,
			                                   const amino_acid_vec &,
			                                   const diff_vec_vec &);

			static bool has_lower_highest_score(const diff_vec &,
			                                    const diff_vec &);
			static ptrdiff_t highest_score_of_scores(const diff_vec &);

			void check_is_symmetric() const;

		public:
			substitution_matrix(const amino_acid_vec &,
			                    const diff_vec_vec &,
			                    const ptrdiff_t &,
			                    const ptrdiff_t &,
			                    const std::string &);

			const std::string & get_name() const;

			ptrdiff_t get_highest_score() const;

			ptrdiff_t get_score(const amino_acid &,
			                    const amino_acid &) const;
		};

		bool operator<(const substitution_matrix &,
		               const substitution_matrix &);

		amino_acid_vec get_standard_substitution_amino_acids();

		substitution_matrix_vec get_all_substitution_matrices();
	}
}

#endif
