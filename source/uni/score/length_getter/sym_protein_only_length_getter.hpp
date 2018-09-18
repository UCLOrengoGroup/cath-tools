/// \file
/// \brief The sym_protein_only_length_getter class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_SCORE_LENGTH_GETTER_SYM_PROTEIN_ONLY_LENGTH_GETTER_HPP
#define _CATH_TOOLS_SOURCE_UNI_SCORE_LENGTH_GETTER_SYM_PROTEIN_ONLY_LENGTH_GETTER_HPP

#include <boost/ptr_container/ptr_vector.hpp>

#include "score/length_getter/protein_only_length_getter.hpp"

namespace cath {
	namespace score {


		/// \brief TODOCUMENT
		class sym_protein_only_length_getter : public protein_only_length_getter {
		private:
			/// \brief TODOCUMENT
			virtual std::unique_ptr<sym_protein_only_length_getter> do_sym_protein_only_clone() const = 0;

			std::unique_ptr<protein_only_length_getter> do_protein_only_clone() const final;

		public:
			std::unique_ptr<sym_protein_only_length_getter> sym_protein_only_clone() const;
		};

		boost::ptr_vector<sym_protein_only_length_getter> get_all_sym_protein_only_length_getters();

		/// \brief Function to make sym_protein_only_length_getter meet the protein_only_length_getter concept (used in ptr_container)
        ///
        /// NOTE: Don't call this yourself. Call the object's clone() method instead because that returns a
        ///       smart pointer rather than the raw pointer this has to return to meet the Clonable concept.
        ///
        /// This gets the smart pointer from the clone() method and then calls release on it.
        ///
        /// \returns A raw pointer to a new copy of the sym_protein_only_length_getter argument, with the same dynamic type.
        ///          The caller is responsible for deleting this new object.
        inline sym_protein_only_length_getter * new_clone(const sym_protein_only_length_getter &prm_sym_protein_only_length_getter ///< The sym_protein_only_length_getter to clone
                                                          ) {
                return prm_sym_protein_only_length_getter.sym_protein_only_clone().release();
        }
	} // namespace score
} // namespace cath

#endif
