/// \file
/// \brief The rbf_model class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION_RBF_MODEL_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION_RBF_MODEL_HPP

#include <boost/filesystem/path.hpp>

#include "cath/common/type_aliases.hpp"
#include "cath/score/score_classification/value_list_scaling.hpp"

namespace cath { namespace file { class prc_scores_entry; } }
namespace cath { namespace file { class ssap_scores_entry; } }
namespace cath { namespace homcheck { class ssap_and_prc; } }

namespace rbf_model_test_suite { struct try_parse; }

namespace cath {
	namespace score {

		/// \brief A convenience type-alias for a tuple of eight doubles
		using double_octuple = std::tuple<double, double, double, double, double, double, double, double>;

		/// \brief Represent a SVM RBF model for a specific combination of SSAP and PRC scores
		///
		/// SVM: Support Vector Machine
		/// RBF: Radial Basis Function
		///
		/// This is currently only used with SVM-light but there is nothing that in principle
		/// ties it to that implementation
		///
		/// \todo Separate out scaling code into a different class
		class rbf_model final {
		private:
			/// \brief The gamma parameter of the SVM RBF model
			double gamma;

			/// \brief The b parameter of the SVM RBF model
			double b;

			/// \brief The support vectors of the SVM RBF model
			std::vector<std::pair<double, double_octuple> > support_vectors;

			/// \brief The scaling to be applied to the log10 of the PRC e-value before it's used in the SVM
			static constexpr value_list_scaling PRC_EVALUE_SCALING     { -0.00440636400685107468,  0.0160543787343339489 };
			/// \brief The scaling to be applied to the PRC reverse score before it's used in the SVM
			static constexpr value_list_scaling PRC_REVERSE_SCALING    {  0.00367107195301027847,  0.0249632892804698935 };
			/// \brief The scaling to be applied to the PRC simple score before it's used in the SVM
			static constexpr value_list_scaling PRC_SIMPLE_SCALING     {  0.00345781466113416350, -0.0055325034578146623 };
			/// \brief The scaling to be applied to the SSAP number of equivalents before it's used in the SVM
			static constexpr value_list_scaling SSAP_NUM_EQUIVS_SCALING{  0.00145985401459854005,  0.0000000000000000000 };
			/// \brief The scaling to be applied to the SSAP overlap percentage before it's used in the SVM
			static constexpr value_list_scaling SSAP_OVERLAP_PC_SCALING{  0.01041666666666666610,  0.0000000000000000000 };
			/// \brief The scaling to be applied to the SSAP RMSD before it's used in the SVM
			static constexpr value_list_scaling SSAP_RMSD_SCALING      { -0.01607975558771506870,  1.0000000000000000000 };
			/// \brief The scaling to be applied to the SSAP sequence ID before it's used in the SVM
			static constexpr value_list_scaling SSAP_SEQ_ID_PC_SCALING {  0.01041666666666666610,  0.0000000000000000000 };
			/// \brief The scaling to be applied to the SSAP score before it's used in the SVM
			static constexpr value_list_scaling SSAP_SCORE_SCALING     {  0.00892379082634302961,  0.2861859718008209490 };

		public:
			rbf_model(const double &,
			          const double &,
			          std::vector<std::pair<double, double_octuple> >);
			static double_octuple make_standard_scores(const homcheck::ssap_and_prc &);

			static double_octuple make_standard_scores(const file::prc_scores_entry &,
			                                           const file::ssap_scores_entry &);

			double get_score(const double_octuple &) const;
		};

		rbf_model parse_rbf_model(std::istream &);
		rbf_model parse_rbf_model(const boost::filesystem::path &);

		double get_score(const rbf_model &,
		                 const homcheck::ssap_and_prc &);

		double get_score(const rbf_model &,
		                 const file::prc_scores_entry &,
		                 const file::ssap_scores_entry &);

	} // namespace score
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION_RBF_MODEL_HPP

