/// \file
/// \brief The cath_assign_domains_options_block class header

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

#ifndef _CATH_TOOLS_SOURCE_CATH_ASSIGN_DOMAINS_OPTIONS_CATH_ASSIGN_DOMAINS_OPTIONS_BLOCK_H
#define _CATH_TOOLS_SOURCE_CATH_ASSIGN_DOMAINS_OPTIONS_CATH_ASSIGN_DOMAINS_OPTIONS_BLOCK_H

#include <boost/ptr_container/ptr_vector.hpp>

#include "options/options_block/options_block.h"

namespace cath { namespace opts { class pdbs_acquirer; } }

namespace cath {
	namespace opts {

		/// \brief Options block for the main cath-assign-domains options
		class cath_assign_domains_options_block final : public cath::opts::options_block {
		private:
			static const std::string PO_SVMLIGHT_RBF_FILE;
			static const std::string PO_FILELIST_FILE;
			static const std::string PO_SF_OF_DOMAIN_FILE;
			static const std::string PO_FORBIDDEN_NODES;

			static const str_vec DEFAULT_FORBIDDEN_NODES;

			/// \brief The SVM-light RBF model file
			boost::filesystem::path rbf_svm_file;

			/// \brief The file containing the list of PRC/SSAP data files
			boost::filesystem::path data_data_file;

			/// \brief The file containing superfamily of domain
			boost::filesystem::path sf_of_dom_file;

			/// \brief The list of CATH nodes forbidden for assignment
			str_vec forbidden_nodes;

			virtual std::unique_ptr<options_block> do_clone() const override final;
			virtual std::string do_get_block_name() const override final;
			virtual void do_add_visible_options_to_description(boost::program_options::options_description &) override final;
			virtual str_opt do_invalid_string(const boost::program_options::variables_map &) const override final;

		public:
			const boost::filesystem::path & get_rbf_svm_file() const;
			const boost::filesystem::path & get_data_data_file() const;
			const boost::filesystem::path & get_sf_of_dom_file() const;
			const str_vec & get_forbidden_nodes() const;
		};

	} // namespace opts
} // namespace cath

#endif
