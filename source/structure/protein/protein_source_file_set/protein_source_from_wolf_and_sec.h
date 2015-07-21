/// \file
/// \brief The protein_source_from_wolf_and_sec class header

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

#ifndef PROTEIN_SOURCE_FROM_WOLF_AND_SEC_H_INCLUDED
#define PROTEIN_SOURCE_FROM_WOLF_AND_SEC_H_INCLUDED

#include "structure/protein/protein_source_file_set/protein_source_file_set.h"

namespace cath {

	/// \brief Concrete protein_source_file_set for reading each protein from a wolf file and a sec file
	class protein_source_from_wolf_and_sec final : public protein_source_file_set {
	private:
		virtual std::unique_ptr<protein_source_file_set> do_clone() const override final;

		virtual opts::data_file_vec do_get_file_set() const override final;

		virtual protein_file_combn do_get_protein_file_combn() const override final;

		virtual protein do_read_files(const opts::data_file_path_map &,
		                              const std::string &,
		                              std::ostream &) const override final;

	public:
		virtual ~protein_source_from_wolf_and_sec() noexcept = default;
	};

}

#endif
