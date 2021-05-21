/// \file
/// \brief The resolve_hits_fixture class definitions

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

#include "resolve_hits_fixture.hpp"

#include <filesystem>
#include <sstream>

#include "cath/common/regex/regex_replace_file.hpp"

using namespace ::cath::common;
using namespace ::cath::rslv;

using ::std::filesystem::path;
using ::std::ostringstream;
using ::std::regex;
using ::std::string;

/// \brief A regular expression for identifying the version region of the standard CRH output
///
/// To use, see the overloads of blank_vrsn()
regex resolve_hits_fixture::output_version_regex() {
	return regex{ R"(by cath-resolve-hits( v\d+\.\d+\.\d+\-\d+\-g?\w{8})?,)" };
}

/// \brief A string for use in blanking out the version in the standard CRH output
///
/// To use, see the overloads of blank_vrsn()
string resolve_hits_fixture::output_version_blankstr() {
	return "by cath-resolve-hits vX.X.X-X-XXXXXXXX,";
}

/// \brief Blank out the version in the standard CRH output in the specified string
///
/// \returns The string in which the version has been blanked out
string resolve_hits_fixture::blank_vrsn(const string &prm_string ///< The string in which the CRH version should be blanked out
                                        ) {
	return regex_replace( prm_string, output_version_regex(), output_version_blankstr() );
}

/// \brief Blank out the version in the standard CRH output in the specified ostringstream
///
/// \returns The ostringstream in which the version has been blanked out
string resolve_hits_fixture::blank_vrsn(const ostringstream &prm_stream ///< The ostringstream in which the CRH version should be blanked out
                                        ) {
	return blank_vrsn( prm_stream.str() );
}

/// \brief Blank out the version in the standard CRH output in the specified file
///
/// \returns The file in which the version has been blanked out
path resolve_hits_fixture::blank_vrsn(const path &prm_file ///< The file in which the CRH version should be blanked out
                                      ) {
	regex_replace_file( prm_file, output_version_regex(), output_version_blankstr() );
	return prm_file;
}
