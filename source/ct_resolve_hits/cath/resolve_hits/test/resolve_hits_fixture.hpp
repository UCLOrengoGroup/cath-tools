/// \file
/// \brief The resolve_hits_fixture header

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

#ifndef _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_TEST_RESOLVE_HITS_FIXTURE_HPP
#define _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_TEST_RESOLVE_HITS_FIXTURE_HPP

#include <filesystem>
#include <regex>
#include <string>
#include <string_view>

namespace cath {
	namespace rslv {

		/// \brief Store data that can be reused for multiple resolve_hits tests
		class resolve_hits_fixture {
		  protected:
			static std::regex              output_version_regex();
			static std::string             output_version_blankstr();
			static std::string             blank_vrsn( const std::string & );
			static std::string             blank_vrsn( const std::ostringstream & );
			static ::std::filesystem::path blank_vrsn( const ::std::filesystem::path & );

			/// \brief Example resolve_hits input
			static constexpr ::std::string_view EXAMPLE_INPUT_RAW =
			  R"(qyikaz cath|current|2i24N00/2-114-i5_1,6.7e-18 2744.56644492722 8-108
qyikaz cath|current|3k3qA00/4-140-i5_1,1.5e-15 2507.21478283538 8-108
qyikaz cath|current|1mfaH01/251-347-i5_1,1.6e-07 1494.83332154362 9-57,135-174
qyikaz cath|current|1yjdC00/1-118-i5_1,1.5e-09 1581.20833423932 9-92
qyikaz cath|current|1hxmA02/122-206-i5_3,9.2e-06 1037.49863991316 839-907
qyikaz cath|current|3gbnL02/124-195-i5_3,1.1e-12 1559.06111935377 839-909
qyikaz cath|current|1hzhH01/24-128-i5_3,8e-05 1212.33426111869 857-910,968-999
qyikaz cath|current|1gxeA00/1-130-i5_4,2.7e-18 3418.51089324429 930-1053
qyikaz cath|current|1mkfA01/12-210-i5_4,3.5e-15 2470.04912752062 953-1053
qyikaz cath|current|1bwmA02/202-392-i5_4,2.2e-10 1651.23649481093 954-1037
qyikaz cath|current|3tv3L01/3-107-i5_7,1.1e-09 1744.19187296544 1424-1515
qyikaz cath|current|4hbcL01/1-110-i5_7,3.3e-09 1700.29671753123 1424-1515
qyikaz cath|current|2y23A02/1245-1358-i5_7,1.4e-13 2239.67945250353 1424-1521
qyikaz cath|current|1mfaH01/251-347-i5_6,0.63 877.256712746992 1425-1472,1601-1638
qyikaz cath|current|2vscA01/2-104-i5_7,8.1e-06 1162.04665354634 1425-1501
qyikaz cath|current|2p45B00/1-121-i5_7,3.8e-07 1346.45774507742 1425-1506
qyikaz cath|current|2qhlD00/3-111-i5_7,6.8e-07 1325.73426915809 1425-1506
qyikaz cath|current|1hzhH01/24-128-i5_8,4.1 572.620184740096 1634-1653,1700-1740
qyikaz cath|current|1bwmA02/202-392-i5_9,1.1e-11 1697.64719250218 1694-1774
qyikaz cath|current|3uzqA02/122-235-i5_9,2.8e-11 1849.7557771792 1697-1786
qyikaz cath|current|2eo1A00/1-102-i5_9,4e-14 2292.99812084986 1697-1794
qyikaz cath|current|1neuA00/1-119-i5_9,2.5e-09 1655.58333922819 1710-1798
qyikaz cath|current|1xedC00/3-112-i5_9,3e-10 1737.53620832995 1710-1798
qyikaz cath|current|1fo0B00/1-116-i5_8,0.0025 441.072099696479 1711-1745
qyikaz cath|current|1mfaH01/251-347-i5_7,2.4e-07 1529.02056576253 1711-1757,1967-2011
qyikaz cath|current|2q8bH01/1-112-i5_9,6.9e-07 1212.08631819471 1711-1785
qyikaz cath|current|2qhlD00/3-111-i5_9,9.9e-06 1140.33172521059 1711-1786
qyikaz cath|current|2p45B00/1-121-i5_9,5.2e-10 1523.43573585285 1711-1789
qyikaz cath|current|1epfA01/1-96-i5_11,0.00046 1080.31661563379 2037-2117
qyikaz cath|current|1nctA00/-6-91-i5_11,4e-07 1492.21254078916 2037-2127
iexvva cath|current|3k3qA00/4-140-i5_1,1.5e-15 2507.21478283538 8-108
iexvva cath|current|1mfaH01/251-347-i5_1,1.5e-07 1497.32787794404 9-57,135-174
iexvva cath|current|1yjdC00/1-118-i5_1,1.4e-09 1583.72524500303 9-92
iexvva cath|current|3q5yA01/3-117-i5_1,1.3e-11 1754.42875840623 9-92
iexvva cath|current|1gxeA00/1-130-i5_3,2.2e-26 3957.99108242873 852-962
iexvva cath|current|3iagC02/201-358-i5_2,1.7 498.247105009708 862-912
iexvva cath|current|2wqrA02/329-435-i5_3,8.8e-18 2597.32966347358 865-960
iexvva cath|current|1bwmA02/202-392-i5_3,8.6e-15 1949.30562544927 867-947
iexvva cath|current|1xiwB00/1-74-i5_3,3.9e-18 2357.16844379572 883-968
iexvva cath|current|1hxmA02/122-206-i5_3,9.1e-06 1037.82614392984 885-953
iexvva cath|current|3gbnL02/124-195-i5_3,1.1e-12 1559.06111935377 885-955
iexvva cath|current|1hzhH01/24-128-i5_3,7.9e-05 1212.80407014902 903-956,1014-1045
iexvva cath|current|2xotA02/278-363-i5_8,6.4e-08 1375.50560208129 1630-1709
iexvva cath|current|1xiwB00/1-74-i5_8,2.3e-07 1397.61486177452 1630-1713
iexvva cath|current|1u58A02/143-242-i5_8,0.032 919.588001734407 1631-1710
iexvva cath|current|4hbqA02/91-210-i5_8,2.3e+02 534.679051478769 1634-1703
iexvva cath|current|2pttB00/2-109-i5_8,1.1e-05 1091.97833398345 1635-1707
iexvva cath|current|4fa8A01/20-123-i5_8,9.7 441.648185020954 1657-1705
iexvva cath|current|1hzhH01/24-128-i5_8,4 573.274340528994 1680-1699,1746-1786
iexvva cath|current|1neuA00/1-119-i5_9,2.4e-09 1657.16119948767 1756-1844
iexvva cath|current|1xedC00/3-112-i5_9,2.9e-10 1738.84657818699 1756-1844
iexvva cath|current|1fo0B00/1-116-i5_8,0.0024 441.692606540094 1757-1791
iexvva cath|current|1mfaH01/251-347-i5_7,5.2e-08 1503.70770910377 1757-1803,2207-2246
iexvva cath|current|2q8bH01/1-112-i5_9,6.8e-07 1212.56183154703 1757-1831
iexvva cath|current|2qhlD00/3-111-i5_9,9.8e-06 1140.66681824737 1757-1832
iexvva cath|current|2p45B00/1-121-i5_9,5.1e-10 1524.10195608826 1757-1835
iexvva cath|current|1epfA01/1-96-i5_11,0.00045 1081.0897863842 2083-2163
iexvva cath|current|1nctA00/-6-91-i5_11,4e-07 1492.21254078916 2083-2173
iexvva cath|current|2kkqA00/1-116-i5_11,4.3e-05 1307.35437054226 2083-2173
iexvva cath|current|1p53A03/368-450-i5_11,0.00016 1103.67040138753 2084-2124,2216-2254
)";

			/// \brief Example resolve_hits output
			static constexpr ::std::string_view EXAMPLE_OUTPUT =
			  R"(# Generated by cath-resolve-hits vX.X.X-X-XXXXXXXX, one of the cath-tools (https://github.com/UCLOrengoGroup/cath-tools)
#FIELDS query-id match-id score boundaries resolved
iexvva cath|current|3k3qA00/4-140-i5_1,1.5e-15 2507.21 8-108 8-108
iexvva cath|current|1gxeA00/1-130-i5_3,2.2e-26 3957.99 852-962 852-962
iexvva cath|current|1xiwB00/1-74-i5_8,2.3e-07 1397.61 1630-1713 1630-1713
iexvva cath|current|1xedC00/3-112-i5_9,2.9e-10 1738.85 1756-1844 1756-1844
iexvva cath|current|1nctA00/-6-91-i5_11,4e-07 1492.21 2083-2173 2083-2173
qyikaz cath|current|2i24N00/2-114-i5_1,6.7e-18 2744.57 8-108 8-108
qyikaz cath|current|3gbnL02/124-195-i5_3,1.1e-12 1559.06 839-909 839-909
qyikaz cath|current|1gxeA00/1-130-i5_4,2.7e-18 3418.51 930-1053 930-1053
qyikaz cath|current|2y23A02/1245-1358-i5_7,1.4e-13 2239.68 1424-1521 1424-1521
qyikaz cath|current|2eo1A00/1-102-i5_9,4e-14 2293 1697-1794 1697-1794
qyikaz cath|current|1nctA00/-6-91-i5_11,4e-07 1492.21 2037-2127 2037-2127
)";
		};

	} // namespace rslv
} // namespace cath
#endif // _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_TEST_RESOLVE_HITS_FIXTURE_HPP
