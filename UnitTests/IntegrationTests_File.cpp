#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include "catch.hpp"

namespace novac
{
    static std::string GetTestDataDirectory()
    {
#ifdef _MSC_VER
        return std::string("../TestData/");
#else
        return std::string("TestData/");
#endif // _MSC_VER 
    }

    static std::string GetHighResolutionSO2CrossSectionFile()
    {
        return GetTestDataDirectory() + std::string("SO2_Bogumil(2003)_293K_239-395nm.txt");
    }

    static std::string GetQDoasConvolvedSO2CrossSectionFile()
    {
        return GetTestDataDirectory() + std::string("SO2_QDOAS.xs");
    }


    TEST_CASE("ReadCrossSectionFile with simple two column cross section file", "[ReadCrossSectionFile][IntegrationTest]")
    {
        CCrossSectionData result;
        ReadCrossSectionFile(GetHighResolutionSO2CrossSectionFile(), result);

        REQUIRE(result.m_crossSection.size() == 1402);
        REQUIRE(result.m_waveLength.size() == 1402);

        REQUIRE(result.m_waveLength.front() == Approx(238.9581));
        REQUIRE(result.m_waveLength.back() == Approx(395.0267));

        REQUIRE(result.m_crossSection.front() == Approx(3.754169E-20));
        REQUIRE(result.m_crossSection.back() == Approx(2.358910E-22));
    }

    // TODO: Implement reading a text file containing comments as well
    /* TEST_CASE("ReadCrossSectionFile with two column reference file from QDOAS", "[ReadCrossSectionFile][IntegrationTest]")
    {
        CCrossSectionData result;
        ReadCrossSectionFile(GetQDoasConvolvedSO2CrossSectionFile(), result);

        REQUIRE(result.m_crossSection.size() == 2048);
        REQUIRE(result.m_waveLength.size() == 2048);

        REQUIRE(result.m_waveLength.front() == Approx(278.463139));
        REQUIRE(result.m_waveLength.back() == Approx(425.206534));

        REQUIRE(result.m_crossSection.front() == Approx(6.94286262e-19));
        REQUIRE(result.m_crossSection.back() == Approx(0.0));
    }    */
}
