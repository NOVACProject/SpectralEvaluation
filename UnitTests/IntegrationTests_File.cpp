#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include "catch.hpp"
#include "TestData.h"

namespace novac
{
    TEST_CASE("ReadCrossSectionFile with simple two column cross section file", "[ReadCrossSectionFile][IntegrationTest][File]")
    {
        CCrossSectionData result;
        ReadCrossSectionFile(TestData::GetHighResolutionSO2CrossSectionFile(), result);

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
        ReadCrossSectionFile(TestData::GetQDoasConvolvedSO2CrossSectionFile(), result);

        REQUIRE(result.m_crossSection.size() == 2048);
        REQUIRE(result.m_waveLength.size() == 2048);

        REQUIRE(result.m_waveLength.front() == Approx(278.463139));
        REQUIRE(result.m_waveLength.back() == Approx(425.206534));

        REQUIRE(result.m_crossSection.front() == Approx(6.94286262e-19));
        REQUIRE(result.m_crossSection.back() == Approx(0.0));
    }    */
}
