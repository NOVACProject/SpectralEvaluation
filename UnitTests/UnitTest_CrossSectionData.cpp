#include "catch.hpp"
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>

using namespace novac;

TEST_CASE("Resample")
{
    novac::CCrossSectionData sut;

    SECTION("Same input and output - data is unchanged.")
    {
        sut.m_waveLength = {100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0};
        sut.m_crossSection = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};

        std::vector<double> result;
        novac::Resample(sut, 1.0, result);

        REQUIRE(result.size() == sut.m_crossSection.size());
        REQUIRE(result[0] == Approx(sut.m_crossSection[0]));
        REQUIRE(result[4] == Approx(sut.m_crossSection[4]));
        REQUIRE(result[9] == Approx(sut.m_crossSection[9]));
    }


}
