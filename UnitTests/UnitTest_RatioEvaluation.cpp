#include "catch.hpp"
#include <SpectralEvaluation/Evaluation/RatioEvaluation.h>
#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>

using namespace novac;

TEST_CASE("IsSuitableScanForRatioEvaluation", "[RatioEvaluation][IsSuitableScanForRatioEvaluation]")
{
    RatioEvaluationSettings settings;
    BasicScanEvaluationResult scanResult;
    CPlumeInScanProperty scanProperty;

    // TODO: Continue here...
    SECTION("Empty result - returns false")
    {
        REQUIRE(false == IsSuitableScanForRatioEvaluation(settings, scanResult, scanProperty));
    }

}