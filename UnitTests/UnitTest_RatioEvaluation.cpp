#include "catch.hpp"
#include <Evaluation/RatioEvaluation.h>
#include <Evaluation/BasicScanEvaluationResult.h>

using namespace Evaluation;

TEST_CASE("IsSuitableScanForRatioEvaluation", "[RatioEvaluation][IsSuitableScanForRatioEvaluation]")
{
    RatioEvaluationSettings settings;
    BasicScanEvaluationResult scanResult;

    // TODO: Continue here...
    SECTION("Empty result - returns false")
    {
        REQUIRE(false == IsSuitableScanForRatioEvaluation(settings, scanResult));
    }

}