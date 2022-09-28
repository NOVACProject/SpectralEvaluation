#include "catch.hpp"
#include <SpectralEvaluation/Evaluation/EvaluationResult.h>
#include <SpectralEvaluation/Evaluation/DoasFit.h>
#include <sstream>

using namespace novac;

#pragma region Private helper methods

void CreateEvaluationResult(CEvaluationResult& result)
{
    result.m_chiSquare = 0.023;
    result.m_delta = 1.23;
    result.m_evaluationStatus = MARK_DELETED;
    result.m_polynomial[0] = 1.0;
    result.m_polynomial[1] = 2.0;
    result.m_polynomial[2] = 3.0;
    result.m_polynomial[3] = 4.0;
    result.m_polynomial[4] = 5.0;
    result.m_polynomial[5] = 6.0;
    result.m_stepNum = 56;

    CReferenceFitResult reference1{ 1.0, 0.1, 2.0, 0.2, 3.0, 0.3 };
    CReferenceFitResult reference2{ 21.0, 0.12, 22.0, 0.22, 23.0, 0.32 };
    CReferenceFitResult reference3{ 31.0, 0.13, 32.0, 0.23, 33.0, 0.33 };
    result.m_referenceResult.push_back(reference1);
    result.m_referenceResult.push_back(reference2);
    result.m_referenceResult.push_back(reference3);
}

DoasResult::ReferenceFitResult CreateReferenceFitResult(double multiplier)
{
    DoasResult::ReferenceFitResult result;
    result.column = multiplier;
    result.columnError = 2 * multiplier;
    result.shift = 3 * multiplier;
    result.shiftError = 4 * multiplier;
    result.squeeze = 5 * multiplier;
    result.squeezeError = 6 * multiplier;

    std::stringstream name;
    name << "name" << multiplier;
    result.name = name.str();

    return result;
}

#pragma endregion

TEST_CASE("CEvaluationResult - Copy constructor, copies relevant parameters", "[CEvaluationResult]")
{
    CEvaluationResult original;
    CreateEvaluationResult(original);

    // Act
    CEvaluationResult copy{ original };

    // Assert
    REQUIRE(original.m_chiSquare == copy.m_chiSquare);
    REQUIRE(original.m_delta == copy.m_delta);
    REQUIRE(original.m_evaluationStatus == copy.m_evaluationStatus);
    REQUIRE(original.m_polynomial[0] == copy.m_polynomial[0]);
    REQUIRE(original.m_polynomial[1] == copy.m_polynomial[1]);
    REQUIRE(original.m_polynomial[2] == copy.m_polynomial[2]);
    REQUIRE(original.m_polynomial[3] == copy.m_polynomial[3]);
    REQUIRE(original.m_polynomial[4] == copy.m_polynomial[4]);
    REQUIRE(original.m_polynomial[5] == copy.m_polynomial[5]);
    REQUIRE(original.m_stepNum == copy.m_stepNum);

    REQUIRE(original.m_referenceResult.size() == copy.m_referenceResult.size());
    for (size_t ii = 0; ii < original.m_referenceResult.size(); ++ii)
    {
        REQUIRE(original.m_referenceResult[ii].m_column == copy.m_referenceResult[ii].m_column);
        REQUIRE(original.m_referenceResult[ii].m_columnError == copy.m_referenceResult[ii].m_columnError);
        REQUIRE(original.m_referenceResult[ii].m_shift == copy.m_referenceResult[ii].m_shift);
        REQUIRE(original.m_referenceResult[ii].m_shiftError == copy.m_referenceResult[ii].m_shiftError);
        REQUIRE(original.m_referenceResult[ii].m_squeeze == copy.m_referenceResult[ii].m_squeeze);
        REQUIRE(original.m_referenceResult[ii].m_squeezeError == copy.m_referenceResult[ii].m_squeezeError);
        REQUIRE(original.m_referenceResult[ii].m_specieName == copy.m_referenceResult[ii].m_specieName);
    }
}

TEST_CASE("CEvaluationResult - Assignment operator, copies relevant parameters", "[CEvaluationResult]")
{
    CEvaluationResult original;
    CreateEvaluationResult(original);
    CEvaluationResult copy;

    // Act
    copy = original;

    // Assert
    REQUIRE(original.m_chiSquare == copy.m_chiSquare);
    REQUIRE(original.m_delta == copy.m_delta);
    REQUIRE(original.m_evaluationStatus == copy.m_evaluationStatus);
    REQUIRE(original.m_polynomial[0] == copy.m_polynomial[0]);
    REQUIRE(original.m_polynomial[1] == copy.m_polynomial[1]);
    REQUIRE(original.m_polynomial[2] == copy.m_polynomial[2]);
    REQUIRE(original.m_polynomial[3] == copy.m_polynomial[3]);
    REQUIRE(original.m_polynomial[4] == copy.m_polynomial[4]);
    REQUIRE(original.m_polynomial[5] == copy.m_polynomial[5]);
    REQUIRE(original.m_stepNum == copy.m_stepNum);

    REQUIRE(original.m_referenceResult.size() == copy.m_referenceResult.size());
    for (size_t ii = 0; ii < original.m_referenceResult.size(); ++ii)
    {
        REQUIRE(original.m_referenceResult[ii].m_column == copy.m_referenceResult[ii].m_column);
        REQUIRE(original.m_referenceResult[ii].m_columnError == copy.m_referenceResult[ii].m_columnError);
        REQUIRE(original.m_referenceResult[ii].m_shift == copy.m_referenceResult[ii].m_shift);
        REQUIRE(original.m_referenceResult[ii].m_shiftError == copy.m_referenceResult[ii].m_shiftError);
        REQUIRE(original.m_referenceResult[ii].m_squeeze == copy.m_referenceResult[ii].m_squeeze);
        REQUIRE(original.m_referenceResult[ii].m_squeezeError == copy.m_referenceResult[ii].m_squeezeError);
        REQUIRE(original.m_referenceResult[ii].m_specieName == copy.m_referenceResult[ii].m_specieName);
    }
}

TEST_CASE("CEvaluationResult - Copy constructor from DoasResult, copies relevant parameters", "[CEvaluationResult]")
{
    DoasResult original;
    original.chiSquare = 0.134;
    original.delta = 0.987;
    original.polynomialCoefficients = std::vector<double>{ 3, 4, 1, 2 };
    original.iterations = 98;
    original.referenceResult = std::vector<DoasResult::ReferenceFitResult>{
        CreateReferenceFitResult(1.0),
        CreateReferenceFitResult(2.0),
        CreateReferenceFitResult(3.0)
    };

    // Act
    CEvaluationResult copy = original;

    // Assert
    REQUIRE(original.chiSquare == copy.m_chiSquare);
    REQUIRE(original.delta == copy.m_delta);
    REQUIRE(original.polynomialCoefficients[0] == copy.m_polynomial[0]);
    REQUIRE(original.polynomialCoefficients[1] == copy.m_polynomial[1]);
    REQUIRE(original.polynomialCoefficients[2] == copy.m_polynomial[2]);
    REQUIRE(original.polynomialCoefficients[3] == copy.m_polynomial[3]);
    REQUIRE(0.0 == copy.m_polynomial[4]);
    REQUIRE(0.0 == copy.m_polynomial[5]);
    REQUIRE(original.iterations == copy.m_stepNum);

    REQUIRE(original.referenceResult.size() == copy.m_referenceResult.size());
    for (size_t ii = 0; ii < original.referenceResult.size(); ++ii)
    {
        REQUIRE(original.referenceResult[ii].column == copy.m_referenceResult[ii].m_column);
        REQUIRE(original.referenceResult[ii].columnError == copy.m_referenceResult[ii].m_columnError);
        REQUIRE(original.referenceResult[ii].shift == copy.m_referenceResult[ii].m_shift);
        REQUIRE(original.referenceResult[ii].shiftError == copy.m_referenceResult[ii].m_shiftError);
        REQUIRE(original.referenceResult[ii].squeeze == copy.m_referenceResult[ii].m_squeeze);
        REQUIRE(original.referenceResult[ii].squeezeError == copy.m_referenceResult[ii].m_squeezeError);
        REQUIRE(original.referenceResult[ii].name == copy.m_referenceResult[ii].m_specieName);
    }
}
