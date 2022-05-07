#include "catch.hpp"
#include <SpectralEvaluation/Calibration/Correspondence.h>
#include <SpectralEvaluation/Math/PolynomialFit.h>

namespace novac
{
    TEST_CASE("Correspondence - AllMeasuredPointsAreUnique", "[Calibration][WavelengthCalibration]")
    {
        Correspondence firstMeasuredFirstTheoreticalIdx = Correspondence(0, 0, 0.0);
        Correspondence firstMeasuredSecondTheoretical = Correspondence(0, 1, 0.0);
        Correspondence secondMeasuredSecondTheoretical = Correspondence(1, 1, 0.0);

        SECTION("Returns true for empty list")
        {
            std::vector<Correspondence> elements;

            bool result = AllMeasuredPointsAreUnique(elements);

            REQUIRE(true == result);
        }

        SECTION("Returns false if duplicate first measured idx")
        {
            std::vector<Correspondence> elements
            {
                firstMeasuredFirstTheoreticalIdx,
                secondMeasuredSecondTheoretical,
                firstMeasuredSecondTheoretical
            };

            bool result = AllMeasuredPointsAreUnique(elements);

            REQUIRE(false == result);
        }

        SECTION("Returns true if no duplicate measured idx")
        {
            std::vector<Correspondence> elements
            {
                firstMeasuredFirstTheoreticalIdx,
                secondMeasuredSecondTheoretical
            };

            bool result = AllMeasuredPointsAreUnique(elements);

            REQUIRE(true == result);
        }
    }

    TEST_CASE("Correspondence - AllTheoreticalPointsAreUnique", "[Calibration][WavelengthCalibration]")
    {
        Correspondence firstMeasuredFirstTheoreticalIdx = Correspondence(0, 0, 0.0);
        Correspondence secondMeasureFirstTheoreticaldIdx = Correspondence(1, 0, 0.0);
        Correspondence firstMeasuredSecondTheoretical = Correspondence(0, 1, 0.0);

        SECTION("Returns true for empty list")
        {
            std::vector<Correspondence> elements;

            bool result = AllTheoreticalPointsAreUnique(elements);

            REQUIRE(true == result);
        }

        SECTION("Returns false if duplicate first theoretical idx")
        {
            std::vector<Correspondence> elements
            {
                firstMeasuredFirstTheoreticalIdx,
                secondMeasureFirstTheoreticaldIdx,
                firstMeasuredSecondTheoretical
            };

            bool result = AllTheoreticalPointsAreUnique(elements);

            REQUIRE(false == result);
        }

        SECTION("Returns true if no duplicate theoretical idx")
        {
            std::vector<Correspondence> elements
            {
                firstMeasuredFirstTheoreticalIdx,
                firstMeasuredSecondTheoretical
            };

            bool result = AllTheoreticalPointsAreUnique(elements);

            REQUIRE(true == result);
        }
    }

    TEST_CASE("Correspondence - ListAllPossibleKeypointCombinations - too few elements", "[Calibration][WavelengthCalibration]")
    {
        Correspondence firstMeasuredFirstTheoreticalIdx = Correspondence(0, 0, 0.0);
        Correspondence secondMeasureFirstTheoreticaldIdx = Correspondence(1, 0, 0.0);
        Correspondence firstMeasuredSecondTheoretical = Correspondence(0, 1, 0.0);
        Correspondence secondMeasuredSecondTheoretical = Correspondence(1, 1, 0.0);

        SECTION("No possible correspondences - returns empty list")
        {
            std::vector<Correspondence> emptyList;

            auto result = ListAllPossibleKeypointCombinations(5, emptyList);

            REQUIRE(result.size() == 0);
        }

        SECTION("Requested combination of five correspondences but only four supplied - returns empty list")
        {
            std::vector<Correspondence> allCorrespondences
            {
                firstMeasuredFirstTheoreticalIdx,
                secondMeasureFirstTheoreticaldIdx,
                firstMeasuredSecondTheoretical,
                secondMeasuredSecondTheoretical
            };

            auto result = ListAllPossibleKeypointCombinations(5, allCorrespondences);

            REQUIRE(result.size() == 0);
        }
    }

    TEST_CASE("Correspondence - ListAllPossibleKeypointCombinations - all combinations possible", "[Calibration][WavelengthCalibration]")
    {
        // Notice that these does not contain any repetitions in measured nor theoretical index
        Correspondence firstMeasuredFirstTheoreticalIdx = Correspondence(0, 0, 0.0);
        Correspondence secondMeasureSecondTheoreticaldIdx = Correspondence(1, 1, 0.0);
        Correspondence thirdMeasuredThirdTheoretical = Correspondence(2, 2, 0.0);
        Correspondence fourthMeasuredFourthTheoretical = Correspondence(3, 3, 0.0);
        Correspondence fifthMeasuredFifthTheoretical = Correspondence(4, 4, 0.0);
        Correspondence sixthMeasuredSixthTheoretical = Correspondence(5, 5, 0.0);

        std::vector<Correspondence> allCorrespondences
        {
            firstMeasuredFirstTheoreticalIdx,
            secondMeasureSecondTheoreticaldIdx,
            thirdMeasuredThirdTheoretical,
            fourthMeasuredFourthTheoretical,
            fifthMeasuredFifthTheoretical,
            sixthMeasuredSixthTheoretical
        };

        SECTION("One element requested, returns each single correspondence")
        {
            auto result = ListAllPossibleKeypointCombinations(1, allCorrespondences);

            REQUIRE(result.size() == allCorrespondences.size());
            for (const auto& singleList : result)
            {
                REQUIRE(singleList.size() == 1);
            }
        }

        SECTION("Two element requested, returns every combination of two elements")
        {
            auto result = ListAllPossibleKeypointCombinations(2, allCorrespondences);

            REQUIRE(result.size() == 15); // actually = N!/(k! - N!)
            for (const auto& singleList : result)
            {
                REQUIRE(singleList.size() == 2);
                REQUIRE(AllMeasuredPointsAreUnique(singleList));
                REQUIRE(AllTheoreticalPointsAreUnique(singleList));
            }
        }

        SECTION("Three element requested, returns every combination of three elements")
        {
            auto result = ListAllPossibleKeypointCombinations(3, allCorrespondences);

            REQUIRE(result.size() == 20);
            for (const auto& singleList : result)
            {
                REQUIRE(singleList.size() == 3);
                REQUIRE(AllMeasuredPointsAreUnique(singleList));
                REQUIRE(AllTheoreticalPointsAreUnique(singleList));
            }
        }

        SECTION("Five element requested, returns every combination of five elements")
        {
            auto result = ListAllPossibleKeypointCombinations(5, allCorrespondences);

            REQUIRE(result.size() == 6);
            for (const auto& singleList : result)
            {
                REQUIRE(singleList.size() == 5);
                REQUIRE(AllMeasuredPointsAreUnique(singleList));
                REQUIRE(AllTheoreticalPointsAreUnique(singleList));
            }
        }
    }

    TEST_CASE("Correspondence - ListAllPossibleKeypointCombinations - few combinations possible", "[Calibration][WavelengthCalibration]")
    {
        Correspondence firstMeasuredFirstTheoreticalIdx = Correspondence(0, 0, 0.0);
        Correspondence secondMeasureFirstTheoreticaldIdx = Correspondence(1, 0, 0.0);
        Correspondence thirdMeasuredThirdTheoreticalIdx = Correspondence(2, 2, 0.0);
        Correspondence firstMeasuredFourthTheoreticalIdx = Correspondence(0, 3, 0.0);

        SECTION("No combination possible, returns no elements")
        {
            std::vector<Correspondence> allCorrespondences
            {
                firstMeasuredFirstTheoreticalIdx,
                secondMeasureFirstTheoreticaldIdx
            };

            auto result = ListAllPossibleKeypointCombinations(2, allCorrespondences);

            REQUIRE(result.size() == 0);
        }

        SECTION("One combination possible, returns one element")
        {
            std::vector<Correspondence> allCorrespondences
            {
                firstMeasuredFirstTheoreticalIdx,
                secondMeasureFirstTheoreticaldIdx,
                firstMeasuredFourthTheoreticalIdx
            };

            auto result = ListAllPossibleKeypointCombinations(2, allCorrespondences);

            REQUIRE(result.size() == 1);
            REQUIRE(result[0].size() == 2);
        }

        SECTION("Some combination possible, returns expected number")
        {
            // here it is possible to combine elements (0, 2), (1, 2), (1, 3), (2, 3)
            std::vector<Correspondence> allCorrespondences
            {
                firstMeasuredFirstTheoreticalIdx,
                secondMeasureFirstTheoreticaldIdx,
                thirdMeasuredThirdTheoreticalIdx,
                firstMeasuredFourthTheoreticalIdx
            };

            auto result = ListAllPossibleKeypointCombinations(2, allCorrespondences);

            REQUIRE(result.size() == 4);
            REQUIRE(result[0].size() == 2);
        }
    }

    TEST_CASE("Correspondence - SelectMaybeInliers - selects correct number and combination", "[Calibration][WavelengthCalibration]")
    {
        // This test actually uses a random number as input and hence have the risk of succeeding sometimes 
        //  even if the algorithm would be flawed... 
        std::random_device r;
        std::mt19937 rnd{ r() };

        // Notice that these does not contain any repetitions in measured nor theoretical index
        Correspondence firstMeasuredFirstTheoreticalIdx = Correspondence(0, 0, 0.0);
        Correspondence secondMeasureSecondTheoreticaldIdx = Correspondence(1, 1, 0.0);
        Correspondence thirdMeasuredThirdTheoretical = Correspondence(2, 2, 0.0);
        Correspondence fourthMeasuredFourthTheoretical = Correspondence(3, 3, 0.0);
        Correspondence fifthMeasuredFifthTheoretical = Correspondence(4, 4, 0.0);
        Correspondence sixthMeasuredSixthTheoretical = Correspondence(5, 5, 0.0);

        std::vector<Correspondence> allCorrespondences
        {
            firstMeasuredFirstTheoreticalIdx,
            secondMeasureSecondTheoreticaldIdx,
            thirdMeasuredThirdTheoretical,
            fourthMeasuredFourthTheoretical,
            fifthMeasuredFifthTheoretical,
            sixthMeasuredSixthTheoretical
        };

        SECTION("Two elements - returns two elements")
        {
            std::vector<Correspondence> result;
            SelectMaybeInliers(2, allCorrespondences, rnd, result);

            REQUIRE(2 == result.size());
            REQUIRE(AllTheoreticalPointsAreUnique(result));
            REQUIRE(AllMeasuredPointsAreUnique(result));
        }

        SECTION("Three elements - returns three elements")
        {
            std::vector<Correspondence> result;
            SelectMaybeInliers(3, allCorrespondences, rnd, result);

            REQUIRE(3 == result.size());
            REQUIRE(AllTheoreticalPointsAreUnique(result));
            REQUIRE(AllMeasuredPointsAreUnique(result));
        }

        SECTION("Five elements - returns five elements")
        {
            std::vector<Correspondence> result;
            SelectMaybeInliers(5, allCorrespondences, rnd, result);

            REQUIRE(5 == result.size());
            REQUIRE(AllTheoreticalPointsAreUnique(result));
            REQUIRE(AllMeasuredPointsAreUnique(result));
        }
    }

    TEST_CASE("Correspondence - FindCorrespondenceWithTheoreticalIdx - finds requested point", "[Calibration][WavelengthCalibration]")
    {
        Correspondence firstMeasuredFirstTheoreticalIdx = Correspondence(0, 0, 0.0);
        Correspondence secondMeasureFirstTheoreticaldIdx = Correspondence(1, 0, 0.0);
        Correspondence thirdMeasuredThirdTheoreticalIdx = Correspondence(2, 2, 0.0);
        Correspondence firstMeasuredFourthTheoreticalIdx = Correspondence(0, 3, 0.0);

        std::vector<Correspondence> allCorrespondences
        {
            firstMeasuredFirstTheoreticalIdx,
            secondMeasureFirstTheoreticaldIdx,
            thirdMeasuredThirdTheoreticalIdx,
            firstMeasuredFourthTheoreticalIdx
        };

        SECTION("First idx returns zero")
        {
            REQUIRE(0 == FindCorrespondenceWithTheoreticalIdx(allCorrespondences, 0));
        }

        SECTION("Second idx returns minus one")
        {
            REQUIRE(-1 == FindCorrespondenceWithTheoreticalIdx(allCorrespondences, 1));
        }

        SECTION("Third idx returns two")
        {
            REQUIRE(2 == FindCorrespondenceWithTheoreticalIdx(allCorrespondences, 2));
        }

        SECTION("Fourth idx returns three")
        {
            REQUIRE(3 == FindCorrespondenceWithTheoreticalIdx(allCorrespondences, 3));
        }
    }

    TEST_CASE("Correspondence - GetMeasuredValueSpan - returns expected range", "[Calibration][WavelengthCalibration]")
    {
        // Notice that these does not contain any repetitions in measured nor theoretical index
        Correspondence firstMeasuredFirstTheoreticalIdx = Correspondence(0, 0, 0.0);
        firstMeasuredFirstTheoreticalIdx.measuredValue = 0.0;
        Correspondence fourthMeasuredFourthTheoretical = Correspondence(3, 3, 0.0);
        fourthMeasuredFourthTheoretical.measuredValue = 3.0;
        Correspondence fifthMeasuredFifthTheoretical = Correspondence(4, 4, 0.0);
        fifthMeasuredFifthTheoretical.measuredValue = 4.0;
        Correspondence sixthMeasuredSixthTheoretical = Correspondence(5, 5, 0.0);
        sixthMeasuredSixthTheoretical.measuredValue = 5.0;

        SECTION("Index span one returns one")
        {
            std::vector<Correspondence> allCorrespondences{
                fifthMeasuredFifthTheoretical,
                fourthMeasuredFourthTheoretical
            };

            double result = GetMeasuredValueSpan(allCorrespondences);

            REQUIRE(result == 1.0);
        }

        SECTION("Index span five returns five")
        {
            std::vector<Correspondence> allCorrespondences{
                fifthMeasuredFifthTheoretical,
                sixthMeasuredSixthTheoretical, // largest measured = 5
                fourthMeasuredFourthTheoretical,
                firstMeasuredFirstTheoreticalIdx // smallest measured = 0
            };

            double result = GetMeasuredValueSpan(allCorrespondences);

            REQUIRE(result == 5.0);
        }
    }

    TEST_CASE("Correspondence - ArrangeByMeasuredKeypoint - correct arrangement", "[Calibration][WavelengthCalibration]")
    {
        Correspondence firstMeasuredFirstTheoreticalIdx = Correspondence(0, 0, 0.0);
        Correspondence secondMeasureFirstTheoreticaldIdx = Correspondence(1, 0, 0.0);
        Correspondence thirdMeasuredThirdTheoreticalIdx = Correspondence(2, 2, 0.0);
        Correspondence firstMeasuredFourthTheoreticalIdx = Correspondence(0, 3, 0.0);
        Correspondence firstMeasuredNinethTheoreticalIdx = Correspondence(0, 8, 0.0);

        SECTION("All unique measured keypoints")
        {
            std::vector<Correspondence> allCorrespondences
            {
                firstMeasuredFirstTheoreticalIdx,
                secondMeasureFirstTheoreticaldIdx,
                thirdMeasuredThirdTheoreticalIdx
            };

            auto result = ArrangeByMeasuredKeypoint(allCorrespondences);

            REQUIRE(result.size() == 3);
            REQUIRE(result[0].size() == 1);
            REQUIRE(result[1].size() == 1);
            REQUIRE(result[2].size() == 1);
        }

        SECTION("All same measured keypoint")
        {
            std::vector<Correspondence> allCorrespondences
            {
                firstMeasuredFirstTheoreticalIdx,
                firstMeasuredFourthTheoreticalIdx,
                firstMeasuredNinethTheoreticalIdx
            };

            auto result = ArrangeByMeasuredKeypoint(allCorrespondences);

            REQUIRE(result.size() == 1);
            REQUIRE(result[0].size() == 3);
        }

        SECTION("Some keypoints duplicates returns expected")
        {
            std::vector<Correspondence> allCorrespondences
            {
                firstMeasuredFirstTheoreticalIdx,
                secondMeasureFirstTheoreticaldIdx,
                thirdMeasuredThirdTheoreticalIdx,
                firstMeasuredFourthTheoreticalIdx,
                firstMeasuredNinethTheoreticalIdx
            };

            auto result = ArrangeByMeasuredKeypoint(allCorrespondences);

            REQUIRE(result.size() == 3);
            REQUIRE(result[0].size() == 3);
            REQUIRE(result[1].size() == 1);
            REQUIRE(result[2].size() == 1);
        }
    }

    TEST_CASE("Correspondence - ListInliers - list all elements in second list", "[Calibration][WavelengthCalibration]")
    {
        Correspondence firstMeasuredFirstTheoreticalIdx = Correspondence(0, 0, 0.0);
        Correspondence secondMeasureSecondTheoreticaldIdx = Correspondence(1, 1, 0.0);
        Correspondence thirdMeasuredThirdTheoretical = Correspondence(2, 2, 0.0);
        Correspondence fourthMeasuredFourthTheoretical = Correspondence(3, 3, 0.0);
        Correspondence fifthMeasuredFifthTheoretical = Correspondence(4, 4, 0.0);
        Correspondence sixthMeasuredSixthTheoretical = Correspondence(5, 5, 0.0);

        std::vector<Correspondence> allCorrespondences
        {
            firstMeasuredFirstTheoreticalIdx,
            secondMeasureSecondTheoreticaldIdx,
            thirdMeasuredThirdTheoretical,
            fourthMeasuredFourthTheoretical,
            fifthMeasuredFifthTheoretical,
            sixthMeasuredSixthTheoretical
        };

        SECTION("No selected element, all return values are false")
        {
            std::vector<Correspondence> selectedElements;

            auto result = ListInliers(selectedElements, allCorrespondences);

            REQUIRE(result.size() == allCorrespondences.size());
            for (bool value : result)
            {
                REQUIRE(false == value);
            }
        }

        SECTION("All selected elements, all return values are true")
        {
            auto result = ListInliers(allCorrespondences, allCorrespondences);

            REQUIRE(result.size() == allCorrespondences.size());
            for (bool value : result)
            {
                REQUIRE(true == value);
            }
        }

        SECTION("Some selected elements, returns expected list")
        {
            std::vector<Correspondence> selectedElements
            {
                allCorrespondences[1],
                allCorrespondences[4],
                allCorrespondences[5]
            };

            auto result = ListInliers(selectedElements, allCorrespondences);

            REQUIRE(result.size() == allCorrespondences.size());
            REQUIRE(false == result[0]);
            REQUIRE(true == result[1]);
            REQUIRE(false == result[2]);
            REQUIRE(false == result[3]);
            REQUIRE(true == result[4]);
            REQUIRE(true == result[5]);
        }
    }

    TEST_CASE("Correspondence - CountInliers - returns expected value for basic linear model", "[Calibration][WavelengthCalibration]")
    {
        SECTION("Model is increasing and fits correspondences exactly, returns all correspondences")
        {
            const std::vector<double> modelPolynomial{ 1.0, 2.0 }; // i.e. y = 1 + 2x

            // Create the set of correspondences
            std::vector<Correspondence> allCorrespondences;
            for (size_t ii = 0; ii < 5; ++ii)
            {
                Correspondence c;
                c.measuredIdx = ii;
                c.measuredValue = (double)ii;
                c.theoreticalIdx = 2 * ii;
                c.theoreticalValue = PolynomialValueAt(modelPolynomial, c.measuredValue);
                allCorrespondences.push_back(c);
            }
            auto correspondencesArrangedByMeasuredKeypoint = ArrangeByMeasuredKeypoint(allCorrespondences);

            std::vector<Correspondence> inliers;
            double meanError;
            bool isMonotonicallyIncreasing;
            size_t result = CountInliers(modelPolynomial, correspondencesArrangedByMeasuredKeypoint, 0.001, inliers, meanError, isMonotonicallyIncreasing);

            REQUIRE(result == allCorrespondences.size());
            REQUIRE(inliers.size() == allCorrespondences.size());
            REQUIRE(std::abs(meanError) < 0.0001);
            REQUIRE(isMonotonicallyIncreasing);
        }

        SECTION("Model is decreasing and fits correspondences exactly, returns all correspondences")
        {
            const std::vector<double> modelPolynomial{ 1.0, -2.0 }; // i.e. y = 1 - 2x

            // Create the set of correspondences
            std::vector<Correspondence> allCorrespondences;
            for (size_t ii = 0; ii < 5; ++ii)
            {
                Correspondence c;
                c.measuredIdx = ii;
                c.measuredValue = (double)ii;
                c.theoreticalIdx = 2 * ii;
                c.theoreticalValue = PolynomialValueAt(modelPolynomial, c.measuredValue);
                allCorrespondences.push_back(c);
            }
            auto correspondencesArrangedByMeasuredKeypoint = ArrangeByMeasuredKeypoint(allCorrespondences);

            std::vector<Correspondence> inliers;
            double meanError;
            bool isMonotonicallyIncreasing;
            size_t result = CountInliers(modelPolynomial, correspondencesArrangedByMeasuredKeypoint, 0.001, inliers, meanError, isMonotonicallyIncreasing);

            REQUIRE(result == allCorrespondences.size());
            REQUIRE(inliers.size() == allCorrespondences.size());
            REQUIRE(std::abs(meanError) < 0.0001);
            REQUIRE(false == isMonotonicallyIncreasing); // the model is not increasing.
        }

        SECTION("Model does not fit exactly, returns no correspondences if tight tolerance")
        {
            const std::vector<double> actualModelPolynomial{ 1.0, 2.0 }; // i.e. y = 1 + 2x
            const std::vector<double> suggestedModelPolynomial{ 2.0, 2.0 }; // i.e. y = 2 + 2x

            // Create the set of correspondences
            std::vector<Correspondence> allCorrespondences;
            for (size_t ii = 0; ii < 5; ++ii)
            {
                Correspondence c;
                c.measuredIdx = ii;
                c.measuredValue = (double)ii;
                c.theoreticalIdx = 2 * ii;
                c.theoreticalValue = PolynomialValueAt(actualModelPolynomial, c.measuredValue);
                allCorrespondences.push_back(c);
            }
            auto correspondencesArrangedByMeasuredKeypoint = ArrangeByMeasuredKeypoint(allCorrespondences);

            std::vector<Correspondence> inliers;
            double meanError;
            bool isMonotonicallyIncreasing;
            size_t result = CountInliers(suggestedModelPolynomial, correspondencesArrangedByMeasuredKeypoint, 0.001, inliers, meanError, isMonotonicallyIncreasing);

            REQUIRE(result == 0);
            REQUIRE(inliers.size() == 0);
        }

        SECTION("Model does not fit exactly, returns all correspondences if tolerance is larger than error")
        {
            // Notice that the two model polynomials here differ by exactly one everywhere. If the tolerance is larger than 1.0 then all correspondences are inliers.
            const std::vector<double> actualModelPolynomial{ 1.0, 2.0 };    // i.e. y = 1 + 2x
            const std::vector<double> suggestedModelPolynomial{ 2.0, 2.0 }; // i.e. y = 2 + 2x

            // Create the set of correspondences
            std::vector<Correspondence> allCorrespondences;
            for (size_t ii = 0; ii < 5; ++ii)
            {
                Correspondence c;
                c.measuredIdx = ii;
                c.measuredValue = (double)ii;
                c.theoreticalIdx = 2 * ii;
                c.theoreticalValue = PolynomialValueAt(actualModelPolynomial, c.measuredValue);
                allCorrespondences.push_back(c);
            }
            auto correspondencesArrangedByMeasuredKeypoint = ArrangeByMeasuredKeypoint(allCorrespondences);

            std::vector<Correspondence> inliers;
            double meanError;
            bool isMonotonicallyIncreasing;
            size_t result = CountInliers(suggestedModelPolynomial, correspondencesArrangedByMeasuredKeypoint, 1.5, inliers, meanError, isMonotonicallyIncreasing);

            REQUIRE(result == allCorrespondences.size());
            REQUIRE(inliers.size() == allCorrespondences.size());
            REQUIRE(std::abs(meanError - 1.0) < 0.0001);
            REQUIRE(isMonotonicallyIncreasing);
        }
    }
}
