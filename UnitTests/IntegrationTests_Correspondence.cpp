#include "catch.hpp"
#include "TestData.h"
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Evaluation/BasicMath.h>
#include <SpectralEvaluation/File/STDFile.h>
#include <SpectralEvaluation/File/TXTFile.h>
#include <SpectralEvaluation/Spectra/SpectrumUtils.h>
#include <SpectralEvaluation/Calibration/Correspondence.h>
#include <SpectralEvaluation/Calibration/WavelengthCalibrationByRansac.h>
#include <numeric>

namespace novac
{
    void RemoveBaseline(CSpectrum& spectrum); // in WavelengthCalibration.cpp

    static SpectrumDataPoint ClosestKeypointInPixels(const std::vector<SpectrumDataPoint>& allKeypoints, double expectedMeasuredPixel)
    {
        SpectrumDataPoint closestInPixels;
        double closestInPixelsValue = std::numeric_limits<double>::max();

        for (const auto& keypoint : allKeypoints)
        {
            const double diffPixel = std::abs(keypoint.pixel - expectedMeasuredPixel);

            if (diffPixel < closestInPixelsValue)
            {
                closestInPixelsValue = diffPixel;
                closestInPixels = keypoint;
            }
        }

        return closestInPixels;
    }

    /** Returns the point out of the points in 'fraunhoferKeypoints' which has the lowest correspondence error with the provided 'pointToSearchFor' */
    static SpectrumDataPoint GetCorrespondenceWithLowestError(const SpectrumDataPoint& pointToSearchFor, const CSpectrum& measuredSpectrum, const CSpectrum& fraunhoferSpectrum, const std::vector<SpectrumDataPoint>& fraunhoferKeypoints)
    {
        SpectrumDataPoint bestMatchingPointInFraunhoferSpectrum;
        double smallestError = std::numeric_limits<double>::max();
        novac::CorrespondenceSelectionSettings correspondenceSelectionSettings;
        for (const auto& pt : fraunhoferKeypoints)
        {
            const double error = novac::MeasureCorrespondenceError(measuredSpectrum, pointToSearchFor.pixel, fraunhoferSpectrum, pt.pixel, correspondenceSelectionSettings);

            if (error < smallestError)
            {
                smallestError = error;
                bestMatchingPointInFraunhoferSpectrum = pt;
            }
        }

        return bestMatchingPointInFraunhoferSpectrum;
    }

    static bool ContainsCorrespondence(const std::vector<Correspondence>& allCorrespondences, double expectedMeasuredPixel, double expectedWavelength, double pixelMargin = 2.0, double wavlengthMargin = 0.7)
    {
        Correspondence closestInPixels;
        double closestInPixelsValue = std::numeric_limits<double>::max();
        Correspondence closestInWavelength;
        double closestInWavelengthValue = std::numeric_limits<double>::max();

        for (const auto& corr : allCorrespondences)
        {
            const double diffPixel = std::abs(corr.measuredValue - expectedMeasuredPixel);
            const double diffWavelength = std::abs(corr.theoreticalValue - expectedWavelength);

            if (diffPixel < pixelMargin && diffWavelength < wavlengthMargin)
            {
                return true;
            }

            if (diffPixel < closestInPixelsValue)
            {
                closestInPixelsValue = diffPixel;
                closestInPixels = corr;
            }

            if (diffWavelength < closestInWavelengthValue)
            {
                closestInWavelengthValue = diffWavelength;
                closestInWavelength = corr;
            }
        }

        std::cout << "Failed to locate correspondence at pixel: " << expectedMeasuredPixel << " and wavelength: " << expectedWavelength << std::endl;
        std::cout << " Closest correspondence in pixels was at: " << closestInPixels.measuredValue << std::endl;
        std::cout << " Closest correspondence in wavelength was at: " << closestInWavelength.measuredValue << std::endl;

        return false;
    }

    TEST_CASE("ListPossibleCorrespondences in measured spectra from FLMS14634", "[WavelengthCalibrationByRansac][Correspondences][Flame]")
    {
        const double minimumPeakIntensityInMeasuredSpectrum = 0.02; // in the normalized units.
        const double minimumPeakIntensityInFraunhoferReference = 0.01; // in the normalized units.
        ::CBasicMath math;

        // Read in the measured spectrum
        CSpectrum measuredSpectrum;
        CSTDFile::ReadSpectrum(measuredSpectrum, TestData::GetMeasuredSpectrumName_FLMS14634());
        {
            CSpectrum darkSpectrum;
            CSTDFile::ReadSpectrum(darkSpectrum, TestData::GetDarkSpectrumName_FLMS14634());
            measuredSpectrum.Sub(darkSpectrum);
        }
        RemoveBaseline(measuredSpectrum);
        Normalize(measuredSpectrum);
        math.HighPassBinomial(measuredSpectrum.m_data, measuredSpectrum.m_length, 500);
        math.LowPassBinomial(measuredSpectrum.m_data, measuredSpectrum.m_length, 5);

        // Find the keypoints of the measured spectrum
        std::vector<SpectrumDataPoint> measuredKeypoints;
        novac::FindKeypointsInSpectrum(measuredSpectrum, minimumPeakIntensityInMeasuredSpectrum, measuredKeypoints);
        REQUIRE(measuredKeypoints.size() > 10);

        // Read in the Fraunhofer spectrum
        CSpectrum fraunhoferSpectrum;
        CTXTFile::ReadSpectrum(fraunhoferSpectrum, TestData::GetSyntheticFraunhoferSpectrumName_FLMS14634());
        Normalize(fraunhoferSpectrum);
        math.HighPassBinomial(fraunhoferSpectrum.m_data, fraunhoferSpectrum.m_length, 500);

        // Find the keypoints of the fraunhofer spectrum
        std::vector<SpectrumDataPoint> fraunhoferKeypoints;
        novac::FindKeypointsInSpectrum(fraunhoferSpectrum, minimumPeakIntensityInFraunhoferReference, fraunhoferKeypoints);
        REQUIRE(fraunhoferKeypoints.size() > 10);

        // Test with listing all the possible correspondences between the two meaured and Fraunhofer spectrum.
        {
            novac::CorrespondenceSelectionSettings correspondenceSelectionSettings;
            const auto allCorrespondences = novac::ListPossibleCorrespondences(measuredKeypoints, measuredSpectrum, fraunhoferKeypoints, fraunhoferSpectrum, correspondenceSelectionSettings);

            // Make sure that (some) of the very reasonable corresponences are included in the list.
            REQUIRE(ContainsCorrespondence(allCorrespondences, 1074, 359.8));
            REQUIRE(ContainsCorrespondence(allCorrespondences, 1585, 393.4));
            REQUIRE(ContainsCorrespondence(allCorrespondences, 1641, 396.8));
        }
    }

    TEST_CASE("MeasureCorrespondenceError returns lowest error for matching point in measured spectra from FLMS14634", "[WavelengthCalibrationByRansac][Correspondences][Flame]")
    {
        const double minimumPeakIntensityInMeasuredSpectrum = 0.02; // in the normalized units.
        const double minimumPeakIntensityInFraunhoferReference = 0.01; // in the normalized units.
        ::CBasicMath math;

        // Read in the measured spectrum
        CSpectrum measuredSpectrum;
        CSTDFile::ReadSpectrum(measuredSpectrum, TestData::GetMeasuredSpectrumName_FLMS14634());
        {
            CSpectrum darkSpectrum;
            CSTDFile::ReadSpectrum(darkSpectrum, TestData::GetDarkSpectrumName_FLMS14634());
            measuredSpectrum.Sub(darkSpectrum);
        }
        RemoveBaseline(measuredSpectrum);
        Normalize(measuredSpectrum);
        math.HighPassBinomial(measuredSpectrum.m_data, measuredSpectrum.m_length, 500);
        math.LowPassBinomial(measuredSpectrum.m_data, measuredSpectrum.m_length, 5);

        // Find the keypoints of the measured spectrum
        std::vector<SpectrumDataPoint> measuredKeypoints;
        novac::FindKeypointsInSpectrum(measuredSpectrum, minimumPeakIntensityInMeasuredSpectrum, measuredKeypoints);
        REQUIRE(measuredKeypoints.size() > 10);

        // Read in the Fraunhofer spectrum
        CSpectrum fraunhoferSpectrum;
        CTXTFile::ReadSpectrum(fraunhoferSpectrum, TestData::GetSyntheticFraunhoferSpectrumName_FLMS14634());
        Normalize(fraunhoferSpectrum);
        math.HighPassBinomial(fraunhoferSpectrum.m_data, fraunhoferSpectrum.m_length, 500);

        // Find the keypoints of the fraunhofer spectrum
        std::vector<SpectrumDataPoint> fraunhoferKeypoints;
        novac::FindKeypointsInSpectrum(fraunhoferSpectrum, minimumPeakIntensityInFraunhoferReference, fraunhoferKeypoints);
        REQUIRE(fraunhoferKeypoints.size() > 10);

        // Verify some points (which have been manually verified to agree)
        {
            const SpectrumDataPoint closestInPixels = ClosestKeypointInPixels(measuredKeypoints, 1585);
            const auto bestMatchingPointInFraunhoferSpectrum = GetCorrespondenceWithLowestError(closestInPixels, measuredSpectrum, fraunhoferSpectrum, fraunhoferKeypoints);
            REQUIRE(bestMatchingPointInFraunhoferSpectrum.wavelength == Approx(393).margin(0.5));
        }

        // Known (i.e. manually checked) that pixel 1585 in the measured spectrum corresponds to the wavelength 393nm in the Fraunhofer spectrum
        {
            const SpectrumDataPoint closestInPixels = ClosestKeypointInPixels(measuredKeypoints, 1640);
            const auto bestMatchingPointInFraunhoferSpectrum = GetCorrespondenceWithLowestError(closestInPixels, measuredSpectrum, fraunhoferSpectrum, fraunhoferKeypoints);
            REQUIRE(bestMatchingPointInFraunhoferSpectrum.wavelength == Approx(397).margin(0.5));
        }

        {
            const SpectrumDataPoint closestInPixels = ClosestKeypointInPixels(measuredKeypoints, 1074);
            const auto bestMatchingPointInFraunhoferSpectrum = GetCorrespondenceWithLowestError(closestInPixels, measuredSpectrum, fraunhoferSpectrum, fraunhoferKeypoints);
            REQUIRE(bestMatchingPointInFraunhoferSpectrum.wavelength == Approx(359.8).margin(0.5));
        }

        {
            const SpectrumDataPoint closestInPixels = ClosestKeypointInPixels(measuredKeypoints, 855);
            const auto bestMatchingPointInFraunhoferSpectrum = GetCorrespondenceWithLowestError(closestInPixels, measuredSpectrum, fraunhoferSpectrum, fraunhoferKeypoints);
            REQUIRE(bestMatchingPointInFraunhoferSpectrum.wavelength == Approx(344.2).margin(0.5));
        }

        {
            const SpectrumDataPoint closestInPixels = ClosestKeypointInPixels(measuredKeypoints, 925);
            const auto bestMatchingPointInFraunhoferSpectrum = GetCorrespondenceWithLowestError(closestInPixels, measuredSpectrum, fraunhoferSpectrum, fraunhoferKeypoints);
            REQUIRE(bestMatchingPointInFraunhoferSpectrum.wavelength == Approx(349.2).margin(0.5));
        }
    }
}
