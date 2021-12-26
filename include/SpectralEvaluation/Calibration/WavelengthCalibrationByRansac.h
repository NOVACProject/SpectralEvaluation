#pragma once

#include <vector>
#include <limits>

// ---------------------------------------------------------------------------------------------------------------------
// - This header contains methods used to perform the wavelength calibration of a spectrometer using a ransac approach -
// ---------------------------------------------------------------------------------------------------------------------

namespace novac
{

    class CSpectrum;
    struct Correspondence;
    struct SpectrumDataPoint;


    // ------------- Keypoint selection and preparation -------------

    /** A collection of the settings necessary to determine which keypoints
        will make up good correspondences between the measured and fraunhofer spectra. */
    struct CorrespondenceSelectionSettings
    {
        /** The width, in pixels, around each keypoint which will be used to gauge the error in the correspondence.
            The default value of 20 is retrieved as 2x the average keypoint distance in tested spectra and should hence cover the entire width of a valley / peak. */
        size_t pixelRegionSizeForCorrespondenceErrorMeasurement = 40;

        /** The maximum pixel error in the initial pixel-to-wavelength calibration guess.
            A smaller value allows for a faster and more precise calibrations whereas
            a larger value makes it possible to converge also for larger initial calibration errors.
            Value here has been determined from experimentation. */
        int maximumPixelDistanceForPossibleCorrespondence = 500;

        /** The relative number of correspondences to select out of the total.
            0.2 corresponds to selecting the 20% correspondences with the lowest error. */
        double percentageOfCorrespondencesToSelect = 0.2;

        /** The first pixel to include in the calibration routine.
            Often do the signal in the spectra decline at short wavelengths and this is a means to disregard points with low intensity. */
        size_t measuredPixelStart = 200;

        /** The last pixel to include in the calibration routine.
            Often do the signal in the spectra decline at short wavelengths and this is a means to disregard points with low intensity.
            This must be larger than measuredPixelStart. */
        size_t measuredPixelStop = 4095;
    };

    /** Measures the similarity between the two spectra at the two indices given by the correspondence.
        The similarity is measured as a sum of squared differences between the two
        spectra in the region around the given pixels.
        A lower return value corresponds to a higher similarity
        @param corr The Correspondence to measure the error of. This will be updated with the measured error.
        @param measuredSpectrum The measured spectrum of the correspondence.
        @param theoreticalSpectrum The theoretical spectrum of the correspondence.*/
    void MeasureCorrespondenceError(novac::Correspondence& corr, const CSpectrum& measuredSpectrum, const CSpectrum& theoreticalSpectrum, const CorrespondenceSelectionSettings& settings);

    struct RansacWavelengthCalibrationSettings;

    /** This should be run as a preparatory step before the Ransac algorithm can be run.
        This generates the list of all reasonable correspondences between the measured and
        fraunhofer spectra based on keypoints found in the two spectra.
        @param measuredKeypoints The keypoints found in the measured spectrum.
        @param measuredSpectrum The measured spectrum itself.
        @param fraunhoferKeypoints The keypoints found in the fraunhofer spectrum.
        @param fraunhoferSpectrum The fraunhofer spectrum itself
        @param settings The settings for the following Ransac wavelength calibration
        @param percentageOfCorrespondencesToSelect A hard limit of how */
    std::vector<novac::Correspondence> ListPossibleCorrespondences(
        const std::vector<novac::SpectrumDataPoint>& measuredKeypoints,
        const CSpectrum& measuredSpectrum,
        const std::vector<novac::SpectrumDataPoint>& fraunhoferKeypoints,
        const CSpectrum& fraunhoferSpectrum,
        const CorrespondenceSelectionSettings& correspondenceSettings);

    // ------------- Wavelength calibration by Ransac  -------------

    struct RansacWavelengthCalibrationSettings
    {
        size_t modelPolynomialOrder = 3;

        int numberOfRansacIterations = 5000000;

        /** The number of correspondences to select in one iteration.
            Default (0) corresponds to (modelPolynomialOrder + 1) */
        size_t sampleSize = 0;

        /** The number of pixels on the detector which generated this spectrum */
        size_t detectorSize = 2048;

        /** how close a keypoint needs to be for it to be considered an inlier. In nm */
        double inlierLimitInWavelength = 0.01;

        bool refine = false;

        /** The number of threads to divide the work up into.
            Special value: 0 corresponds to automatic.
            Default is 0 */
        size_t numberOfThreads = 0;
    };

    struct RansacWavelengthCalibrationResult
    {
        RansacWavelengthCalibrationResult();
        RansacWavelengthCalibrationResult(size_t polynomialOrder);

        RansacWavelengthCalibrationResult(const RansacWavelengthCalibrationResult& other);
        RansacWavelengthCalibrationResult& operator=(const RansacWavelengthCalibrationResult& other);

        /** The best estimation of the pixel-to-wavelength mapping polynomial that we have.
            The coefficients make up a polynomial and are stored with the 0th order coefficient first. */
        std::vector<double> bestFittingModelCoefficients;

        /** The order of the 'bestFittingModelCoefficients' */
        size_t modelPolynomialOrder = 3;

        /** The number of inliers which were achieved with this model */
        size_t highestNumberOfInliers = 0U;

        /** Lists which of the incoming correspondences is an inlier.
            The number of true elements in this vector equals 'highestNumberOfInliers' */
        std::vector<bool> correspondenceIsInlier;

        /** The smallest error in the model, by using the inliers */
        double smallestError = std::numeric_limits<double>::max();

        /** The largest pixel span over which there are inliers in this model */
        double largestPixelSpan = 0.0;

        /** The total number of possible correlations, the maximum number for 'highestNumberOfInliers' */
        size_t numberOfPossibleCorrelations = 0U;
    };

    /** RansacWavelengthCalibrationSetup is the setup of a calibration run
            and contains all necessary elements to perform the calibration
            by using a set of correspondences between a measured and an already calibrated spectrum.*/
    class RansacWavelengthCalibrationSetup
    {
    public:
        RansacWavelengthCalibrationSetup(RansacWavelengthCalibrationSettings calibrationSettings);

        /** This performs the actual calibration of a measured spectrum against a
              high resolution fraunhofer spectrum.
            @return The result of the calibration */
        RansacWavelengthCalibrationResult DoWavelengthCalibration(const std::vector<Correspondence>& possibleCorrespondences);

    private:

        RansacWavelengthCalibrationSettings settings;

        /** Runs a part of the ransac calibrations. Used for splitting up the total workload on multiple threads, each performing a small part of the operations.
            @param possibleCorrespondences
            @return The best fitted model and the number of inliers. */
        RansacWavelengthCalibrationResult RunRansacCalibrations(const std::vector<Correspondence>& possibleCorrespondences, const std::vector<std::vector<Correspondence>>& possibleCorrespondencesOrderedByMeasuredKeypoint, int numberOfIterations) const;

        /**
            Runs a basic calibration without random sampling, testing every possible combination.
            Only possible to use if the total number of combinations is low.
            @param possibleCorrespondences
            @return The best fitted model and the number of inliers. */
        RansacWavelengthCalibrationResult RunDeterministicCalibration(const std::vector<Correspondence>& possibleCorrespondences, const std::vector<std::vector<Correspondence>>& possibleCorrespondencesOrderedByMeasuredKeypoint) const;

    };

}