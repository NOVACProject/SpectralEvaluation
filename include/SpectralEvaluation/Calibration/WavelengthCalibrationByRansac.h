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

    /// <summary>
    /// A collection of the settings necessary to determine which keypoints
    /// will make up good correspondences between the measured and fraunhofer spectra.
    /// </summary>
    struct CorrespondenceSelectionSettings
    {
        /// <summary>
        /// The width, in pixels, around each keypoint which will be used to gauge the error in the correspondence.
        /// The default value of 20 is retrieved as 2x the average keypoint distance in tested spectra and should hence cover the entire width of a valley / peak.
        /// </summary>
        size_t pixelRegionSizeForCorrespondenceErrorMeasurement = 40;

        /// <summary>
        /// The maximum pixel error in the initial pixel-to-wavelength calibration guess.
        /// A smaller value allows for a faster and more precise calibrations whereas 
        /// a larger value makes it possible to converge also for larger initial calibration errors.
        /// Value here has been determined from experimentation.
        /// </summary>
        int maximumPixelDistanceForPossibleCorrespondence = 150;

        /// <summary>
        /// The relative number of correspondences to select out of the total.
        /// 0.2 corresponds to selecting the 20% correspondences with the lowest error.
        /// </summary>
        double percentageOfCorrespondencesToSelect = 0.2;

        /// <summary>
        /// The first pixel to include in the calibration routine. 
        /// Often do the signal in the spectra decline at short wavelengths and this is a means to disregard points with low intensity.
        /// </summary>
        size_t measuredPixelStart = 200;

        /// <summary>
        /// The last pixel to include in the calibration routine. 
        /// Often do the signal in the spectra decline at short wavelengths and this is a means to disregard points with low intensity.
        /// This must be larger than measuredPixelStart.
        /// </summary>
        size_t measuredPixelStop = 4095;
    };

    /// <summary>
    /// Measures the similarity between the two spectra at the two indices given by the correspondence.
    /// The similarity is measured as a sum of squared differences between the two 
    /// spectra in the region around the given pixels.
    /// A lower return value corresponds to a higher similarity
    /// <param name="corr">The Correspondence to measure the error of. This will be updated with the measured error.</param>
    /// <param name="measuredSpectrum">The measured spectrum of the correspondence.</param>
    /// <param name="theoreticalSpectrum">The theoretical spectrum of the correspondence.</param>
    /// </summary>
    void MeasureCorrespondenceError(novac::Correspondence& corr, const CSpectrum& measuredSpectrum, const CSpectrum& theoreticalSpectrum, const CorrespondenceSelectionSettings& settings);

    struct RansacWavelengthCalibrationSettings;
    /// <summary>
    /// This should be run as a preparatory step before the Ransac algorithm can be run.
    /// This generates the list of all reasonable correspondences between the measured and 
    /// fraunhofer spectra based on keypoints found in the two spectra.
    /// </summary>
    /// <param name="measuredKeypoints">The keypoints found in the measured spectrum.</param>
    /// <param name="measuredSpectrum">The measured spectrum itself.</param>
    /// <param name="fraunhoferKeypoints">The keypoints found in the fraunhofer spectrum.</param>
    /// <param name="fraunhoferSpectrum">The fraunhofer spectrum itself</param>
    /// <param name="settings">The settings for the following Ransac wavelength calibration</param>
    /// <param name="percentageOfCorrespondencesToSelect">A hard limit of how</param>
    /// <returns></returns>
    std::vector<novac::Correspondence> ListPossibleCorrespondences(
        const std::vector<novac::SpectrumDataPoint>& measuredKeypoints,
        const CSpectrum& measuredSpectrum,
        const std::vector<novac::SpectrumDataPoint>& fraunhoferKeypoints,
        const CSpectrum& fraunhoferSpectrum,
        const CorrespondenceSelectionSettings& correspondenceSettings);


    // TODO: Move
    double PolynomialValueAt(const std::vector<double>& coefficients, double x);

    // ------------- Wavelength calibration by Ransac  -------------

    struct RansacWavelengthCalibrationSettings
    {
        size_t modelPolynomialOrder = 3;

        int numberOfRansacIterations = 5000000;

        /// <summary>
        /// The number of correspondences to select in one iteration.
        /// Default (0) corresponds to (modelPolynomialOrder + 1)
        /// </summary>
        size_t sampleSize = 0;

        /// <summary>
        /// The number of pixels on the detector which generated this spectrum
        /// </summary>
        size_t detectorSize = 2048;

        /// <summary>
        /// how close a keypoint needs to be for it to be considered an inlier. In nm
        /// </summary>
        double inlierLimitInWavelength = 0.01;

        bool refine = false;

        /// <summary>
        /// The number of threads to divide the work up into.
        /// Special value: 0 corresponds to automatic.
        /// Default is 0
        /// </summary>
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

        /// <summary>
        /// Lists which of the incoming correspondences is an inlier.
        /// The number of true elements in this vector equals 'highestNumberOfInliers'
        /// </summary>
        std::vector<bool> correspondenceIsInlier;

        /** The smallest error in the model, by using the inliers */
        double smallestError = std::numeric_limits<double>::max();

        /** The largest pixel span over which there are inliers in this model */
        double largestPixelSpan = 0.0;

        /** The total number of possible correlations, the maximum number for 'highestNumberOfInliers' */
        size_t numberOfPossibleCorrelations = 0U;
    };

    /// <summary>
    /// RansacWavelengthCalibrationSetup is the setup of a calibration run
    ///     and contains all necessary elements to perform the calibration
    ///     by using a set of correspondences between a measured and an already calibrated spectrum.
    /// </summary>
    class RansacWavelengthCalibrationSetup
    {
    public:
        RansacWavelengthCalibrationSetup(RansacWavelengthCalibrationSettings calibrationSettings);

        /// <summary>
        /// This performs the actual calibration of a measured spectrum against a 
        ///   high resolution fraunhofer spectrum.
        /// </summary>
        /// <param name="possibleCorrespondences"></param>
        /// <param name="measuredKeypoints"></param>
        /// <param name="fraunhoferKeypoints"></param>
        /// <returns>The result of the calibration</returns>
        RansacWavelengthCalibrationResult DoWavelengthCalibration(const std::vector<Correspondence>& possibleCorrespondences);

    private:

        RansacWavelengthCalibrationSettings settings;

        /// <summary>
        /// Runs a part of the ransac calibrations. Used for splitting up the total workload on multiple threads, each performing a small part of the operations.
        /// </summary>
        /// <param name="possibleCorrespondences"></param>
        /// <returns>The best fitted model and the number of inliers.</returns>
        RansacWavelengthCalibrationResult RunRansacCalibrations(const std::vector<Correspondence>& possibleCorrespondences, const std::vector<std::vector<Correspondence>>& possibleCorrespondencesOrderedByMeasuredKeypoint, int numberOfIterations) const;

        /// <summary>
        /// Runs a basic calibration without random sampling, testing every possible combination.
        //  Only possible to use if the total number of combinations is low.
        /// </summary>
        /// <param name="possibleCorrespondences"></param>
        /// <returns>The best fitted model and the number of inliers.</returns>
        RansacWavelengthCalibrationResult RunDeterministicCalibration(const std::vector<Correspondence>& possibleCorrespondences, const std::vector<std::vector<Correspondence>>& possibleCorrespondencesOrderedByMeasuredKeypoint) const;

    };

}