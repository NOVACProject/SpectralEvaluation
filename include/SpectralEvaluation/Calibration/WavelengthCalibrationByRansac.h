#pragma once

#include <vector>

// ---------------------------------------------------------------------------------------------------------------------
// - This header contains methods used to perform the wavelength calibration of a spectrometer using a ransac approach -
// ---------------------------------------------------------------------------------------------------------------------

namespace novac
{

/// <summary>
/// A Correspondence is an essential part in the ransac calibration routine, it
/// represents a connection between a point in the measured spectrum and a point
/// in the theoretical fraunhofer spectrum.
/// </summary>
struct Correspondence
{
    Correspondence()
        : measuredIdx(0), theoreticalIdx(0), error(0.0)
    { }

    Correspondence(size_t measured, size_t theoretical, double error)
        : measuredIdx(measured), theoreticalIdx(theoretical), error(error)
    { }

    /// <summary>
    /// The index of the keypoint in the measured spectrum (for reference)
    /// </summary>
    size_t measuredIdx = 0;

    /// <summary>
    /// The value of the keypoint in the measured spectrum (pixel, in the case of wavelength calibration)
    /// </summary>
    double measuredValue = 0.0;

    /// <summary>
    /// The index of the keypoint in the theoretical spectrum (for reference)
    /// </summary>
    size_t theoreticalIdx = 0;

    /// <summary>
    /// The value of the keypoint in the theoretical spectrum (wavelength, in the case of wavelength calibration)
    /// </summary>
    double theoreticalValue = 0.0;

    /// <summary>
    /// An error measure between the keypoints in the two spectra, lower is better
    /// </summary>
    double error = 0.0;
};


// TODO: Move
double PolynomialValueAt(const std::vector<double>& coefficients, double x);

struct RansacWavelengthCalibrationSettings
{
    size_t modelPolynomialOrder = 3;
    int numberOfRansacIterations = 500000;
    size_t sampleSize = 4; // the number of correspondences to select in one iteration
    double inlierLimitInWavelength = 0.5; // how close a keypoint needs to be for it to be considered an inlier. In nm
    int maximumPixelDistanceForPossibleCorrespondence = 150; // maximum pixel error in the initial clb file.
    bool refine = true;
};

struct RansacWavelengthCalibrationResult
{
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

    /** The smallest error in the model, by using the inliers */
    double smallestError = std::numeric_limits<double>::max();

    /** The total number of possible correlations, the maximum number for 'highestNumberOfInliers' */
    size_t numberOfPossibleCorrelations = 0U;
};

/// <summary>
/// RansacWavelengthCalibrationSetup is the setup of a calibration run
///     and contains all necessary elements to perform the calibration.
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

};

}