#pragma once

#include <utility>
#include <memory>
#include <vector>
#include <string>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Spectra/SpectrumUtils.h>

// ---------------------------------------------------------------------------------------------------------------
// ----------- This header contains methods used to perform wavelength calibration of a spectrometer -------------
// ---------------------------------------------------------------------------------------------------------------

namespace novac
{

    class CSpectrum;
    struct SpectrumDataPoint;
    struct Correspondence;
    struct RansacWavelengthCalibrationResult;
    class FraunhoferSpectrumGeneration;
    class ICrossSectionSpectrumGenerator;
    class ParametricInstrumentLineShape;

    enum class InstrumentLineshapeEstimationOption
    {
        /** Do not make any estimate of the instrument line shape. */
        None = 0,

        /** Estimate the instrument line shape as an approximate Gaussian function,
        based on the general distance between peaks/valleys. */
        ApproximateGaussian,

        /** Estimate the instrument line shape as a Gaussian function,
        by minimizing the residual between the generated fraunhofer spectrum and the measured. */
        SuperGaussian
    };

    /** This is the result of running the wavelength calibration routine. */
    struct SpectrometerCalibrationResult
    {
        /** The final estimate for the pixel to wavelength mapping. */
        std::vector<double> pixelToWavelengthMapping;

        /** The coefficients of the pixel-to-wavelength mapping polynomial */
        std::vector<double> pixelToWavelengthMappingCoefficients;

        /** A basic error estimate for the pixel-to-wavelength mapping.
            This is the R2 of the polynomial fit to the included inliers. */
        double pixelToWavelengthMappingError;

        /** The maximum pixel difference between any two inliers into the
            pixel-to-wavelength mapping polynomial. A large value shows that the result
            is valid over a larger pixel range. */
        double pixelToWavelengthMappingPixelRange = 0.0;

        /** The total number of inliers into the pixel-to-wavelength mapping polynomial. */
        size_t pixelToWavelengthMappingInliers = 0;

        /** The estimation of the instrument line shape.
            This is only set if the instrument line shape is set to be estimated
            in the process, otherwise this is empty. */
        novac::CCrossSectionData estimatedInstrumentLineShape;

        /** An estimation of the error in the produced line shape.
            This is the chi2 of the DOAS fit where a synthetically generated Fraunhofer spectrum (and a Ring spectrum derived from it)
            were fitted to the measured spectrum.  */
        double estimatedInstrumentLineShapeError = 0.0;

        /** The shift applied while fitting the instrument line shape.
            This is the shift of the DOAS fit where a synthetically generated Fraunhofer spectrum (and a Ring spectrum derived from it)
            were fitted to the measured spectrum.
            A high value here indicates errors in the wavelength calibration. */
        double estimatedInstrumentLineShapeShift = 0.0;

        /** The Range of pixels over which the instrument line shape was estimated. */
        std::pair<size_t, size_t> estimatedInstrumentLineShapePixelRange;

        /** The parameterization of the estimated instrument line shape.
            This is only set if the instrument line shape is set to be estimated
            in the process, otherwise this is empty. */
        std::unique_ptr<novac::ParametricInstrumentLineShape> estimatedInstrumentLineShapeParameters;
    };

    struct WavelengthCalibrationSettings
    {
        /** The intial estimate for the pixel to wavelength mapping. */
        std::vector<double> initialPixelToWavelengthMapping;

        /** The initial estimate for the instrument line shape (measured or estimated). */
        novac::CCrossSectionData initialInstrumentLineShape;

        /** The path to the high resolved solar atlas to be used.
             Assumed that this has the wavelength unit of nm air. */
        std::string highResSolarAtlas;

        /** A set of high resolved cross sections to include into the generation
            of the fraunhofer spectrum. Each pair makes up a path to the cross section
            file on disk and a total column of the cross section to use.
            Only cross sections with a total column != 0 will be included. */
        std::vector<std::pair<std::string, double>> crossSections;

        /** A set of high resolved cross sections to include into the fit when
            fitting an instrument line shape. This should include references which are
            strongly absorbing in the wavelength range used to estimate the instrument line shape.
            NOTICE: Only the first cross section here is used so far.
            This is primarily used to remove the impact of Ozone on the fitted instrument line shape. */
        std::vector<std::string> crossSectionsForInstrumentLineShapeFitting;

        /** The option for how, and if, the instrument line shape should also be estimated
            during the wavelength calibration procedure. */
        InstrumentLineshapeEstimationOption estimateInstrumentLineShape = InstrumentLineshapeEstimationOption::None;

        /** The wavelength region in which the instrument line shape is estimated (if estimateInstrumentLineShape != None) */
        std::pair<double, double> estimateInstrumentLineShapeWavelengthRegion;
    };

    class WavelengthCalibrationFailureException : public std::exception
    {
    private:
        const char* const m_msg = "";

    public:
        WavelengthCalibrationFailureException(const char* msg) :
            m_msg(msg)
        {}

        const char* what() const noexcept override final { return m_msg; }
    };

    /**  Reads the pixel-to-wavelength mapping from a:
            1. .clb calibration data file - containing a single column of wavelength data OR
            2. .txt/.xs data file containing two columns of wavelength + cross section data
        If the file contains two columns then this will return the first. This cannot handle files with three or more columns of data. */
    std::vector<double> GetPixelToWavelengthMappingFromFile(const std::string& clbFile);

    /** Generates a pixel-to-wavelength mapping for each pixel on the detector given
        the wavelength calibration polynomial and the number of pixels. */
    std::vector<double> GetPixelToWavelengthMapping(const std::vector<double>& polynomialCoefficients, size_t detectorSize);

    /** Returns true if the provided polynomial can be a possible calibration for a spectrometer with the given detector size.
        This verifies that the provided polynomial is monotonically increasing in the pixel interval [0, numberOfPixels-1] */
    bool IsPossiblePixelToWavelengthCalibrationPolynomial(const std::vector<double>& candidatePolynomial, size_t numberOfPixels);

    /** Returns true if the provided vector is be a possible calibration for a spectrometer.
        This checks the mappings to make sure that they are monotonically increasing. */
    bool IsPossiblePixelToWavelengthCalibration(const std::vector<double>& pixelToWavelengthMapping);

    /** This is a helper structure used to extract the internal state of the pixel-to-wavelength calibration
        of a spectrometer from a measured mercury spectrum. */
    struct MercurySpectrumCalibrationState
    {
        /** This lists all the peaks found in the mercury spectrum (defined in pixels). */
        std::vector<SpectrumDataPoint> peaks;

        /** This is a list of peaks found in the spectrum but for which we didn't find a
            suitable emission line in the theoretical spectrum, or the line seems to be not fully resolved. */
        std::vector<SpectrumDataPoint> rejectedPeaks;

        /** If the Mercury calibration failed, then this will be filled in with the reason why */
        std::string errorMessage;
    };

    /**
     * @brief Performs a wavelength calibration of a spectrometer using a measured mercury spectrum.
     *  This will identify the peaks in the measured spectrum and fit a polynomial to them such that the
     *  pixel-to-wavelength mapping for the entire spectrum can be calculated.
     * @param measuredMercurySpectrum The measured spectrum
     * @param polynomialOrder The order of the polynomial to fit (in the range 1-3).
     * @param initialPixelToWavelength The initial guess for the pixel to wavelength calibration of the spectrum,
     *  This will be used to constrain the identification of the mercury lines.
     * @param result Will on successful return be filled with the resulting pixel-to-wavelength calibration.
     * @param state If not null, then this will be filled with information on how the calibration did perform.
     * @return True if the calibration was successful. */
    bool MercuryCalibration(
        const CSpectrum& measuredMercurySpectrum,
        int polynomialOrder,
        const std::vector<double>& initialPixelToWavelength,
        SpectrometerCalibrationResult& result,
        MercurySpectrumCalibrationState* state = nullptr);

    /** WavelengthCalibrationSetup is the setup of a calibration run
            and contains all necessary elements to perform the calibration. */
    class WavelengthCalibrationSetup
    {
    public:
        WavelengthCalibrationSetup(const WavelengthCalibrationSettings& calibrationSettings);

        /** This performs the actual calibration of a measured spectrum against a
              high resolution fraunhofer spectrum assuming that the provided instrument line shape
              is the correct line shape for the instrument.
            @throws std::invalid_argument if any of the incoming parameters is invalid.
            @throws WavelengthCalibrationFailureException if the calibration fails. */
        SpectrometerCalibrationResult DoWavelengthCalibration(const CSpectrum& measuredSpectrum);

        /** Simple structure used to save the internal state of the wavelength calibration. For inspection and debugging */
        struct SpectrumeterCalibrationState
        {
            std::unique_ptr<CSpectrum> measuredSpectrum;
            std::unique_ptr<CSpectrum> fraunhoferSpectrum;
            std::unique_ptr<CSpectrum> originalFraunhoferSpectrum;
            std::vector<novac::SpectrumDataPoint> measuredKeypoints;
            std::vector<novac::SpectrumDataPoint> fraunhoferKeypoints;
            std::vector<double> measuredSpectrumEnvelopePixels;
            std::vector<double> measuredSpectrumEnvelopeIntensities;
            std::vector<novac::Correspondence> allCorrespondences;
            std::vector<bool> correspondenceIsInlier;
        };

        /** Helper function which retrieves the state and result of the last call to 'DoWavelengthCalibration', for debugging. */
        const WavelengthCalibrationSetup::SpectrumeterCalibrationState& GetLastCalibrationSetup() { return this->calibrationState; }

    private:

        WavelengthCalibrationSettings settings;

        SpectrumeterCalibrationState calibrationState;

        /** Creates a wavelength to intensity spline of the measured spectrum using the current wavelength calibration result
            and correct the current generated Fraunhofer spectrum with it. This improves the accuracy of finding the spectrum peaks/valleys. */
        void UpdateFraunhoferSpectrumWithApparentSensitivity(novac::RansacWavelengthCalibrationResult& ransacResult);

        /** Creates an estimate of the instrument line shape of the measured spectrum as an approximate gaussian by judging the average distance between keypoints. */
        void EstimateInstrumentLineShapeAsApproximateGaussian(novac::SpectrometerCalibrationResult& result, novac::FraunhoferSpectrumGeneration& fraunhoferSetup);

        /** Creates an estimate of the instrument line shape of the measured spectrum by fitting a SuperGaussian to the measured spectrum. */
        void EstimateInstrumentLineShapeAsSuperGaussian(novac::SpectrometerCalibrationResult& result, novac::FraunhoferSpectrumGeneration& fraunhoferSetup, novac::ICrossSectionSpectrumGenerator* ozoneSetup = nullptr);

    };

}
