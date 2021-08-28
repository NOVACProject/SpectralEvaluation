#pragma once

#include <vector>
#include <string>

namespace novac
{
class CSpectrum;
class CFitWindow;
class CReferenceFitResult;

class DoasFitException : public std::exception
{
private:
    const char* const m_msg = "";

public:
    DoasFitException(const char* msg) :
        m_msg(msg)
    {}

    const char* what() const noexcept override final { return m_msg; }
};

/// <summary>
/// Representation of the result of one DOAS evaluation.
/// </summary>
struct DoasResult
{
    /// <summary>
    /// The first pixel to include in the DOAS fit.
    /// </summary>
    int fitLow = 0;

    /// <summary>
    /// The first pixel after the end of the DOAS fit region, must be larger than m_fitLow.
    /// </summary>
    int fitHigh = 0;

    /// <summary>
    /// The chi2 of the fit, measure of the quality.
    /// </summary>
    double chiSquare = 0.0;

    /// <summary>
    /// The peak-to-peak amplitude of the residual, measure of the quality.
    /// </summary>
    double delta = 0.0;

    /// <summary>
    /// The number of iterations required to converge to this solution.
    /// </summary>
    int iterations = 0;

    /// <summary>
    /// The residual.
    /// Length equals the length of the fit region used.
    /// </summary>
    std::vector<double> residual;

    /// <summary>
    /// The valeus of the fitted polynomial.
    /// Length equals the length of the fit region used.
    /// </summary>
    std::vector<double> polynomialValues;

    /// <summary>
    /// The coefficients of the fitted polynomial.
    /// Saved with the 0th order coefficient first.
    /// </summary>
    std::vector<double> polynomialCoefficients;

    struct ReferenceFitResult
    {
        /// <summary>
        /// The resulting column
        /// </summary>
        double column = 0.0;

        /// <summary>
        /// The error in the resulting column
        /// </summary>
        double columnError = 0.0;

        /// <summary>
        /// The shift that was applied to the reference in the evaluation
        /// </summary>
        double shift = 0.0;

        /// <summary>
        /// The uncertainty in the applied shift
        /// </summary>
        double shiftError = 0.0;

        /// <summary>
        /// The squeeze that was applied to the reference in the evaluation
        /// </summary>
        double squeeze = 0.0;

        /// <summary>
        /// The uncertainty in the applied squeeze
        /// </summary>
        double squeezeError = 0.0;

        /// <summary>
        /// The name of the specie that the reference identifies
        /// </summary>
        std::string name = "";

        /// <summary>
        /// The scaled vales of the reference.
        /// This basically equals the reference's values multiplied by the column and adjusted for shift / squeeze.
        /// </summary>
        std::vector<double> scaledValues;
    };

    /// <summary>
    /// The evaluation result of each reference included.
    /// </summary>
    std::vector<ReferenceFitResult> referenceResult;
};

/** The DoasFit is a basic class for performing DOAS fits,
*   intended to simplify the rather complex setup of the EvaluationBase class.
*   May someday replace the EvaluationBase class altogether (if practical) */
class DoasFit
{
public:
    DoasFit();

    ~DoasFit();

    // Do not copy this object as it contains a member pointer.
    DoasFit(const DoasFit&) = delete;
    DoasFit& operator=(const DoasFit&) = delete;

    /** Sets up the parameters required to do a Doas fit */
    void Setup(const CFitWindow& setup);

    /** Runs the actual Doas fit.
    *   This assumes that the measuredData is already in OpticalDepth and will not do anything further processing with this.
    *   This also assumes that the CFitWindow contains the sky-spectrum / fraunhofer-reference-spectrum to use.
    *   @throws std::invalid_argument if Setup hasn't been called or if the input spectra are invalid.
    *   @throws DoasFitException if the fit itself failed for some reason. */
    void Run(const double* measuredData, size_t measuredLength, DoasResult& result);

private:

    /// <summary>
    /// The first pixel to include in the DOAS fit.
    /// </summary>
    int m_fitLow = 0;

    /// <summary>
    /// The first pixel after the end of the DOAS fit region, must be larger than m_fitLow.
    /// </summary>
    int m_fitHigh = 0;

    /// <summary>
    /// The order of the polynomial to include in the DOAS fit.
    /// </summary>
    int m_polynomialOrder = 3;

    /// <summary>
    /// The maximum number of steps in the DOAS evaluation.
    /// </summary>
    int m_maximumNumberOfSteps = 1000;

    /// <summary>
    /// Pointer to an internal structure holding each of the references which are to be included in the DOAS fit.
    /// This is setup when providing a FitWindow to use.
    /// </summary>
    void* m_referenceSetup = nullptr;
};

}