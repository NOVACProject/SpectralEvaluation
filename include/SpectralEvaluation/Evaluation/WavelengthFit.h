#pragma once

#include "BasicMath.h"
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Evaluation/FitWindow.h>
#include "EvaluationResult.h"
#include "../Fit/ReferenceSpectrumFunction.h"

class CSpectrum;

namespace MathFit
{
class CStandardFit;
}

namespace Evaluation
{

// TODO: Move this to a separate file!
/** The CDoasFit performs a general DOAS fit using an already created set of
    references (with associated shift/squeeze options). */
class CDoasFit
{
public:
    CDoasFit() = default;

    int fitLow = 0;

    int fitHigh = 0;

};


/** The CWavelengthFit is a specialized DOAS fit used to perform a wavelength
    calibration of measured spectra using a Kurucz solar atlas + some additional references.
    The solar atlas and references are assumed to be well calibrated and this routine
    estimates a shift between the solar atlas and the measured spectrum.
    This is similar to the CEvaluationBase class. */
class CWavelengthFit : public CBasicMath
{
public:
    CWavelengthFit();


};
}