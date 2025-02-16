#pragma once

#include <SpectralEvaluation/DialogControllers/WavelengthCalibrationController.h>

class NovacProgramWavelengthCalibrationController : public WavelengthCalibrationController
{
public:
    NovacProgramWavelengthCalibrationController(novac::ILogger& log)
        : WavelengthCalibrationController(log)
    {
        // NovacProgram will always add spectra together
        m_spectraAreAverages = false;
    }

    /** The full path to the spectrum to calibrate (should be a .pak file) */
    std::string m_inputSpectrumFile;

protected:

    virtual void ReadInput(novac::CSpectrum& measurement) override;

};