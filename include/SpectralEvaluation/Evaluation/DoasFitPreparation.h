#pragma once

#include <vector>
#include <SpectralEvaluation/Evaluation/DoasFitEnumDeclarations.h>
#include <SpectralEvaluation/Math/IndexRange.h>

namespace novac
{
    class CSpectrum;

    /** The DoasFitPreparation class is a helper class for preparing spectra and references
        for a DOAS fit. */
    class DoasFitPreparation
    {
    public:
        /** Prepares the sky spectrum for inclusion into the DOAS fit.
            This should only be called with FIT_TYPE::FIT_HP_SUB or FIT_TYPE::FIT_POLY as the
            sky spectrum will only be included into the actual fit for these types.
            @throws std::invalid_argument if doasFitType is FIT_TYPE::FIT_HP_DIV. */
        static std::vector<double> PrepareSkySpectrum(const CSpectrum& skySpectrum, FIT_TYPE doasFitType);
        static std::vector<double> PrepareSkySpectrum(const CSpectrum& skySpectrum, FIT_TYPE doasFitType, const IndexRange& offsetRemovalRange);

        /** Prepares the measured spectrum for inclusion into the DOAS fit.
            This can be called for any type of fit, but does require that skySpectrum.m_length == measuredSpectrum.m_length */
        static std::vector<double> PrepareMeasuredSpectrum(const CSpectrum& measuredSpectrum, const CSpectrum& skySpectrum, FIT_TYPE doasFitType);
        static std::vector<double> PrepareMeasuredSpectrum(const CSpectrum& measuredSpectrum, const CSpectrum& skySpectrum, FIT_TYPE doasFitType, const IndexRange& offsetRemovalRange);

        /** Calculates a ring spectrum from the provided sky spectrum (which must have a wavelength calibration) 
            and prepares it for a DOAS fit of the given type */
        static std::vector<double> PrepareRingSpectrum(const CSpectrum& skySpectrum, FIT_TYPE doasFitType);

        /** Creates a reference which can be used as an intensity space polynomial (aka offset polynomial) 
            in a DOAS fit. This currently only supports the type FIT_TYPE::FIT_POLY and order == 0 of the polynomial. */
        static std::vector<double> PrepareIntensitySpacePolynomial(const CSpectrum& skySpectrum);

        /** Removes the electronic offset from the provided measured spectrum by calculating the average intensity
            in the pixel interval [startIndex, endIndex[ and then subtracting that value from all data points in the spectrum. */
        static void RemoveOffset(CSpectrum& spectrum, int startIndex = 50, int endIndex = 200);
        static void RemoveOffset(std::vector<double>& spectrum, int startIndex = 50, int endIndex = 200);

    };
}

