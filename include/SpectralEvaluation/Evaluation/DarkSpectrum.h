#pragma once
#include <string>

// This file contains common code used to extract a darks spectrum from a scan using different methods.

class CSpectrum;
namespace FileHandler
{
    class CScanFileHandler;
}
namespace Configuration
{
    struct CDarkSettings;
    struct CSkySettings;
}

namespace Evaluation
{
    /** Retrieves the dark spectrum to use for the provided spectrum given the settings from the user.
        @param scan The scan file to retrieve the spectrum from
        @param spec The spectrum for which the dark current spectrum should be retrieved.
        @param darkSettings The settings for how the dark current spectrum should be retrieved.
        @param dark Will on successful return be filled with the dark current spectrum.
        @param errorMessage Will be set to an error message if the result was not successful.
        @return true on success. */
    bool GetDark(FileHandler::CScanFileHandler& scan, const CSpectrum &spec, const Configuration::CDarkSettings& darkSettings, CSpectrum &dark, std::string& errorMessage);

    /** Models a dark spectrum to use for the provided spectrum from an offset and a dark-current spectrum which should exist
        in the provided scan.
        @param scan The scan file to retrieve the offset and dark-current spectra from.
        @param spec The spectrum for which the dark current spectrum should be modelled.
        @param darkSettings The settings for how the dark current spectrum should be retrieved.
        @param dark Will on successful return be filled with the dark current spectrum.
        @param errorMessage Will be set to an error message if the result was not successful.
        @return true on success. */
    bool ModelDarkSpectrum(FileHandler::CScanFileHandler& scan, const CSpectrum& spec, const Configuration::CDarkSettings& darkSettings, CSpectrum &dark, std::string& errorMessage);

    /** Reads the offset spectrum from the given scan and stores the result in offsetSpectrum.
        @return true on success. */
    bool GetOffsetSpectrum(FileHandler::CScanFileHandler& scan, const Configuration::CDarkSettings& darkSettings, CSpectrum& offsetSpectrum);

    /** Reads the dark-current spectrum from the given scan and stores the result in offsetSpectrum.
        @return true on success. */
    bool GetDarkCurrentSpectrum(FileHandler::CScanFileHandler& scan, const Configuration::CDarkSettings& darkSettings, CSpectrum& darkCurrent, bool& needsOffsetCorrection);

}
