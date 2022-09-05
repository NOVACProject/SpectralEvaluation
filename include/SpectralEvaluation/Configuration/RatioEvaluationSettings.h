#pragma once

#include <string>

namespace Configuration
{
    // Defining how to select spectra for performing Ratio evaluations (i.e. calculating BrO/SO2 ratios in plumes).
    // These relates mostly to how to select spectra from the scan.
    struct RatioEvaluationSettings
    {
        // The minimum number of spectra which needs to be selected in the plume for the ratio calculation to be successful.
        int minNumberOfSpectraInPlume = 4;

        // The maximum number of spectra which are to be selected in the plume.
        // If the number of selected values are higher than this value then the spectrum to the left or to the right 
        //  which has the lowest SO2 column value is rejected until the total number of spectra selected equals this value.
        int maxNumberOfSpectraInPlume = 10;

        // The minimum number of spectra which needs to be averaged outside of the plume for the calculation to be successful.
        // int minNumberOfReferenceSpectra = 7;

        int numberOfSpectraOutsideOfPlume = 10;

        // The minimum (SO2) column for the selected spectra in the plume. Defaults to 40ppmm = 1e17 molec/cm2.
        double minInPlumeColumn = 1e17;

        // Lowest allowed angle to include, in degrees from zenith. Used to exclude spectra too close to the horizon.
        double minimumScanAngle = -75.0;

        // Highest allowed angle to include, in degrees from zenith. Used to exclude spectra too close to the horizon.
        double maximumScanAngle = +75.0;

        // The maximum saturation ratio (intensity / maximum intensity of spectrometer)
        double maxSaturationRatio = 0.88;

        // The minimum saturation ratio (intensity / maximum intensity of spectrometer)
        double minSaturationRatio = 0.12;

        // The minimum plume completeness for a scan to be considered
        double minimumPlumeCompleteness = 0.7;
    };
}