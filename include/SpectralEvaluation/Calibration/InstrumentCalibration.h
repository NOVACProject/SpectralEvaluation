#pragma once

#include <vector>

namespace novac
{
class ParametricInstrumentLineShape;

/** This class is an attempt to represent everything there is to know about
    the calibration of an instrument into one place.
    The instrument line shape is considered to be not wavelength dependent. */
class InstrumentCalibration
{
public:
    InstrumentCalibration() { };

    ~InstrumentCalibration();

    void Clear();

    // This object is not copyable due to the contained pointer.
    //  This could be implemented in the future if necessary
    InstrumentCalibration(const InstrumentCalibration&) = delete;
    InstrumentCalibration& operator=(const InstrumentCalibration&) = delete;

    /** The coefficients of the pixel to wavelength calibration.
        With the zero:th order coefficient first. */
    std::vector<double> pixelToWavelengthPolynomial;

    /** The mapping from pixels to wavelength.
        Length equals number of pixels on detector. */
    std::vector<double> pixelToWavelengthMapping;

    /** A sampling of the instrument line shape.
        These values should be normalized, i.e. fall in the range [0, 1] */
    std::vector<double> instrumentLineShape;

    /** A sampling grid for the instrument line shape.
        Must have equal length as 'instrumentLineShape' */
    std::vector<double> instrumentLineShapeGrid;

    /** The wavelength where the current instrument line shape is representative. */
    double instrumentLineShapeCenter = 0.0;

    /** A parametrization of the instrument line shape */
    ParametricInstrumentLineShape* instrumentLineShapeParameter = nullptr;
};

}