#include <SpectralEvaluation/Calibration/InstrumentCalibration.h>
#include <SpectralEvaluation/Calibration/InstrumentLineShape.h>

namespace novac
{
InstrumentCalibration::~InstrumentCalibration()
{
    Clear();
}

void InstrumentCalibration::Clear()
{
    if (instrumentLineShapeParameter != nullptr)
    {
        delete instrumentLineShapeParameter;
        instrumentLineShapeParameter = nullptr;
    }
    pixelToWavelengthPolynomial.clear();
    pixelToWavelengthMapping.clear();
    instrumentLineShape.clear();
    instrumentLineShapeGrid.clear();
    instrumentLineShapeCenter = 0.0;
}

}
