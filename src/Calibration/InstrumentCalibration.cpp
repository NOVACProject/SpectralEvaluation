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

    InstrumentCalibration::InstrumentCalibration(const InstrumentCalibration& other)
        : pixelToWavelengthPolynomial(other.pixelToWavelengthPolynomial),
        pixelToWavelengthMapping(other.pixelToWavelengthMapping),
        instrumentLineShape(other.instrumentLineShape),
        instrumentLineShapeGrid(other.instrumentLineShapeGrid),
        instrumentLineShapeCenter(other.instrumentLineShapeCenter)
    {
        if (other.instrumentLineShapeParameter != nullptr)
        {
            this->instrumentLineShapeParameter = other.instrumentLineShapeParameter->Clone();
        }
    }

    InstrumentCalibration& InstrumentCalibration::operator=(const InstrumentCalibration& other)
    {
        this->pixelToWavelengthPolynomial = other.pixelToWavelengthPolynomial;
        this->pixelToWavelengthMapping = other.pixelToWavelengthMapping;
        this->instrumentLineShape = other.instrumentLineShape;
        this->instrumentLineShapeGrid = other.instrumentLineShapeGrid;
        this->instrumentLineShapeCenter = other.instrumentLineShapeCenter;

        if (other.instrumentLineShapeParameter != nullptr)
        {
            this->instrumentLineShapeParameter = other.instrumentLineShapeParameter->Clone();
        }

        return *this;
    }

    InstrumentCalibration::InstrumentCalibration(InstrumentCalibration&& other)
    {
        this->pixelToWavelengthPolynomial = std::move(other.pixelToWavelengthPolynomial);
        this->pixelToWavelengthMapping = std::move(other.pixelToWavelengthMapping);
        this->instrumentLineShape = std::move(other.instrumentLineShape);
        this->instrumentLineShapeGrid = std::move(other.instrumentLineShapeGrid);
        this->instrumentLineShapeCenter = std::move(other.instrumentLineShapeCenter);

        this->instrumentLineShapeParameter = other.instrumentLineShapeParameter;
        other.instrumentLineShapeParameter = nullptr;
    }

    InstrumentCalibration& InstrumentCalibration::operator=(InstrumentCalibration&& other)
    {
        this->pixelToWavelengthPolynomial = std::move(other.pixelToWavelengthPolynomial);
        this->pixelToWavelengthMapping = std::move(other.pixelToWavelengthMapping);
        this->instrumentLineShape = std::move(other.instrumentLineShape);
        this->instrumentLineShapeGrid = std::move(other.instrumentLineShapeGrid);
        this->instrumentLineShapeCenter = std::move(other.instrumentLineShapeCenter);

        this->instrumentLineShapeParameter = other.instrumentLineShapeParameter;
        other.instrumentLineShapeParameter = nullptr;

        return *this;
    }

}
