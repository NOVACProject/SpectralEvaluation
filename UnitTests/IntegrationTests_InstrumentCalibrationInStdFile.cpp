#include <SpectralEvaluation/Calibration/InstrumentCalibration.h>
#include <SpectralEvaluation/Calibration/WavelengthCalibration.h>
#include <SpectralEvaluation/File/File.h>
#include "catch.hpp"

namespace novac
{
std::string GetInstrumentCalibrationStdFileName()
{
#ifdef _MSC_VER
    return std::string("../TestData/InstrumentCalibration.std");
#else
    return std::string("TestData/InstrumentCalibration.std");
#endif // _MSC_VER
}

// -------- Doing a calibration from a mercury spectrum --------
TEST_CASE("InstrumentCalibration in extended std file", "[StdFile][InstrumentCalibration][IntegrationTest]")
{
    InstrumentCalibration originalCalibration;

    originalCalibration.pixelToWavelengthPolynomial = {280.277915, 0.052285394, -2.3381128E-06, -9.55167214E-11};
    originalCalibration.pixelToWavelengthMapping = GetPixelToWavelengthMapping(originalCalibration.pixelToWavelengthPolynomial, 2047);
    originalCalibration.instrumentLineShape = 
        {1.787112205e-16, 7.845422578e-14, 2.000832953e-11,
        0.000000002, 0.000000159, 0.000005429,
        0.000104902, 0.001202028, 0.008582179,
        0.040129658, 0.129208140, 0.301303033,
        0.535461911, 0.763515297, 0.920273384,
        0.988724318, 1, 0.988724318,
        0.920273384, 0.763515297, 0.535461911,
        0.301303033, 0.129208140, 0.040129658,
        0.008582179, 0.001202028, 0.000104902,
        0.000005429, 0.000000159, 0.000000002,
        2.000832953e-11, 7.845422578e-14, 1.787112205e-16};  // 33 values, symmetricallly arranged around zero
    originalCalibration.instrumentLineShapeGrid.resize(originalCalibration.instrumentLineShape.size());
    size_t centerPixel = originalCalibration.instrumentLineShape.size() / 2;
    for(size_t ii = 0; ii < originalCalibration.instrumentLineShapeGrid.size(); ++ii)
    {
        originalCalibration.instrumentLineShapeGrid[ii] = originalCalibration.pixelToWavelengthMapping[410 + ii] - originalCalibration.pixelToWavelengthMapping[410 + centerPixel];
    }
    originalCalibration.instrumentLineShapeCenter = originalCalibration.pixelToWavelengthMapping[410 + centerPixel];

    // Save the calibration above to file
    bool savedSuccessfully = SaveInstrumentCalibration(GetInstrumentCalibrationStdFileName(), originalCalibration);
    REQUIRE(savedSuccessfully);

    // Now read in the calibration again from file and make sure that we get the original contents back 
    InstrumentCalibration readCalibration;
    bool readSuccessfully = ReadInstrumentCalibration(GetInstrumentCalibrationStdFileName(), readCalibration);
    REQUIRE(readSuccessfully);

    SECTION("Correct pixel to wavelength mapping of spectrum read")
    {
        REQUIRE(readCalibration.pixelToWavelengthMapping.size() == originalCalibration.pixelToWavelengthMapping.size());
        REQUIRE(readCalibration.pixelToWavelengthMapping.front() == Approx(originalCalibration.pixelToWavelengthMapping.front()));
        REQUIRE(readCalibration.pixelToWavelengthMapping.back() ==  Approx(originalCalibration.pixelToWavelengthMapping.back()));
    }

    SECTION("Correct instrument line shape")
    {
        REQUIRE(readCalibration.instrumentLineShape.size() == originalCalibration.instrumentLineShape.size());
        REQUIRE(readCalibration.instrumentLineShape.front() ==  Approx(originalCalibration.instrumentLineShape.front()));
        REQUIRE(readCalibration.instrumentLineShape.back() ==  Approx(originalCalibration.instrumentLineShape.back()));
    }

    SECTION("Correct instrument line shape grid")
    {
        REQUIRE(readCalibration.instrumentLineShapeGrid.size() == originalCalibration.instrumentLineShapeGrid.size());
        REQUIRE(readCalibration.instrumentLineShapeGrid.front() ==  Approx(originalCalibration.instrumentLineShapeGrid.front()));
        REQUIRE(readCalibration.instrumentLineShapeGrid.back() ==  Approx(originalCalibration.instrumentLineShapeGrid.back()));
    }

}


}