#include <SpectralEvaluation/Calibration/InstrumentCalibration.h>
#include <SpectralEvaluation/Calibration/InstrumentLineShape.h>
#include <SpectralEvaluation/Calibration/WavelengthCalibration.h>
#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/VectorUtils.h>
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

// -------- Saving and reading instrument calibrations in .STD files  --------
TEST_CASE("InstrumentCalibration in extended std file", "[StdFile][InstrumentCalibration][IntegrationTest]")
{
    InstrumentCalibration originalCalibration;

    originalCalibration.pixelToWavelengthPolynomial = { 280.277915, 0.052285394, -2.3381128E-06, -9.55167214E-11 };
    originalCalibration.pixelToWavelengthMapping = GetPixelToWavelengthMapping(originalCalibration.pixelToWavelengthPolynomial, 2047);
    originalCalibration.instrumentLineShape =
    { 1.787112205e-16, 7.845422578e-14, 2.000832953e-11,
    0.000000002, 0.000000159, 0.000005429,
    0.000104902, 0.001202028, 0.008582179,
    0.040129658, 0.129208140, 0.301303033,
    0.535461911, 0.763515297, 0.920273384,
    0.988724318, 1, 0.988724318,
    0.920273384, 0.763515297, 0.535461911,
    0.301303033, 0.129208140, 0.040129658,
    0.008582179, 0.001202028, 0.000104902,
    0.000005429, 0.000000159, 0.000000002,
    2.000832953e-11, 7.845422578e-14, 1.787112205e-16 };  // 33 values, symmetricallly arranged around zero
    originalCalibration.instrumentLineShapeGrid.resize(originalCalibration.instrumentLineShape.size());
    size_t centerPixel = originalCalibration.instrumentLineShape.size() / 2;
    for (size_t ii = 0; ii < originalCalibration.instrumentLineShapeGrid.size(); ++ii)
    {
        originalCalibration.instrumentLineShapeGrid[ii] = originalCalibration.pixelToWavelengthMapping[410 + ii] - originalCalibration.pixelToWavelengthMapping[410 + centerPixel];
    }
    originalCalibration.instrumentLineShapeCenter = originalCalibration.pixelToWavelengthMapping[410 + centerPixel];

    const double originalW = 0.235;
    const double originalK = 2.954;
    {
        auto supergauss = new SuperGaussianLineShape();
        supergauss->k = originalK;
        supergauss->w = originalW;
        originalCalibration.instrumentLineShapeParameter = supergauss;
    }

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
        REQUIRE(readCalibration.pixelToWavelengthMapping.back() == Approx(originalCalibration.pixelToWavelengthMapping.back()));
    }

    SECTION("Correct pixel to wavelength polynomial  of spectrum read")
    {
        REQUIRE(readCalibration.pixelToWavelengthPolynomial.size() == originalCalibration.pixelToWavelengthPolynomial.size());
        REQUIRE(readCalibration.pixelToWavelengthPolynomial.front() == Approx(originalCalibration.pixelToWavelengthPolynomial.front()));
        REQUIRE(readCalibration.pixelToWavelengthPolynomial.back() == Approx(originalCalibration.pixelToWavelengthPolynomial.back()));
    }

    SECTION("Correct instrument line shape")
    {
        REQUIRE(readCalibration.instrumentLineShape.size() == originalCalibration.instrumentLineShape.size());
        REQUIRE(readCalibration.instrumentLineShape.front() == Approx(originalCalibration.instrumentLineShape.front()));
        REQUIRE(readCalibration.instrumentLineShape.back() == Approx(originalCalibration.instrumentLineShape.back()));
    }

    SECTION("Correct instrument line shape grid")
    {
        REQUIRE(readCalibration.instrumentLineShapeGrid.size() == originalCalibration.instrumentLineShapeGrid.size());
        REQUIRE(readCalibration.instrumentLineShapeGrid.front() == Approx(originalCalibration.instrumentLineShapeGrid.front()));
        REQUIRE(readCalibration.instrumentLineShapeGrid.back() == Approx(originalCalibration.instrumentLineShapeGrid.back()));
    }

    SECTION("Correct instrument line shape center")
    {
        REQUIRE(readCalibration.instrumentLineShapeCenter == Approx(originalCalibration.instrumentLineShapeCenter));
    }

    SECTION("Correct instrument line shape parametrization")
    {
        REQUIRE(readCalibration.instrumentLineShapeParameter != nullptr);
        REQUIRE(readCalibration.instrumentLineShapeParameter != originalCalibration.instrumentLineShapeParameter);

        const SuperGaussianLineShape* readParameters = dynamic_cast<SuperGaussianLineShape*>(readCalibration.instrumentLineShapeParameter);
        REQUIRE(readParameters != nullptr);
        REQUIRE(readParameters->w == Approx(originalW));
        REQUIRE(readParameters->k == Approx(originalK));
    }
}

TEST_CASE("InstrumentCalibration in extended std file - SuperGaussianInstrumentLineShape", "[StdFile][InstrumentCalibration][IntegrationTest]")
{
    InstrumentCalibration originalCalibration;

    originalCalibration.pixelToWavelengthPolynomial = { 280.277915, 0.052285394, -2.3381128E-06, -9.55167214E-11 };
    originalCalibration.pixelToWavelengthMapping = GetPixelToWavelengthMapping(originalCalibration.pixelToWavelengthPolynomial, 2047);

    SuperGaussianLineShape originalParameterizedLineShape;
    originalParameterizedLineShape.w = 0.324;
    originalParameterizedLineShape.k = 2.449;

    // Here the input instrument line shape does not have the same sampling density as the pixel-to-wavelength mapping above.
    {
        const auto originalSampledLineShape = SampleInstrumentLineShape(originalParameterizedLineShape);

        originalCalibration.instrumentLineShape = originalSampledLineShape.m_crossSection;
        originalCalibration.instrumentLineShapeGrid = originalSampledLineShape.m_waveLength;
    }
    originalCalibration.instrumentLineShapeCenter = originalCalibration.pixelToWavelengthMapping[410];

    {
        auto supergauss = new SuperGaussianLineShape();
        supergauss->k = originalParameterizedLineShape.k;
        supergauss->w = originalParameterizedLineShape.w;
        originalCalibration.instrumentLineShapeParameter = supergauss;
    }

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
        REQUIRE(readCalibration.pixelToWavelengthMapping.back() == Approx(originalCalibration.pixelToWavelengthMapping.back()));
    }

    SECTION("Correct pixel to wavelength polynomial  of spectrum read")
    {
        REQUIRE(readCalibration.pixelToWavelengthPolynomial.size() == originalCalibration.pixelToWavelengthPolynomial.size());
        REQUIRE(readCalibration.pixelToWavelengthPolynomial.front() == Approx(originalCalibration.pixelToWavelengthPolynomial.front()));
        REQUIRE(readCalibration.pixelToWavelengthPolynomial.back() == Approx(originalCalibration.pixelToWavelengthPolynomial.back()));
    }

    SECTION("Correct instrument line shape")
    {
        REQUIRE(readCalibration.instrumentLineShape.front() == Approx(originalCalibration.instrumentLineShape.front()).margin(0.001));
        REQUIRE(readCalibration.instrumentLineShape.back() == Approx(originalCalibration.instrumentLineShape.back()).margin(0.001));

        REQUIRE(Max(readCalibration.instrumentLineShape) == Approx(Max(originalCalibration.instrumentLineShape)).margin(0.001));
    }

    SECTION("Correct instrument line shape grid")
    {
        REQUIRE(readCalibration.instrumentLineShapeGrid.front() == Approx(originalCalibration.instrumentLineShapeGrid.front()).margin(0.1));
        REQUIRE(readCalibration.instrumentLineShapeGrid.back() == Approx(originalCalibration.instrumentLineShapeGrid.back()).margin(0.1));
    }

    SECTION("Correct instrument line shape center")
    {
        REQUIRE(readCalibration.instrumentLineShapeCenter == Approx(originalCalibration.instrumentLineShapeCenter));
    }

    SECTION("Correct instrument line shape parametrization")
    {
        REQUIRE(readCalibration.instrumentLineShapeParameter != nullptr);
        REQUIRE(readCalibration.instrumentLineShapeParameter != originalCalibration.instrumentLineShapeParameter);

        const SuperGaussianLineShape* readParameters = dynamic_cast<SuperGaussianLineShape*>(readCalibration.instrumentLineShapeParameter);
        REQUIRE(readParameters != nullptr);
        REQUIRE(readParameters->w == Approx(originalParameterizedLineShape.w));
        REQUIRE(readParameters->k == Approx(originalParameterizedLineShape.k));
    }
}

TEST_CASE("InstrumentCalibration in extended std file - SuperGaussianInstrumentLineShape - missing center information", "[StdFile][InstrumentCalibration][IntegrationTest]")
{
    InstrumentCalibration originalCalibration;

    originalCalibration.pixelToWavelengthPolynomial = { 280.277915, 0.052285394, -2.3381128E-06, -9.55167214E-11 };
    originalCalibration.pixelToWavelengthMapping = GetPixelToWavelengthMapping(originalCalibration.pixelToWavelengthPolynomial, 2047);

    SuperGaussianLineShape originalParameterizedLineShape;
    originalParameterizedLineShape.w = 0.324;
    originalParameterizedLineShape.k = 2.449;

    // Here the input instrument line shape does not have the same sampling density as the pixel-to-wavelength mapping above.
    {
        const auto originalSampledLineShape = SampleInstrumentLineShape(originalParameterizedLineShape);

        originalCalibration.instrumentLineShape = originalSampledLineShape.m_crossSection;
        originalCalibration.instrumentLineShapeGrid = originalSampledLineShape.m_waveLength;
    }

    originalCalibration.instrumentLineShapeCenter = 0.0; // i.e. we don't really know where in wavelength the instrument line shape is measured

    {
        auto supergauss = new SuperGaussianLineShape();
        supergauss->k = originalParameterizedLineShape.k;
        supergauss->w = originalParameterizedLineShape.w;
        originalCalibration.instrumentLineShapeParameter = supergauss;
    }

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
        REQUIRE(readCalibration.pixelToWavelengthMapping.back() == Approx(originalCalibration.pixelToWavelengthMapping.back()));
    }

    SECTION("Correct pixel to wavelength polynomial  of spectrum read")
    {
        REQUIRE(readCalibration.pixelToWavelengthPolynomial.size() == originalCalibration.pixelToWavelengthPolynomial.size());
        REQUIRE(readCalibration.pixelToWavelengthPolynomial.front() == Approx(originalCalibration.pixelToWavelengthPolynomial.front()));
        REQUIRE(readCalibration.pixelToWavelengthPolynomial.back() == Approx(originalCalibration.pixelToWavelengthPolynomial.back()));
    }

    SECTION("Correct instrument line shape")
    {
        REQUIRE(readCalibration.instrumentLineShape.front() == Approx(originalCalibration.instrumentLineShape.front()).margin(0.001));
        REQUIRE(readCalibration.instrumentLineShape.back() == Approx(originalCalibration.instrumentLineShape.back()).margin(0.001));

        REQUIRE(Max(readCalibration.instrumentLineShape) == Approx(Max(originalCalibration.instrumentLineShape)).margin(0.001));
    }

    SECTION("Correct instrument line shape grid")
    {
        REQUIRE(readCalibration.instrumentLineShapeGrid.front() == Approx(originalCalibration.instrumentLineShapeGrid.front()).margin(0.1));
        REQUIRE(readCalibration.instrumentLineShapeGrid.back() == Approx(originalCalibration.instrumentLineShapeGrid.back()).margin(0.1));
    }

    SECTION("Correct instrument line shape parametrization")
    {
        REQUIRE(readCalibration.instrumentLineShapeParameter != nullptr);
        REQUIRE(readCalibration.instrumentLineShapeParameter != originalCalibration.instrumentLineShapeParameter);

        const SuperGaussianLineShape* readParameters = dynamic_cast<SuperGaussianLineShape*>(readCalibration.instrumentLineShapeParameter);
        REQUIRE(readParameters != nullptr);
        REQUIRE(readParameters->w == Approx(originalParameterizedLineShape.w));
        REQUIRE(readParameters->k == Approx(originalParameterizedLineShape.k));
    }
}

}