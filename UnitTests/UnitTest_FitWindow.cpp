#include "catch.hpp"
#include <SpectralEvaluation/Evaluation/FitWindow.h>
#include <SpectralEvaluation/VectorUtils.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>

using namespace novac;

TEST_CASE("FitWindow - Default constructor", "[CFitWindow]")
{
    CFitWindow sut;

    SECTION("No references are defined")
    {
        REQUIRE(sut.ref[0].m_data == nullptr);
        REQUIRE(sut.nRef == 0);
    }
}

CCrossSectionData* CreateCrossSection(int startValue)
{
    CCrossSectionData* obj = new CCrossSectionData();

    for (int k = 0; k < 100; ++k)
    {
        obj->m_waveLength.push_back(startValue + k);
        obj->m_crossSection.push_back(k * startValue);
    }

    return obj;
}

TEST_CASE("FitWindow - Copy constructor", "[CFitWindow]")
{
    CFitWindow original;

    SECTION("Original has no references, then no references are copied.")
    {
        REQUIRE(original.nRef == 0);
        REQUIRE(original.ref[0].m_data == nullptr);

        CReferenceFile ref1;
        ref1.m_specieName = "SO2";
        ref1.m_path = "C:/Novac/So2.txt";
        ref1.m_data.reset(CreateCrossSection(1));

        CReferenceFile ref2;
        ref2.m_specieName = "O3";
        ref2.m_path = "C:/Novac/O3.txt";
        ref2.m_data.reset(CreateCrossSection(2));

        original.ref[0] = ref1;
        original.ref[1] = ref2;
        original.nRef = 2;

        CFitWindow copy{ original };

        REQUIRE(copy.nRef == 2);
        REQUIRE(copy.ref[0].m_data != nullptr);
        REQUIRE(copy.ref[1].m_data != nullptr);
    }

}

TEST_CASE("FitWindow - AddAsReference", "[CFitWindow]")
{
    CFitWindow window;
    REQUIRE(window.nRef == 0);
    REQUIRE(window.ref[0].m_specieName.empty());
    std::vector<double> referenceData = GenerateVector(0.0, 1.0, window.specLength);
    const std::string referenceName = "ReferenceName";

    // Act
    AddAsReference(window, referenceData, referenceName);

    // Assert
    REQUIRE(window.nRef == 1);
    REQUIRE(window.ref[0].m_specieName == referenceName);
    REQUIRE((size_t)window.ref[0].m_data->GetSize() == referenceData.size());
    REQUIRE(window.ref[0].m_data->m_crossSection[10] == referenceData[10]);

    REQUIRE(window.ref[0].m_columnOption == SHIFT_TYPE::SHIFT_FREE);
    REQUIRE(window.ref[0].m_shiftOption == SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(window.ref[0].m_shiftValue == 0.0);
    REQUIRE(window.ref[0].m_squeezeOption == SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(window.ref[0].m_squeezeValue == 1.0);
}

TEST_CASE("FitWindow - AddAsSky", "[CFitWindow]")
{
    CFitWindow window;
    REQUIRE(window.nRef == 0);
    REQUIRE(window.ref[0].m_specieName.empty());
    std::vector<double> referenceData = GenerateVector(0.0, 1.0, window.specLength);

    // Act
    AddAsSky(window, referenceData);

    // Assert
    REQUIRE(window.nRef == 1);
    REQUIRE(window.ref[0].m_specieName == "sky");
    REQUIRE((size_t)window.ref[0].m_data->GetSize() == referenceData.size());
    REQUIRE(window.ref[0].m_data->m_crossSection[10] == referenceData[10]);

    REQUIRE(window.ref[0].m_columnOption == SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(window.ref[0].m_columnValue == -1.0);
    REQUIRE(window.ref[0].m_shiftOption == SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(window.ref[0].m_squeezeOption == SHIFT_TYPE::SHIFT_FIX);
}

TEST_CASE("FitWindow - AddAsSky - shift free", "[CFitWindow]")
{
    CFitWindow window;
    REQUIRE(window.nRef == 0);
    REQUIRE(window.ref[0].m_specieName.empty());
    std::vector<double> referenceData = GenerateVector(0.0, 1.0, window.specLength);

    // Act
    AddAsSky(window, referenceData, SHIFT_TYPE::SHIFT_FREE);

    // Assert
    REQUIRE(window.nRef == 1);
    REQUIRE(window.ref[0].m_specieName == "sky");
    REQUIRE((size_t)window.ref[0].m_data->GetSize() == referenceData.size());
    REQUIRE(window.ref[0].m_data->m_crossSection[10] == referenceData[10]);

    REQUIRE(window.ref[0].m_columnOption == SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(window.ref[0].m_columnValue == -1.0);
    REQUIRE(window.ref[0].m_shiftOption == SHIFT_TYPE::SHIFT_FREE);
    REQUIRE(window.ref[0].m_squeezeOption == SHIFT_TYPE::SHIFT_FIX);
}