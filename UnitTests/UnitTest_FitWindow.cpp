#include "catch.hpp"
#include <Evaluation/FitWindow.h>
#include <Evaluation/CrossSectionData.h>

TEST_CASE("FitWindow - Default constructor", "[CFitWindow]")
{
    Evaluation::CFitWindow sut;

    SECTION("No references are defined")
    {
        REQUIRE(sut.ref[0].m_data == nullptr);
        REQUIRE(sut.nRef == 0);
    }
}

Evaluation::CCrossSectionData* CreateCrossSection(int startValue)
{
    Evaluation::CCrossSectionData* obj = new Evaluation::CCrossSectionData();

    for (int k = 0; k < 100; ++k)
    {
        obj->m_waveLength.push_back(startValue + k);
        obj->m_crossSection.push_back(k * startValue);
    }

    return obj;
}

TEST_CASE("FitWindow - Copy constructor", "[CFitWindow]")
{
    Evaluation::CFitWindow original;

    SECTION("Original has no references, then no references are copied.")
    {
        REQUIRE(original.nRef == 0);
        REQUIRE(original.ref[0].m_data == nullptr);

        Evaluation::CReferenceFile ref1;
        ref1.m_specieName = "SO2";
        ref1.m_path = "C:/Novac/So2.txt";
        ref1.m_data = CreateCrossSection(1);
        
        Evaluation::CReferenceFile ref2;
        ref2.m_specieName = "O3";
        ref2.m_path = "C:/Novac/O3.txt";
        ref2.m_data = CreateCrossSection(2);

        original.ref[0] = ref1;
        original.ref[1] = ref2;
        original.nRef = 2;

        Evaluation::CFitWindow copy{original};

        REQUIRE(copy.nRef == 2);
        REQUIRE(copy.ref[0].m_data != nullptr);
        REQUIRE(copy.ref[1].m_data != nullptr);
    }

}