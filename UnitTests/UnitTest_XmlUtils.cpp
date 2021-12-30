#include "catch.hpp"
#include <SpectralEvaluation/File/XmlUtil.h>


TEST_CASE("ParseXmlInteger", "[ParseXmlInteger][XmlUtils]")
{
    SECTION("Correct input - returns expected value")
    {
        const int expectedResult = 12345;
        const char* xmlTag = "SomeXmlTag";
        const std::string entireLine = "   <SomeXmlTag>12345</SomeXmlTag>  ";

        int result = novac::ParseXmlInteger(xmlTag, entireLine);

        REQUIRE(expectedResult == result);
    }

    SECTION("Input does not contain tag - returns default value")
    {
        const int defaultValue = 233;
        const int expectedResult = defaultValue;
        const char* xmlTag = "SomeXmlTag";
        const std::string entireLine = "   <AnotherXmlTag>12345</AnotherXmlTag>  ";

        int result = novac::ParseXmlInteger(xmlTag, entireLine, defaultValue);

        REQUIRE(expectedResult == result);
    }

    SECTION("Input does not contain end tag - returns default value")
    {
        const int defaultValue = 233;
        const int expectedResult = defaultValue;
        const char* xmlTag = "SomeXmlTag";
        const std::string entireLine = "   <SomeXmlTag>12345  ";

        int result = novac::ParseXmlInteger(xmlTag, entireLine, defaultValue);

        REQUIRE(expectedResult == result);
    } 

    SECTION("Input empty - returns default value")
    {
        const int defaultValue = 233;
        const int expectedResult = defaultValue;
        const char* xmlTag = "SomeXmlTag";
        const std::string entireLine = "";

        int result = novac::ParseXmlInteger(xmlTag, entireLine, defaultValue);

        REQUIRE(expectedResult == result);
    }
}

TEST_CASE("ParseXmlString", "[ParseXmlString][XmlUtils]")
{
    SECTION("Correct input - returns expected value")
    {
        const std::string expectedResult = "another file name";
        const char* xmlTag = "SomeXmlTag";
        const std::string entireLine = "   <SomeXmlTag>another file name</SomeXmlTag>  ";

        const std::string result = novac::ParseXmlString(xmlTag, entireLine);

        REQUIRE(expectedResult.compare(result) == 0);
    }

    SECTION("Input contains filename - returns expected value")
    {
        const std::string expectedResult = "C:\\References\\SO2_Bogumil(2003)_293K_239-395nm.xs";
        const char* xmlTag = "SomeXmlTag";
        const std::string entireLine = "   <SomeXmlTag>C:\\References\\SO2_Bogumil(2003)_293K_239-395nm.xs</SomeXmlTag>  ";

        const std::string result = novac::ParseXmlString(xmlTag, entireLine);

        REQUIRE(expectedResult.compare(result) == 0);
    }    

    SECTION("Input does not contain tag - returns empty string")
    {
        const std::string expectedResult = "another file name";
        const char* xmlTag = "SomeXmlTag";
        const std::string entireLine = "   <AnotherXmlTag>another file name</AnotherXmlTag>  ";

        const std::string result = novac::ParseXmlString(xmlTag, entireLine);

        REQUIRE(result.size() == 0);
    }

    SECTION("Input does not contain end tag - returns empty string")
    {
        const std::string expectedResult = "another file name";
        const char* xmlTag = "SomeXmlTag";
        const std::string entireLine = "   <SomeXmlTag>another file name  ";

        const std::string result = novac::ParseXmlString(xmlTag, entireLine);

        REQUIRE(result.size() == 0);
    } 

    SECTION("Input empty - returns empty string")
    {
        const char* xmlTag = "SomeXmlTag";
        const std::string entireLine = "";

        const std::string result = novac::ParseXmlString(xmlTag, entireLine);

        REQUIRE(result.size() == 0);
    }
}