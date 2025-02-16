#include <SpectralEvaluation/File/File.h>
#include "catch.hpp"

namespace novac
{
TEST_CASE("GetFileName Returns expected filename", "[GetFileName][File]")
{
    SECTION("Windows path with regular file extension")
    {
        const std::string input = "C:\\Windows\\System32\\wordfile.docx";
        const std::string expected = "wordfile.docx";

        const std::string result = GetFileName(input);

        REQUIRE(expected == result);
    }

    SECTION("Windows path without file extension")
    {
        const std::string input = "C:\\Windows\\System32\\hosts";
        const std::string expected = "hosts";

        const std::string result = GetFileName(input);

        REQUIRE(expected == result);
    }

    SECTION("Windows path without file, returns exmpty string")
    {
        const std::string input = "C:\\Windows\\System32\\";
        const std::string expected = "";

        const std::string result = GetFileName(input);

        REQUIRE(expected == result);
    }

    SECTION("Linux path without file extension")
    {
        const std::string input = "/home/user/myfile.sh";
        const std::string expected = "myfile.sh";

        const std::string result = GetFileName(input);

        REQUIRE(expected == result);
    }
}

TEST_CASE("GetFileExtension Returns expected extension", "[GetFileName][File]")
{
    SECTION("Windows path with regular file extension")
    {
        const std::string input = "C:\\Windows\\System32\\wordfile.docx";
        const std::string expected = ".docx";

        const std::string result = GetFileExtension(input);

        REQUIRE(expected == result);
    }

    SECTION("Filename only, with lowercase file extension")
    {
        const std::string input = "wordfile.docx";
        const std::string expected = ".docx";

        const std::string result = GetFileExtension(input);

        REQUIRE(expected == result);
    }

    SECTION("Filename only, with uppercase file extension")
    {
        const std::string input = "WORDFILE.DOCX";
        const std::string expected = ".DOCX";

        const std::string result = GetFileExtension(input);

        REQUIRE(expected == result);
    }

    SECTION("Windows path without file extension")
    {
        const std::string input = "C:\\Windows\\System32\\hosts";
        const std::string expected = "";

        const std::string result = GetFileExtension(input);

        REQUIRE(expected == result);
    }

    SECTION("Windows path without file, returns empty string")
    {
        const std::string input = "C:\\Windows\\System32\\";
        const std::string expected = "";

        const std::string result = GetFileExtension(input);

        REQUIRE(expected == result);
    }

    SECTION("Linux path with file extension")
    {
        const std::string input = "/home/user/myfile.sh";
        const std::string expected = ".sh";

        const std::string result = GetFileExtension(input);

        REQUIRE(expected == result);
    }

    SECTION("Linux path without file extension")
    {
        const std::string input = "/home/user/myfile";
        const std::string expected = "";

        const std::string result = GetFileExtension(input);

        REQUIRE(expected == result);
    }
}
}
