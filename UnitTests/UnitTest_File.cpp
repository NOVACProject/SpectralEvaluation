#include <SpectralEvaluation/File/File.h>
#include "catch.hpp"

namespace novac
{
    TEST_CASE("GetFileName Returns expected filename", "[GetFileName][File]")
    {
        SECTION("Windows path with regular file extension")
        {
            std::string input = "C:\\Windows\\System32\\wordfile.docx";
            std::string expected = "wordfile.docx";

            std::string result = GetFileName(input);

            REQUIRE(expected == result);
        }

        SECTION("Windows path without file extension")
        {
            std::string input = "C:\\Windows\\System32\\hosts";
            std::string expected = "hosts";

            std::string result = GetFileName(input);

            REQUIRE(expected == result);
        }

        SECTION("Windows path without file, returns exmpty string")
        {
            std::string input = "C:\\Windows\\System32\\";
            std::string expected = "";

            std::string result = GetFileName(input);

            REQUIRE(expected == result);
        }

        SECTION("Linux path without file extension")
        {
            std::string input = "/home/user/myfile.sh";
            std::string expected = "myfile.sh";

            std::string result = GetFileName(input);

            REQUIRE(expected == result);
        }
    }

    TEST_CASE("GetFileExtension Returns expected extension", "[GetFileName][File]")
    {
        SECTION("Windows path with regular file extension")
        {
            std::string input = "C:\\Windows\\System32\\wordfile.docx";
            std::string expected = ".docx";

            std::string result = GetFileExtension(input);

            REQUIRE(expected == result);
        }

        SECTION("Windows path without file extension")
        {
            std::string input = "C:\\Windows\\System32\\hosts";
            std::string expected = "";

            std::string result = GetFileExtension(input);

            REQUIRE(expected == result);
        }

        SECTION("Windows path without file, returns exmpty string")
        {
            std::string input = "C:\\Windows\\System32\\";
            std::string expected = "";

            std::string result = GetFileExtension(input);

            REQUIRE(expected == result);
        }

        SECTION("Linux path without file extension")
        {
            std::string input = "/home/user/myfile.sh";
            std::string expected = ".sh";

            std::string result = GetFileExtension(input);

            REQUIRE(expected == result);
        }
    }
}
