#pragma once

// This file contains basic Xml tag reading and writing used by the novac software

#include <string>

namespace novac
{

    /** Attempts to parse an integer from an XML tag. The tag and value must be on one single line.
    *   E.g. to parse the line:
    *       <MinimumColumn>5300</MinimumColumn>
    *   pass in xmlTag="MinimumColumn" and entireLine="   <MinimumColumn>5300</MinimumColumn>"
    *   @param xmlTag The expected tag, do not include the '<' or '>' of the Xml tag.
    *   @param entireLine The entire line to parse.
    *   @param defaultValue The default value to return if the parsing fails. */
    int ParseXmlInteger(const char* xmlTag, const std::string& entireLine, int defaultValue = 0);

    /** Attempts to parse a string from the contents inside one XML tag.
    *   The tag and value must be on one single line.
    *   E.g. to parse the line:
    *       <File>SO2_Bogumil(2003)_293K_239-395nm.xs</File>
    *   pass in xmlTag="File" and entireLine="    <File>SO2_Bogumil(2003)_293K_239-395nm.xs</File>"
    *   @param xmlTag The expected tag, do not include the '<' or '>' of the Xml tag.
    *   @param entireLine The entire line to parse.
    *   @return The string between the start tag and the end tag. Returns an empty string if the parse fails. */
    std::string ParseXmlString(const char* xmlTag, const std::string& entireLine);

}