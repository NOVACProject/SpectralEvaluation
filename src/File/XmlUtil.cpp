#include <SpectralEvaluation/File/XmlUtil.h>

#include <string.h>

namespace novac
{

std::string GetStartTag(const char* elementName)
{
    std::string tag = "<";
    tag.append(elementName);
    tag.append(">");
    return tag;
}
std::string GetStopTag(const char* elementName)
{
    std::string tag = "</";
    tag.append(elementName);
    tag.append(">");
    return tag;
}

int ParseXmlInteger(const char* xmlTag, const std::string& entireLine, int defaultValue)
{
    if (xmlTag == nullptr || strlen(xmlTag) == 0)
    {
        return defaultValue;
    }

    const std::string startTag = GetStartTag(xmlTag);
    const std::string stopTag = GetStopTag(xmlTag);

    const size_t firstIdx = entireLine.find(startTag);
    const size_t start = firstIdx + startTag.length();
    const size_t stop = entireLine.find(stopTag);
    if (stop > start && firstIdx != entireLine.npos && stop != entireLine.npos)
    {
        std::string valueStr = entireLine.substr(start, stop - start);
        return std::atoi(valueStr.c_str());
    }
    return defaultValue; // parse failure, return default value
}

std::string ParseXmlString(const char* xmlTag, const std::string& entireLine)
{
    if (xmlTag == nullptr || strlen(xmlTag) == 0)
    {
        return std::string();
    }

    const std::string startTag = GetStartTag(xmlTag);
    const std::string stopTag = GetStopTag(xmlTag);

    const size_t firstIdx = entireLine.find(startTag);
    const size_t start = firstIdx + startTag.length();
    const size_t stop = entireLine.find(stopTag);
    if (stop > start && firstIdx != entireLine.npos && stop != entireLine.npos)
    {
        return std::string(entireLine.c_str() + start, static_cast<int>(stop - start));
    }

    return std::string(); // parse failure, return empty string.
}

}
