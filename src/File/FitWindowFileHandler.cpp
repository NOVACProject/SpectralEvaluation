#include <SpectralEvaluation/File/FitWindowFileHandler.h>
#include <SpectralEvaluation/StringUtils.h>
#include <SpectralEvaluation/File/File.h>

#include <rapidxml.hpp>

#include <iostream>

using namespace novac;

// TODO: Move
std::string FormatBoolean(bool value)
{
    return value ? "true" : "false";
}

std::string FormatEnum(novac::RING_CALCULATION_OPTION value)
{
    switch (value)
    {
    case novac::RING_CALCULATION_OPTION::CALCULATE_RING: return "calculate";
    case novac::RING_CALCULATION_OPTION::CALCULATE_RING_X2: return "calculatex2";
    default: return "none";
    }
}

template<class T> bool TagNameEqualsIgnoringCase(const T* element, const char* name)
{
    return EqualsIgnoringCase(element->name(), name);
}

template<class T> bool TagValueEqualsIgnoringCase(const T* element, const char* name)
{
    return EqualsIgnoringCase(element->value(), name);
}

bool ToBooleanOrDefault(const char* str, bool defaultValue)
{
    if (EqualsIgnoringCase(str, "true"))
    {
        return true;
    }
    else if (EqualsIgnoringCase(str, "false"))
    {
        return false;
    }
    return defaultValue;
}


bool ParseReference(rapidxml::xml_node<>* referenceNode, novac::CReferenceFile& reference)
{
    // Get the name of the reference itself.
    auto attr = referenceNode->first_attribute();
    while (attr != nullptr)
    {
        if (TagNameEqualsIgnoringCase(attr, "Name"))
        {
            reference.m_specieName = attr->value();
            break;
        }
        attr = attr->next_attribute();
    }

    auto childNode = referenceNode->first_node();
    while (childNode != nullptr)
    {
        if (TagNameEqualsIgnoringCase(childNode, "path"))
        {
            reference.m_path = childNode->value();
        }
        else if (TagNameEqualsIgnoringCase(childNode, "shiftOption"))
        {
            const int value = std::atoi(childNode->value());
            switch (value)
            {
            case 0: reference.m_shiftOption = novac::SHIFT_TYPE::SHIFT_FREE; break;
            case 1: reference.m_shiftOption = novac::SHIFT_TYPE::SHIFT_FIX; break;
            case 2: reference.m_shiftOption = novac::SHIFT_TYPE::SHIFT_LINK; break;
            case 3: reference.m_shiftOption = novac::SHIFT_TYPE::SHIFT_LIMIT; break;
            }
        }
        else if (TagNameEqualsIgnoringCase(childNode, "shiftValue"))
        {
            reference.m_shiftValue = std::atof(childNode->value());
        }
        else if (TagNameEqualsIgnoringCase(childNode, "squeezeOption"))
        {
            const int value = std::atoi(childNode->value());
            switch (value)
            {
            case 0: reference.m_squeezeOption = novac::SHIFT_TYPE::SHIFT_FREE; break;
            case 1: reference.m_squeezeOption = novac::SHIFT_TYPE::SHIFT_FIX; break;
            case 2: reference.m_squeezeOption = novac::SHIFT_TYPE::SHIFT_LINK; break;
            case 3: reference.m_squeezeOption = novac::SHIFT_TYPE::SHIFT_LIMIT; break;
            }
        }
        else if (TagNameEqualsIgnoringCase(childNode, "squeezeValue"))
        {
            reference.m_squeezeValue = std::atof(childNode->value());
        }
        else if (TagNameEqualsIgnoringCase(childNode, "columnOption"))
        {
            const int value = std::atoi(childNode->value());
            switch (value)
            {
            case 0: reference.m_columnOption = novac::SHIFT_TYPE::SHIFT_FREE; break;
            case 1: reference.m_columnOption = novac::SHIFT_TYPE::SHIFT_FIX; break;
            case 2: reference.m_columnOption = novac::SHIFT_TYPE::SHIFT_LINK; break;
            case 3: reference.m_columnOption = novac::SHIFT_TYPE::SHIFT_LIMIT; break;
            }
        }
        else if (TagNameEqualsIgnoringCase(childNode, "columnValue"))
        {
            reference.m_columnValue = std::atof(childNode->value());
        }

        childNode = childNode->next_sibling();
    }

    return true; // TODO: When to return false?
}

static bool ParseFitWindow(rapidxml::xml_node<>* fitWindowNode, novac::CFitWindow& window)
{
    // Get the name of the window itself.
    auto attr = fitWindowNode->first_attribute();
    while (attr != nullptr)
    {
        if (TagNameEqualsIgnoringCase(attr, "Name"))
        {
            window.name = attr->value();
            break;
        }
        attr = attr->next_attribute();
    }


    auto childNode = fitWindowNode->first_node();
    while (childNode != nullptr)
    {
        if (TagNameEqualsIgnoringCase(childNode, "fitLow"))
        {
            window.fitLow = std::atoi(childNode->value());
        }
        else if (TagNameEqualsIgnoringCase(childNode, "fitHigh"))
        {
            window.fitHigh = std::atoi(childNode->value());
        }
        else if (TagNameEqualsIgnoringCase(childNode, "polyOrder"))
        {
            window.polyOrder = std::atoi(childNode->value());
        }
        else if (TagNameEqualsIgnoringCase(childNode, "includeIntensitySpacePolyominal"))
        {
            window.includeIntensitySpacePolyominal = ToBooleanOrDefault(childNode->value(), false);
        }
        else if (TagNameEqualsIgnoringCase(childNode, "ringCalculation"))
        {
            if (TagValueEqualsIgnoringCase(childNode, "calculatex2"))
            {
                window.ringCalculation = RING_CALCULATION_OPTION::CALCULATE_RING_X2;
            }
            else if (TagValueEqualsIgnoringCase(childNode, "calculate"))
            {
                window.ringCalculation = RING_CALCULATION_OPTION::CALCULATE_RING;
            }
            else
            {
                window.ringCalculation = RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING;
            }
        }
        else if (TagNameEqualsIgnoringCase(childNode, "fitType"))
        {
            window.fitType = (FIT_TYPE)std::atoi(childNode->value());
        }
        else if (TagNameEqualsIgnoringCase(childNode, "channel"))
        {
            window.channel = std::atoi(childNode->value());
        }
        else if (TagNameEqualsIgnoringCase(childNode, "specLength"))
        {
            window.specLength = std::atoi(childNode->value());
        }
        else if (TagNameEqualsIgnoringCase(childNode, "fOptShift"))
        {
            window.findOptimalShift = std::atoi(childNode->value());
        }
        else if (TagNameEqualsIgnoringCase(childNode, "UV"))
        {
            const int UV = std::atoi(childNode->value());
            if (UV)
            {
                window.offsetRemovalRange = novac::CFitWindow::StandardUvOffsetRemovalRange();
            }
            else
            {
                window.offsetRemovalRange = novac::CFitWindow::StandardUSB2000OffsetRemovalRange();
            }
        }
        else if (TagNameEqualsIgnoringCase(childNode, "shiftSky"))
        {
            window.shiftSky = std::atoi(childNode->value());
        }
        else if (TagNameEqualsIgnoringCase(childNode, "skyShift"))
        {
            window.skyShift = std::atof(childNode->value());
        }
        else if (TagNameEqualsIgnoringCase(childNode, "skySqueeze"))
        {
            window.skySqueeze = std::atof(childNode->value());
        }
        else if (TagNameEqualsIgnoringCase(childNode, "interlaceStep"))
        {
            window.interlaceStep = std::atoi(childNode->value());
        }
        else if (TagNameEqualsIgnoringCase(childNode, "interlaced"))
        {
            window.interlaceStep = 1;
        }
        else if (TagNameEqualsIgnoringCase(childNode, "solarSpectrum"))
        {
            window.fraunhoferRef.m_path = childNode->value();
            window.fraunhoferRef.m_specieName = "SolarSpec";
        }
        else if (TagNameEqualsIgnoringCase(childNode, "ref"))
        {
            CReferenceFile reference;
            if (ParseReference(childNode, reference))
            {
                window.reference.push_back(reference);
            }
        }

        childNode = childNode->next_sibling();
    }

    return true; // when to return false??
}

std::vector<novac::CFitWindow> CFitWindowFileHandler::ReadFitWindowFile(const std::string& fileName)
{
    const bool caseSensitive = false;

    std::vector<novac::CFitWindow> allWindowsRead;
    if (!IsExistingFile(fileName))
    {
        return allWindowsRead;
    }

    std::string fileContents = ReadEntireFile(fileName);

    // Parse the contents of the file using rapidxml
    rapidxml::xml_document<> document;
    document.parse<0>(&fileContents[0]);

    auto fitWindowNode = document.first_node("fitWindow", 0, caseSensitive);
    if (nullptr == fitWindowNode)
    {
        return allWindowsRead;
    }

    while (fitWindowNode != nullptr)
    {
        novac::CFitWindow window;
        if (ParseFitWindow(fitWindowNode, window))
        {
            allWindowsRead.push_back(window);
        }

        fitWindowNode = fitWindowNode->next_sibling();
    }

    return allWindowsRead;
}

bool CFitWindowFileHandler::WriteFitWindow(const novac::CFitWindow& window, const std::string& fileName, bool overWrite)
{
    FILE* f = nullptr;
    std::string indent;

    // Open the file
    if (overWrite)
        f = fopen(fileName.c_str(), "w");
    else
        f = fopen(fileName.c_str(), "a+");

    // Check so that we could open the file.
    if (nullptr == f)
        return false;

    fprintf(f, "<fitWindow name=\"%s\">\n", window.name.c_str());
    indent = "\t";

    fprintf(f, "%s<fitLow>%d</fitLow>\n", indent.c_str(), window.fitLow);
    fprintf(f, "%s<fitHigh>%d</fitHigh>\n", indent.c_str(), window.fitHigh);
    fprintf(f, "%s<polyOrder>%d</polyOrder>\n", indent.c_str(), window.polyOrder);

    std::string entryValue = FormatBoolean(window.includeIntensitySpacePolyominal);
    fprintf(f, "%s<includeIntensitySpacePolyominal>%s</includeIntensitySpacePolyominal>\n", indent.c_str(), entryValue.c_str());

    entryValue = FormatEnum(window.ringCalculation);
    fprintf(f, "%s<ringCalculation>%s</ringCalculation>\n", indent.c_str(), entryValue.c_str());

    fprintf(f, "%s<fitType>%d</fitType>\n", indent.c_str(), static_cast<int>(window.fitType));

    fprintf(f, "%s<channel>%d</channel>\n", indent.c_str(), window.channel);
    fprintf(f, "%s<specLength>%d</specLength>\n", indent.c_str(), window.specLength);

    fprintf(f, "%s<fOptShift>%d</fOptShift>\n", indent.c_str(), window.findOptimalShift);
    fprintf(f, "%s<offsetFrom>%zd</offsetFrom>\n", indent.c_str(), window.offsetRemovalRange.from);
    fprintf(f, "%s<offsetTo>%zd</offsetTo>\n", indent.c_str(), window.offsetRemovalRange.to);
    fprintf(f, "%s<shiftSky>%d</shiftSky>\n", indent.c_str(), window.shiftSky);
    if (window.shiftSky == 2)
    {
        fprintf(f, "%s<skyShift>%lf</skyShift>\n", indent.c_str(), window.skyShift);
        fprintf(f, "%s<skySqueeze>%lf</skySqueeze>\n", indent.c_str(), window.skySqueeze);
    }
    fprintf(f, "%s<interlaceStep>%d</interlaceStep>\n", indent.c_str(), window.interlaceStep);

    if (window.fraunhoferRef.m_path.size() > 4)
    {
        fprintf(f, "%s<solarSpectrum>%s</solarSpectrum>\n", indent.c_str(), window.fraunhoferRef.m_path.c_str());
    }

    fprintf(f, "%s<nRef>%zd</nRef>\n", indent.c_str(), window.reference.size());

    for (size_t i = 0; i < window.reference.size(); ++i)
    {
        fprintf(f, "%s<ref name=\"%s\">\n", indent.c_str(), window.reference[i].m_specieName.c_str());
        fprintf(f, "%s\t<path>%s</path>\n", indent.c_str(), window.reference[i].m_path.c_str());

        fprintf(f, "%s\t<shiftOption>%d</shiftOption>\n", indent.c_str(), (int)window.reference[i].m_shiftOption);
        if (window.reference[i].m_shiftOption != novac::SHIFT_TYPE::SHIFT_FREE)
        {
            fprintf(f, "%s\t<shiftValue>%lf</shiftValue>\n", indent.c_str(), window.reference[i].m_shiftValue);
        }

        fprintf(f, "%s\t<squeezeOption>%d</squeezeOption>\n", indent.c_str(), (int)window.reference[i].m_squeezeOption);
        if (window.reference[i].m_squeezeOption != novac::SHIFT_TYPE::SHIFT_FREE)
        {
            fprintf(f, "%s\t<squeezeValue>%lf</squeezeValue>\n", indent.c_str(), window.reference[i].m_squeezeValue);
        }

        fprintf(f, "%s\t<columnOption>%d</columnOption>\n", indent.c_str(), (int)window.reference[i].m_columnOption);
        if (window.reference[i].m_columnOption != novac::SHIFT_TYPE::SHIFT_FREE)
        {
            fprintf(f, "%s\t<columnValue>%lf</columnValue>\n", indent.c_str(), window.reference[i].m_columnValue);
        }

        fprintf(f, "%s</ref>\n", indent.c_str());
    }

    fprintf(f, "</fitWindow>\n");

    // close the file again
    fclose(f);

    return true;
}

