#include <SpectralEvaluation/Calibration/StandardCrossSectionSetup.h>
#include <SpectralEvaluation/File/XmlUtil.h>
#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/StringUtils.h>

#include <fstream>
#include <sstream>
#include <string.h>

namespace novac
{

    StandardCrossSectionSetup::StandardCrossSectionSetup(const std::string& directory)
    {
        // check if there is a set of standard references in the expected format in the provided directory.

        try
        {
            StandardReference currentReference;

            std::string xmlFileLocation = directory;
            if (xmlFileLocation.size() > 0 && xmlFileLocation.back() != '/' && xmlFileLocation.back() != '\\')
            {
                xmlFileLocation.append("/");
            }
            xmlFileLocation.append("StandardCrossSections/");

            std::string xmlFileName = xmlFileLocation;
            xmlFileName.append("StandardCrossSections.xml");

            // Super basic xml parsing, without using any additional libraries.
            std::ifstream file(xmlFileName, std::ios::in);
            std::string line;
            while (std::getline(file, line))
            {
                if (line.find("/CrossSection") != std::string::npos)
                {
                    // Add the currently read reference to the end, if it points to a valid location.
                    if (IsExistingFile(currentReference.fileName))
                    {
                        m_standardReferences.push_back(currentReference);
                    }
                }
                else if (line.find("CrossSection") != std::string::npos)
                {
                    // Start the reading over with a new reference
                    currentReference = StandardReference();
                }
                else if (line.find("/FraunhoferSpectrum") != std::string::npos)
                {
                    // Add the currently read reference as the Fraunhofer reference to use, if it points to a valid location.
                    if (IsExistingFile(currentReference.fileName))
                    {
                        m_fraunhoferReference.fileName = currentReference.fileName;
                        m_fraunhoferReference.isVacuum = currentReference.isVacuum;
                    }
                }
                else if (line.find("FraunhoferSpectrum") != std::string::npos)
                {
                    // Start the reading over with a new reference
                    currentReference = StandardReference();
                }
                else if (line.find("File") != std::string::npos)
                {
                    std::string fileName = ParseXmlString("File", line);
                    if (fileName.size() > 0)
                    {
                        currentReference.fileName = xmlFileLocation;
                        currentReference.fileName.append(fileName);
                    }
                }
                else if (line.find("Name") != std::string::npos)
                {
                    currentReference.specieName = ParseXmlString("Name", line);
                }
                else if (line.find("Medium") != std::string::npos)
                {
                    std::string mediumName = ParseXmlString("Medium", line);

                    if (mediumName.length() > 0 && EqualsIgnoringCase(mediumName.c_str(), "air", 3))
                    {
                        currentReference.isVacuum = false;
                    }
                    else
                    {
                        currentReference.isVacuum = true;
                    }
                }
            }
        }
        catch (std::exception&)
        {
        }
    }

    std::vector<std::string> StandardCrossSectionSetup::ListReferences() const
    {
        std::vector<std::string> result;

        for (const auto& reference : m_standardReferences)
        {
            std::string name = reference.specieName;
            if (reference.isVacuum)
            {
                name.append(" (vac)");
            }
            result.push_back(name);
        }
        return result;
    }

    void ThrowInvalidIndexException(size_t requestedIndex, size_t vectorSize)
    {
        std::stringstream message;
        message << "Invalid index (" << requestedIndex << ") into the list of references, containing " << vectorSize << " elements.";
        throw std::invalid_argument(message.str());
    }

    std::string StandardCrossSectionSetup::ReferenceFileName(size_t index) const
    {
        if (index >= m_standardReferences.size())
        {
            ThrowInvalidIndexException(index, m_standardReferences.size());
        }

        return m_standardReferences[index].fileName;
    }

    std::string StandardCrossSectionSetup::ReferenceSpecieName(size_t index) const
    {
        if (index >= m_standardReferences.size())
        {
            ThrowInvalidIndexException(index, m_standardReferences.size());
        }

        return m_standardReferences[index].specieName;
    }

    bool StandardCrossSectionSetup::IsReferenceInVacuum(size_t index) const
    {
        if (index >= m_standardReferences.size())
        {
            ThrowInvalidIndexException(index, m_standardReferences.size());
        }

        return m_standardReferences[index].isVacuum;
    }

    std::string StandardCrossSectionSetup::FraunhoferReferenceFileName() const
    {
        return m_fraunhoferReference.fileName;
    }

}
