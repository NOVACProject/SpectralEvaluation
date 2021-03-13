#include <SpectralEvaluation/Evaluation/FitWindow.h>
#include <SpectralEvaluation/StringUtils.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>

#include <iostream>

namespace novac
{

CFitWindow::CFitWindow(const CFitWindow &other)
    : fitLow(other.fitLow),
    fitHigh(other.fitHigh),
    channel(other.channel),
    nRef(other.nRef),
    fraunhoferRef(other.fraunhoferRef),
    polyOrder(other.polyOrder),
    includeIntensitySpacePolyominal(other.includeIntensitySpacePolyominal),
    ringCalculation(other.ringCalculation),
    specLength(other.specLength),
    startChannel(other.startChannel),
    name(other.name),
    fitType(other.fitType),
    shiftSky(other.shiftSky),
    skyShift(other.skyShift),
    skySqueeze(other.skySqueeze),
    interlaceStep(other.interlaceStep),
    UV(other.UV),
    findOptimalShift(other.findOptimalShift),
    child(begin(other.child), end(other.child))
{
    for (int i = 0; i < other.nRef; ++i)
    {
        this->ref[i] = other.ref[i];
    }
}

CFitWindow::CFitWindow(CFitWindow&& other)
    : fitLow(other.fitLow),
    fitHigh(other.fitHigh),
    channel(other.channel),
    nRef(other.nRef),
    fraunhoferRef(other.fraunhoferRef),
    polyOrder(other.polyOrder),
    includeIntensitySpacePolyominal(other.includeIntensitySpacePolyominal),
    ringCalculation(other.ringCalculation),
    specLength(other.specLength),
    startChannel(other.startChannel),
    name(other.name),
    fitType(other.fitType),
    shiftSky(other.shiftSky),
    skyShift(other.skyShift),
    skySqueeze(other.skySqueeze),
    interlaceStep(other.interlaceStep),
    UV(other.UV),
    findOptimalShift(other.findOptimalShift),
    child(begin(other.child), end(other.child))
{
    for (int i = 0; i < other.nRef; ++i)
    {
        this->ref[i] = std::move(other.ref[i]);
    }
}

CFitWindow &CFitWindow::operator=(const CFitWindow &other)
{
    this->channel = other.channel;
    this->fitHigh = other.fitHigh;
    this->fitLow = other.fitLow;
    this->fitType = other.fitType;
    this->shiftSky = other.shiftSky;
    this->skyShift = other.skyShift;
    this->skySqueeze = other.skySqueeze;
    this->interlaceStep = other.interlaceStep;
    this->name = other.name;
    this->nRef = other.nRef;
    this->polyOrder = other.polyOrder;
    this->includeIntensitySpacePolyominal = other.includeIntensitySpacePolyominal;
    this->ringCalculation = other.ringCalculation;
    this->UV = other.UV;
    this->specLength = other.specLength;
    this->startChannel = other.startChannel;

    for (int i = 0; i < other.nRef; ++i)
    {
        this->ref[i] = other.ref[i];
    }
    this->fraunhoferRef = other.fraunhoferRef;
    this->findOptimalShift = other.findOptimalShift;
    this->child = std::vector<CFitWindow>(begin(other.child), end(other.child));
    return *this;
}

CFitWindow &CFitWindow::operator=(CFitWindow&& other)
{
    this->channel = other.channel;
    this->fitHigh = other.fitHigh;
    this->fitLow = other.fitLow;
    this->fitType = other.fitType;
    this->shiftSky = other.shiftSky;
    this->skyShift = other.skyShift;
    this->skySqueeze = other.skySqueeze;
    this->interlaceStep = other.interlaceStep;
    this->name = std::move(other.name);
    this->nRef = other.nRef;
    this->polyOrder = other.polyOrder;
    this->includeIntensitySpacePolyominal = other.includeIntensitySpacePolyominal;
    this->ringCalculation = other.ringCalculation;
    this->UV = other.UV;
    this->specLength = other.specLength;
    this->startChannel = other.startChannel;

    for (int i = 0; i < other.nRef; ++i)
    {
        this->ref[i] = std::move(other.ref[i]);
    }
    this->fraunhoferRef = other.fraunhoferRef;
    this->findOptimalShift = other.findOptimalShift;
    this->child = std::move(other.child);

    return *this;
}

void CFitWindow::Clear()
{
    fitHigh = 460;
    fitLow = 320;
    channel = 0;
    specLength = 2048;
    startChannel = 0;
    fitType = FIT_HP_DIV;
    shiftSky = true;
    interlaceStep = 1;
    name = "SO2";
    nRef = 0;
    polyOrder = 5;
    ringCalculation = RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING;
    includeIntensitySpacePolyominal = false;
    UV = true;
    for (int i = 0; i < MAX_N_REFERENCES; ++i)
    {
        ref[i].m_path = "";
        ref[i].m_specieName = "";
    }

    child.clear();

    fraunhoferRef.m_path = "";
    fraunhoferRef.m_specieName = "SolarSpec";
    findOptimalShift = false;
}

bool ReadReferences(CFitWindow& window)
{
    // For each reference in the fit-window, read it in and make sure that it exists...
    for (int referenceIndex = 0; referenceIndex < window.nRef; ++referenceIndex)
    {
        // Read in the cross section
        if (window.ref[referenceIndex].ReadCrossSectionDataFromFile())
        {
            std::cout << "Failed to read cross section file: " << window.ref[referenceIndex].m_path.c_str() << std::endl;
            return false;
        }
    }

    if (window.fraunhoferRef.m_path.size() > 4)
    {
        if (window.fraunhoferRef.ReadCrossSectionDataFromFile())
        {
            std::cout << "Failed to read Fraunhofer reference file: " << window.fraunhoferRef.m_path.c_str() << std::endl;
            return false;
        }
    }

    // If children are defined, then read them as well
    for (CFitWindow& c : window.child)
    {
        ReadReferences(c);
    }

    return true;
}

void HighPassFilterReferences(CFitWindow& window)
{
    // If children are defined, then handle them as well
    for (CFitWindow& c : window.child)
    {
        HighPassFilterReferences(c);
    }

    if (window.fitType != novac::FIT_HP_DIV && window.fitType != novac::FIT_HP_SUB)
    {
        return;
    }

    for (int referenceIndex = 0; referenceIndex < window.nRef; ++referenceIndex)
    {
        // Local handle for more convenient syntax.
        CReferenceFile& thisReference = window.ref[referenceIndex];

        if (!thisReference.m_isFiltered)
        {
            if (EqualsIgnoringCase(thisReference.m_specieName, "ring"))
            {
                HighPassFilter_Ring(*thisReference.m_data);
            }
            else
            {
                HighPassFilter(*thisReference.m_data);
            }
        }
    }
}

void ScaleReferencesToMolecCm2(CFitWindow& window)
{
    // If children are defined, then handle them as well
    for (CFitWindow& c : window.child)
    {
        ScaleReferencesToMolecCm2(c);
    }

    for (int referenceIndex = 0; referenceIndex < window.nRef; ++referenceIndex)
    {
        // Local handle for more convenient syntax.
        CReferenceFile& thisReference = window.ref[referenceIndex];

        if (thisReference.m_isFiltered)
        {
            // Convert from ppmm to moleculues / cm2
            Multiply(*thisReference.m_data, (1.0 / 2.5e15));
        }
    }
}

}