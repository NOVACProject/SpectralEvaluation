#include <SpectralEvaluation/Evaluation/FitWindow.h>
#include <SpectralEvaluation/StringUtils.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>

#include <iostream>

namespace novac
{

CFitWindow::CFitWindow(const CFitWindow& other)
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

CFitWindow& CFitWindow::operator=(const CFitWindow& other)
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

CFitWindow& CFitWindow::operator=(CFitWindow&& other)
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
    fitType = FIT_TYPE::FIT_HP_DIV;
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

void ReadReferences(CFitWindow& window)
{
    // For each reference in the fit-window, read it in and make sure that it exists...
    for (int referenceIndex = 0; referenceIndex < window.nRef; ++referenceIndex)
    {
        window.ref[referenceIndex].ReadCrossSectionDataFromFile();
    }

    if (window.fraunhoferRef.m_path.size() > 4)
    {
        window.fraunhoferRef.ReadCrossSectionDataFromFile();
    }

    // If children are defined, then read them as well
    for (CFitWindow& c : window.child)
    {
        ReadReferences(c);
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

void AddAsReference(CFitWindow& window, const std::vector<double>& referenceData, const std::string& name, int linkShiftToIdx)
{
    window.ref[window.nRef].m_data = std::make_unique<CCrossSectionData>(referenceData);
    window.ref[window.nRef].m_specieName = name;
    window.ref[window.nRef].m_columnOption = novac::SHIFT_TYPE::SHIFT_FREE;
    window.ref[window.nRef].m_columnValue = 1.0;
    window.ref[window.nRef].m_squeezeOption = novac::SHIFT_TYPE::SHIFT_FIX;
    window.ref[window.nRef].m_squeezeValue = 1.0;

    if (linkShiftToIdx >= 0)
    {
        window.ref[window.nRef].m_shiftOption = SHIFT_TYPE::SHIFT_LINK;
        window.ref[window.nRef].m_shiftValue = linkShiftToIdx;
    }
    else
    {
        window.ref[window.nRef].m_shiftOption = SHIFT_TYPE::SHIFT_FIX;
        window.ref[window.nRef].m_shiftValue = 0.0;
    }

    window.nRef += 1;
}

int AddAsSky(CFitWindow& window, const std::vector<double>& referenceData, SHIFT_TYPE shiftOption)
{
    int indexOfSkySpectrum = window.nRef;

    window.ref[window.nRef].m_data = std::make_unique<novac::CCrossSectionData>(referenceData);
    window.ref[window.nRef].m_specieName = "sky";
    window.ref[window.nRef].m_columnOption = novac::SHIFT_TYPE::SHIFT_FIX;
    window.ref[window.nRef].m_columnValue = -1.0;
    window.ref[window.nRef].m_squeezeOption = novac::SHIFT_TYPE::SHIFT_FIX;
    window.ref[window.nRef].m_squeezeValue = 1.0;
    window.ref[window.nRef].m_shiftOption = shiftOption;
    window.ref[window.nRef].m_shiftValue = 0.0;
    window.nRef += 1;

    return indexOfSkySpectrum;
}
}
