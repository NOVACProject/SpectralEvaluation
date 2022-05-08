#include <SpectralEvaluation/File/ScanEvaluationLogFileHandler.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>
#include <SpectralEvaluation/StringUtils.h>

#include <algorithm>
#include <numeric> // std::iota

using namespace novac;

CScanEvaluationLogFileHandler::CScanEvaluationLogFileHandler()
{
    // Defining which column contains which information
    for (int i = 0; i < MAX_N_REFERENCES_IN_EVALLOG; ++i) {
        m_tableColumnMapping.column[i] = 7;
        m_tableColumnMapping.columnError[i] = 8;
        m_tableColumnMapping.shift[i] = 9;
        m_tableColumnMapping.shiftError[i] = 10;
        m_tableColumnMapping.squeeze[i] = 11;
        m_tableColumnMapping.squeezeError[i] = 12;
    }
}

CScanEvaluationLogFileHandler::~CScanEvaluationLogFileHandler(void)
{
    m_scan.clear();
}

void CScanEvaluationLogFileHandler::ParseScanHeader(const char szLine[8192]) {
    // reset some old information
    ResetColumns();

    char str[8192];
    if (szLine[0] == '#')
        strncpy(str, szLine + 1, 8191 * sizeof(char));
    else
        strncpy(str, szLine, 8192 * sizeof(char));

    char* szToken = str;
    int curCol = -1;
    char elevation[] = "elevation";
    char scanAngle[] = "scanangle";
    char obsAngle[] = "observationangle";
    char azimuth[] = "azimuth";
    char column[] = "column";
    char columnError[] = "columnerror";
    char intensity[] = "intensity";	// peak intensity
    char fitIntensity[] = "intens(fitregion)"; // fit-region intensity
    char fitIntensity2[] = "fitintensity"; // fit-region intensity
    char peakSat[] = "specsaturation";	// maximum saturation ratio of the whole spectrum
    char fitSat[] = "fitsaturation";	// maximum saturation ratio in the fit region
    char delta[] = "delta";
    char chiSquare[] = "chisquare";
    char shift[] = "shift";
    char shiftError[] = "shifterror";
    char squeeze[] = "squeeze";
    char squeezeError[] = "squeezeerror";
    char exposureTime[] = "exposuretime";
    char numSpec[] = "numSpec";
    char offset[] = "offset";
    char starttime[] = "starttime";
    char stoptime[] = "stoptime";
    char nameStr[] = "name";

    while (true)
    {
        szToken = strtok(szToken, "\t");
        if (szToken == nullptr)
        {
            break;
        }

        ++curCol;

        // The scan-angle (previously known as elevation)
        if (0 == _strnicmp(szToken, elevation, strlen(elevation))) {
            m_tableColumnMapping.position = curCol;
            szToken = nullptr;
            continue;
        }

        // The scan-angle (previously known as elevation)
        if (0 == _strnicmp(szToken, scanAngle, strlen(scanAngle))) {
            m_tableColumnMapping.position = curCol;
            szToken = nullptr;
            continue;
        }

        // The observation-angle (the scan-angle for the heidelberg instrument)
        if (0 == _strnicmp(szToken, obsAngle, strlen(obsAngle))) {
            m_tableColumnMapping.position = curCol;
            szToken = nullptr;
            continue;
        }

        // The azimuth-angle (defined for the heidelberg instrument)
        if (0 == _strnicmp(szToken, azimuth, strlen(azimuth))) {
            m_tableColumnMapping.position2 = curCol;
            szToken = nullptr;
            continue;
        }

        // The exposure time
        if (0 == _strnicmp(szToken, exposureTime, strlen(exposureTime))) {
            m_tableColumnMapping.expTime = curCol;
            szToken = nullptr;
            continue;
        }

        // The start time
        if (0 == _strnicmp(szToken, starttime, strlen(starttime))) {
            m_tableColumnMapping.starttime = curCol;
            szToken = nullptr;
            continue;
        }

        // The stop time
        if (0 == _strnicmp(szToken, stoptime, strlen(stoptime))) {
            m_tableColumnMapping.stoptime = curCol;
            szToken = nullptr;
            continue;
        }

        // The name of the spectrum
        if (0 == _strnicmp(szToken, nameStr, strlen(nameStr))) {
            m_tableColumnMapping.name = curCol;
            szToken = nullptr;
            continue;
        }

        // The number of co-added spectra
        if (0 == _strnicmp(szToken, numSpec, strlen(numSpec))) {
            m_tableColumnMapping.nSpec = curCol;
            szToken = nullptr;
            continue;
        }

        // The offset
        if (0 == _strnicmp(szToken, offset, strlen(offset))) {
            m_tableColumnMapping.offset = curCol;
            szToken = nullptr;
            continue;
        }

        // The column error (must be looked for before 'column')
        if (0 == _strnicmp(szToken, columnError, strlen(columnError))) {
            m_tableColumnMapping.columnError[m_evResult.NumberOfSpecies() - 1] = curCol;
            szToken = nullptr;
            continue;
        }

        // The column
        if (0 == _strnicmp(szToken, column, strlen(column))) {
            m_tableColumnMapping.column[m_evResult.NumberOfSpecies()] = curCol;
            char* pt = szToken + strlen(column) + 1;
            szToken[strlen(szToken) - 1] = 0;
            std::string specieName(pt);
            m_evResult.InsertSpecie(specieName);
            ++m_tableColumnMapping.nSpecies;
            szToken = nullptr;
            continue;
        }

        // The shift error (must be checked before 'shift')
        if (0 == _strnicmp(szToken, shiftError, strlen(shiftError))) {
            m_tableColumnMapping.shiftError[m_evResult.NumberOfSpecies() - 1] = curCol;
            szToken = nullptr;
            continue;
        }

        // The shift
        if (0 == _strnicmp(szToken, shift, strlen(shift))) {
            m_tableColumnMapping.shift[m_evResult.NumberOfSpecies() - 1] = curCol;
            szToken = nullptr;
            continue;
        }

        // The squeeze error (must be checked before 'squeeze')
        if (0 == _strnicmp(szToken, squeezeError, strlen(squeezeError))) {
            m_tableColumnMapping.squeezeError[m_evResult.NumberOfSpecies() - 1] = curCol;
            szToken = nullptr;
            continue;
        }

        // The squeeze
        if (0 == _strnicmp(szToken, squeeze, strlen(squeeze))) {
            m_tableColumnMapping.squeeze[m_evResult.NumberOfSpecies() - 1] = curCol;
            szToken = nullptr;
            continue;
        }

        // The spectrum peak-intensity
        if (0 == _strnicmp(szToken, intensity, strlen(intensity))) {
            m_tableColumnMapping.intensity = curCol;
            szToken = nullptr;
            continue;
        }

        // The spectrum fit-intensity
        if (0 == _strnicmp(szToken, fitIntensity, strlen(fitIntensity)) ||
            0 == _strnicmp(szToken, fitIntensity2, strlen(fitIntensity2))) {

            m_tableColumnMapping.fitIntensity = curCol;
            szToken = nullptr;
            continue;
        }

        // The spectrum maximum saturation ratio of the whole spectrum
        if (0 == _strnicmp(szToken, peakSat, strlen(peakSat))) {
            m_tableColumnMapping.peakSaturation = curCol;
            szToken = nullptr;
            continue;
        }

        // The spectrum maximum saturation ratio in the fit region
        if (0 == _strnicmp(szToken, fitSat, strlen(fitSat))) {
            m_tableColumnMapping.fitSaturation = curCol;
            szToken = nullptr;
            continue;
        }

        // The delta of the fit
        if (0 == _strnicmp(szToken, delta, strlen(delta))) {
            m_tableColumnMapping.delta = curCol;
            szToken = nullptr;
            continue;
        }

        // The chi-square of the fit
        if (0 == _strnicmp(szToken, chiSquare, strlen(chiSquare))) {
            m_tableColumnMapping.chiSquare = curCol;
            szToken = nullptr;
            continue;
        }

        szToken = nullptr;
    }

    m_specie.resize(m_evResult.NumberOfSpecies());
    for (size_t k = 0; k < m_evResult.NumberOfSpecies(); ++k)
    {
        m_specie[k] = m_evResult.m_referenceResult[k].m_specieName;
    }

    return;
}

bool CScanEvaluationLogFileHandler::ReadEvaluationLog(const std::string& evaluationLogFile)
{
    char  expTimeStr[] = "exposuretime";         // this string only exists in the header line.
    char  scanInformation[] = "<scaninformation>";    // this string only exists in the scan-information section before the scan-data
    char  fluxInformation[] = "<fluxinfo>";           // this string only exists in the flux-information section before the scan-data
    char  spectralData[] = "<spectraldata>";
    char  endofSpectralData[] = "</spectraldata>";
    std::string str;
    char szLine[8192];
    int measNr = 0;
    double fValue;
    bool fReadingScan = false;
    double flux = 0.0;

    // If no evaluation log selected, quit
    if (evaluationLogFile.size() <= 1)
    {
        return false;
    }

    m_evaluationLog = evaluationLogFile;

    // First count the number of scans in the file.
    // This to speed up the initialization of the arrays
    std::vector<CDateTime> allStartTimes;
    GetScanStartTimes(allStartTimes);
    if (allStartTimes.size() == 0)
    {
        return false;
    }
    m_scan.resize(allStartTimes.size());
    // m_windField.SetSize(nScans + 1);

    std::vector<unsigned int> sortOrder;
    SortScanStartTimes(allStartTimes, sortOrder);

    int m_scanNum = 0; // TODO: Rename

    // Open the evaluation log. (Notice that the NovacProgram did use a CriticalSection here for locking..)
    {
        FILE* f = fopen(m_evaluationLog.c_str(), "r");
        if (nullptr == f) {
            return false;
        }

        // Reset the column- and spectrum info
        ResetColumns();
        ResetScanInformation();

        // Read the file, one line at a time
        while (fgets(szLine, 8192, f)) {

            // ignore empty lines
            if (strlen(szLine) < 2) {
                if (fReadingScan) {
                    fReadingScan = false;
                    // Reset the column- and spectrum-information
                    ResetColumns();
                    ResetScanInformation();
                }
                continue;
            }

            // convert the string to all lower-case letters
            for (unsigned int it = 0; it < strlen(szLine); ++it) {
                szLine[it] = (char)tolower(szLine[it]);
            }

            // find the next scan-information section
            if (nullptr != strstr(szLine, scanInformation)) {
                ParseScanInformation(m_specInfo, flux, f);
                continue;
            }

            // find the next flux-information section
            if (nullptr != strstr(szLine, fluxInformation)) {
                // ParseFluxInformation(m_windField[sortOrder[m_scanNum + 1]], flux, f);
                continue;
            }

            if (nullptr != strstr(szLine, spectralData)) {
                fReadingScan = true;
                continue;
            }
            else if (nullptr != strstr(szLine, endofSpectralData)) {
                fReadingScan = false;
                continue;
            }

            // find the next start of a scan 
            if (nullptr != strstr(szLine, expTimeStr))
            {
                // check so that there was some information in the last scan read
                //	if not the re-use the memory space
                if ((measNr > 0) || (measNr == 0 && m_scanNum < 0)) {

                    // The current measurement position inside the scan
                    measNr = 0;

                    // before we start the next scan, calculate some information about
                    // the old one

                    // 1. If the sky and dark were specified, remove them from the measurement
                    if (m_scanNum >= 0 && fabs(m_scan[sortOrder[m_scanNum]].m_specInfo[1].m_scanAngle - 180.0) < 1) {
                        m_scan[sortOrder[m_scanNum]].RemoveResult(0); // remove sky
                        m_scan[sortOrder[m_scanNum]].RemoveResult(0); // remove dark
                    }

                    // 2. Calculate the offset
                    if (m_scanNum >= 0 && m_scan[sortOrder[m_scanNum]].m_spec.size() > 0) {
                        // const CEvaluationResult& evResult = m_scan[sortOrder[m_scanNum]].m_spec.front();
                        // if (evResult.m_referenceResult.size() > 0) {
                        //     m_scan[sortOrder[m_scanNum]].CalculateOffset(evResult.m_referenceResult.front().m_specieName);
                        // }
                    }

                    // start the next scan.
                    ++m_scanNum;
                }

                // This line is the header line which says what each column represents.
                //  Read it and parse it to find out how to interpret the rest of the 
                //  file. 
                ParseScanHeader(szLine);

                // start parsing the lines
                fReadingScan = true;

                // read the next line, which is the first line in the scan
                continue;
            }

            // ignore comment lines
            if (szLine[0] == '#')
                continue;

            // if we're not reading a scan, let's read the next line
            if (!fReadingScan)
                continue;

            // Split the scan information up into tokens and parse them. 
            char* szToken = (char*)szLine;
            int curCol = -1;
            while (true)
            {
                szToken = strtok(szToken, " \t");
                if (szToken == nullptr)
                {
                    break;
                }

                ++curCol;

                // First check the starttime
                if (curCol == m_tableColumnMapping.starttime) {
                    int fValue1, fValue2, fValue3, ret;
                    if (strstr(szToken, ":")) {
                        ret = sscanf(szToken, "%d:%d:%d", &fValue1, &fValue2, &fValue3);
                    }
                    else {
                        ret = sscanf(szToken, "%d.%d.%d", &fValue1, &fValue2, &fValue3);
                    }
                    if (ret == 3) {
                        m_specInfo.m_startTime.hour = (unsigned char)fValue1;
                        m_specInfo.m_startTime.minute = (unsigned char)fValue2;
                        m_specInfo.m_startTime.second = (unsigned char)fValue3;
                        szToken = nullptr;
                    }
                    continue;
                }

                // Then check the stoptime
                if (curCol == m_tableColumnMapping.stoptime) {
                    int fValue1, fValue2, fValue3, ret;
                    if (strstr(szToken, ":")) {
                        ret = sscanf(szToken, "%d:%d:%d", &fValue1, &fValue2, &fValue3);
                    }
                    else {
                        ret = sscanf(szToken, "%d.%d.%d", &fValue1, &fValue2, &fValue3);
                    }
                    if (ret == 3) {
                        m_specInfo.m_stopTime.hour = (unsigned char)fValue1;
                        m_specInfo.m_stopTime.minute = (unsigned char)fValue2;
                        m_specInfo.m_stopTime.second = (unsigned char)fValue3;
                        szToken = nullptr;
                    }
                    continue;
                }

                // Also check the name...
                if (curCol == m_tableColumnMapping.name) {
                    m_specInfo.m_name = std::string(szToken);
                    szToken = nullptr;
                    continue;
                }

                // ignore columns whose value cannot be parsed into a float
                if (1 != sscanf(szToken, "%lf", &fValue)) {
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_tableColumnMapping.position) {
                    m_specInfo.m_scanAngle = (float)fValue;
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_tableColumnMapping.position2) {
                    m_specInfo.m_scanAngle2 = (float)fValue;
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_tableColumnMapping.intensity) {
                    m_specInfo.m_peakIntensity = (float)fValue;
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_tableColumnMapping.fitIntensity) {
                    m_specInfo.m_fitIntensity = (float)fValue;
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_tableColumnMapping.fitSaturation) {
                    m_specInfo.m_fitIntensity = (float)fValue;
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_tableColumnMapping.peakSaturation) {
                    m_specInfo.m_peakIntensity = (float)fValue;
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_tableColumnMapping.offset) {
                    m_specInfo.m_offset = (float)fValue;
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_tableColumnMapping.delta) {
                    m_evResult.m_delta = (float)fValue;
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_tableColumnMapping.chiSquare) {
                    m_evResult.m_chiSquare = (float)fValue;
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_tableColumnMapping.nSpec) {
                    m_specInfo.m_numSpec = (long)fValue;
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_tableColumnMapping.expTime) {
                    m_specInfo.m_exposureTime = (long)fValue;
                    szToken = nullptr;
                    continue;
                }

                for (int k = 0; k < m_tableColumnMapping.nSpecies; ++k) {
                    if (curCol == m_tableColumnMapping.column[k]) {
                        m_evResult.m_referenceResult[k].m_column = (float)fValue;
                        break;
                    }
                    if (curCol == m_tableColumnMapping.columnError[k]) {
                        m_evResult.m_referenceResult[k].m_columnError = (float)fValue;
                        break;
                    }
                    if (curCol == m_tableColumnMapping.shift[k]) {
                        m_evResult.m_referenceResult[k].m_shift = (float)fValue;
                        break;
                    }
                    if (curCol == m_tableColumnMapping.shiftError[k]) {
                        m_evResult.m_referenceResult[k].m_shiftError = (float)fValue;
                        break;
                    }
                    if (curCol == m_tableColumnMapping.squeeze[k]) {
                        m_evResult.m_referenceResult[k].m_squeeze = (float)fValue;
                        break;
                    }
                    if (curCol == m_tableColumnMapping.squeezeError[k]) {
                        m_evResult.m_referenceResult[k].m_squeezeError = (float)fValue;
                        break;
                    }
                }
                szToken = nullptr;
            }

            // start reading the next line in the evaluation log (i.e. the next
            //  spectrum in the scan). Insert the data from this spectrum into the 
            //  CScanResult structure

            // If this is the first spectrum in the new scan, then make
            //	an initial guess for how large the arrays are going to be...
            if (measNr == 0 && m_scanNum > 1) {
                // If this is the first spectrum in a new scan, then initialize the 
                //	size of the arrays, to save some time on re-allocating memory
                m_scan[sortOrder[m_scanNum]].InitializeArrays((long)m_scan[m_scanNum - 1].m_spec.size());
            }

            m_specInfo.m_scanIndex = (short)measNr;
            if (EqualsIgnoringCase(m_specInfo.m_name, "sky")) {
                m_scan[sortOrder[m_scanNum]].m_skySpecInfo = m_specInfo;
            }
            else if (EqualsIgnoringCase(m_specInfo.m_name, "dark")) {
                m_scan[sortOrder[m_scanNum]].m_darkSpecInfo = m_specInfo;
            }
            else if (EqualsIgnoringCase(m_specInfo.m_name, "offset")) {
                m_scan[sortOrder[m_scanNum]].m_offsetSpecInfo = m_specInfo;
            }
            else if (EqualsIgnoringCase(m_specInfo.m_name, "dark_cur")) {
                m_scan[sortOrder[m_scanNum]].m_darkCurSpecInfo = m_specInfo;
            }
            else {
                m_scan[sortOrder[m_scanNum]].AppendResult(m_evResult, m_specInfo);
                // m_scan[sortOrder[m_scanNum]].SetFlux(flux);
            }

            double dynamicRange = 1.0; // <-- unknown
            if (m_tableColumnMapping.peakSaturation != -1) { // If the intensity is specified as a saturation ratio...
                dynamicRange = CSpectrometerDatabase::GetInstance().GetModel(m_specInfo.m_specModelName).maximumIntensityForSingleReadout;
            }
            // m_scan[sortOrder[m_scanNum]].CheckGoodnessOfFit(m_specInfo);
            ++measNr;
        }

        // close the evaluation log
        fclose(f);
    }

    // If the sky and dark were specified, remove them from the measurement
    if (m_scanNum >= 0 && fabs(m_scan[sortOrder[m_scanNum]].m_specInfo[1].m_scanAngle - 180.0) < 1) {
        m_scan[sortOrder[m_scanNum]].RemoveResult(0); // remove sky
        m_scan[sortOrder[m_scanNum]].RemoveResult(0); // remove dark
    }

    // Calculate the offset
    if (m_scanNum >= 0) {
        // m_scan[sortOrder[m_scanNum]].CalculateOffset(m_evResult.m_referenceResult[0].m_specieName);
    }

    // make sure that scan num is correct
    ++m_scanNum;

    // Sort the scans in order of collection
    SortScans();

    return true;
}

/** Makes a quick scan through the evaluation-log
    to count the number of scans in it */
long CScanEvaluationLogFileHandler::CountScansInFile() {
    char  expTimeStr[] = "exposuretime"; // this string only exists in the header line.
    char szLine[8192];
    long  nScans = 0;

    // If no evaluation log selected, quit
    if (m_evaluationLog.size() <= 1)
    {
        return 0;
    }

    // Open the evaluation log. (Notice that the NovacProgram did use a CriticalSection here for locking..)
    {

        // Open the evaluation log
        FILE* f = fopen(m_evaluationLog.c_str(), "r");
        if (nullptr == f) {
            return 0;
        }

        // Read the file, one line at a time
        while (fgets(szLine, 8192, f)) {
            // convert the string to all lower-case letters
            for (unsigned int it = 0; it < strlen(szLine); ++it) {
                szLine[it] = (char)tolower(szLine[it]);
            }

            // find the next start of a scan 
            if (nullptr != strstr(szLine, expTimeStr)) {
                ++nScans;
            }
        }

        fclose(f);

    }

    // Return the number of scans found in the file
    return nScans;
}

void CScanEvaluationLogFileHandler::GetScanStartTimes(std::vector<CDateTime>& startTimes)
{
    char  scanInfoStart[] = "<scaninformation>"; // this string indicates the beginning of a 'scanInformation' section
    char  scanInfoStop[] = "</scaninformation>"; // this string indicates the beginning of a 'scanInformation' section
    char szLine[8192];
    long  nScans = 0;
    bool inScanInfoSection = false;
    CDateTime scanStartTime;
    int tmpInt[3];
    bool foundDate = false, foundTime = false;

    startTimes.clear();

    // If no evaluation log selected, quit
    if (m_evaluationLog.size() <= 1)
    {
        return;
    }

    {
        // Open the evaluation log
        FILE* f = fopen(m_evaluationLog.c_str(), "r");
        if (nullptr == f) {
            return;
        }

        // Read the file, one line at a time
        while (fgets(szLine, 8192, f)) {
            // convert the string to all lower-case letters
            for (unsigned int it = 0; it < strlen(szLine); ++it) {
                szLine[it] = (char)tolower(szLine[it]);
            }

            // find the next start of a scan 
            if (nullptr != strstr(szLine, scanInfoStart)) {
                inScanInfoSection = true;
            }

            if (nullptr != strstr(szLine, scanInfoStop)) {
                if (foundDate && foundTime) {
                    startTimes.push_back(CDateTime(scanStartTime));
                }
                inScanInfoSection = false;
                foundDate = false;
                foundTime = false;
                ++nScans;
            }

            if (inScanInfoSection) {
                const char* pt = strstr(szLine, "date=");
                if (pt != nullptr)
                {
                    if (3 == sscanf(pt + 5, "%d.%d.%d", &tmpInt[0], &tmpInt[1], &tmpInt[2])) {
                        scanStartTime.day = (unsigned char)tmpInt[0];
                        scanStartTime.month = (unsigned char)tmpInt[1];
                        scanStartTime.year = (unsigned short)tmpInt[2];
                        foundDate = true;
                        continue;
                    }
                }
                pt = strstr(szLine, "starttime=");
                if (pt != nullptr)
                {
                    if (3 == sscanf(pt + 10, "%d:%d:%d", &tmpInt[0], &tmpInt[1], &tmpInt[2])) {
                        scanStartTime.hour = (unsigned char)tmpInt[0];
                        scanStartTime.minute = (unsigned char)tmpInt[1];
                        scanStartTime.second = (unsigned char)tmpInt[2];
                        foundTime = true;
                    }
                    continue;
                }
            }
        }
        fclose(f);
    }
}

/** Reads and parses the 'scanInfo' header before the scan */
void CScanEvaluationLogFileHandler::ParseScanInformation(CSpectrumInfo& scanInfo, double& flux, FILE* f) {
    char szLine[8192];
    int tmpInt[3];
    double tmpDouble;

    // Reset the column- and spectrum info
    // ResetColumns();
    ResetScanInformation();

    // read the additional scan-information, line by line
    while (fgets(szLine, 8192, f)) {

        // convert to lower-case
        for (unsigned int it = 0; it < strlen(szLine); ++it) {
            szLine[it] = (char)tolower(szLine[it]);
        }

        const char* pt = strstr(szLine, "</scaninformation>");
        if (pt != nullptr)
        {
            break;
        }

        pt = strstr(szLine, "compiledate=");
        if (pt != nullptr)
        {
            continue;
        }

        pt = strstr(szLine, "site=");
        if (pt != nullptr)
        {
            scanInfo.m_site = std::string(pt + 5);
            Remove(scanInfo.m_site, '\n'); // Remove newline characters
            continue;
        }

        pt = strstr(szLine, "date=");
        if (pt != nullptr)
        {
            if (3 == sscanf(pt + 5, "%d.%d.%d", &tmpInt[0], &tmpInt[1], &tmpInt[2])) {
                scanInfo.m_startTime.year = (unsigned short)tmpInt[2];
                scanInfo.m_startTime.month = (unsigned char)tmpInt[1];
                scanInfo.m_startTime.day = (unsigned char)tmpInt[0];

                scanInfo.m_stopTime.year = scanInfo.m_startTime.year;
                scanInfo.m_stopTime.month = scanInfo.m_startTime.month;
                scanInfo.m_stopTime.day = scanInfo.m_startTime.day;
            }
            continue;
        }

        pt = strstr(szLine, "starttime=");
        if (pt != nullptr)
        {
            if (3 == sscanf(pt + 10, "%d:%d:%d", &tmpInt[0], &tmpInt[1], &tmpInt[2])) {
                scanInfo.m_startTime.hour = (unsigned char)tmpInt[0];
                scanInfo.m_startTime.minute = (unsigned char)tmpInt[1];
                scanInfo.m_startTime.second = (unsigned char)tmpInt[2];
            }
            continue;
        }

        pt = strstr(szLine, "stoptime=");
        if (pt != nullptr)
        {
            if (3 == sscanf(pt + 9, "%d.%d.%d", &tmpInt[0], &tmpInt[1], &tmpInt[2])) {
                scanInfo.m_stopTime.hour = (unsigned char)tmpInt[0];
                scanInfo.m_stopTime.minute = (unsigned char)tmpInt[1];
                scanInfo.m_stopTime.second = (unsigned char)tmpInt[2];
            }
            continue;
        }

        pt = strstr(szLine, "compass=");
        if (pt != nullptr)
        {
            if (sscanf(pt + 8, "%lf", &tmpDouble) == 1) {
                scanInfo.m_compass = (float)fmod(tmpDouble, 360.0);
            }
            continue;
        }

        pt = strstr(szLine, "tilt=");
        if (pt != nullptr)
        {
            if (sscanf(pt + 5, "%lf", &tmpDouble) == 1) {
                scanInfo.m_pitch = (float)tmpDouble;
            }
            continue;
        }

        pt = strstr(szLine, "lat=");
        if (pt != nullptr)
        {
            if (sscanf(pt + 4, "%lf", &tmpDouble) == 1) {
                scanInfo.m_gps.m_latitude = tmpDouble;
            }
            continue;
        }

        pt = strstr(szLine, "long=");
        if (pt != nullptr)
        {
            if (sscanf(pt + 5, "%lf", &tmpDouble) == 1) {
                scanInfo.m_gps.m_longitude = tmpDouble;
            }
            continue;
        }

        pt = strstr(szLine, "alt=");
        if (pt != nullptr)
        {
            if (sscanf(pt + 4, "%lf", &tmpDouble) == 1) {
                scanInfo.m_gps.m_altitude = (long)tmpDouble;
            }
            continue;
        }

        pt = strstr(szLine, "serial=");
        if (pt != nullptr)
        {
            scanInfo.m_device = std::string(pt + 7);
            Remove(scanInfo.m_device, '\n'); // remove remaining strange things in the serial-number
            MakeUpper(scanInfo.m_device);	// Convert the serial-number to all upper case letters

            // Extract the spectrometer-model from the serial-number of the spectrometer
            SpectrometerModel spectrometer = CSpectrometerDatabase::GetInstance().GuessModelFromSerial(scanInfo.m_device);
            scanInfo.m_specModelName = spectrometer.modelName;
            scanInfo.m_average = spectrometer.averagesSpectra;
            continue;
        }

        pt = strstr(szLine, "spectrometer=");
        if (pt != nullptr)
        {
            // TODO:: read in the spectrometer model somewhere
            continue;
        }

        pt = strstr(szLine, "volcano=");
        if (pt != nullptr)
        {
            scanInfo.m_volcano = std::string(pt + 8);
            Remove(scanInfo.m_volcano, '\n'); // Remove newline characters
            continue;
        }

        pt = strstr(szLine, "observatory=");
        if (pt != nullptr)
        {
            scanInfo.m_observatory = std::string(pt + 12);
            Remove(scanInfo.m_observatory, '\n'); // Remove newline characters
            continue;
        }

        pt = strstr(szLine, "channel=");
        if (pt != nullptr)
        {
            if (sscanf(pt + 8, "%lf", &tmpDouble) == 1) {
                scanInfo.m_channel = (unsigned char)tmpDouble;
            }
        }

        pt = strstr(szLine, "coneangle=");
        if (pt != nullptr)
        {
            if (sscanf(pt + 10, "%lf", &tmpDouble) == 1) {
                scanInfo.m_coneAngle = (float)tmpDouble;
            }
        }

        pt = strstr(szLine, "flux=");
        if (pt != nullptr)
        {
            if (sscanf(pt + 5, "%lf", &tmpDouble) == 1) {
                flux = tmpDouble;
            }
        }

        pt = strstr(szLine, "battery=");
        if (pt != nullptr)
        {
            (void)sscanf(pt + 8, "%f", &scanInfo.m_batteryVoltage);
        }

        pt = strstr(szLine, "temperature");
        if (pt != nullptr)
        {
            (void)sscanf(pt + 12, "%f", &scanInfo.m_temperature);
        }
    }
}

/* void CScanEvaluationLogFileHandler::ParseFluxInformation(CWindField& windField, double& flux, FILE* f) {
    char szLine[8192];
    char* pt = nullptr;
    double windSpeed = 10, windDirection = 0, plumeHeight = 1000;
    MET_SOURCE windSpeedSource = MET_USER;
    MET_SOURCE windDirectionSource = MET_USER;
    MET_SOURCE plumeHeightSource = MET_USER;
    char source[512];

    // read the additional scan-information, line by line
    while (fgets(szLine, 8192, f)) {
        if (pt = strstr(szLine, "</fluxinfo>")) {
            // save all the values
            windField.SetPlumeHeight(plumeHeight, plumeHeightSource);
            windField.SetWindDirection(windDirection, windDirectionSource);
            windField.SetWindSpeed(windSpeed, windSpeedSource);
            break;
        }

        if (pt = strstr(szLine, "flux=")) {
            int ret = sscanf(pt + 5, "%lf", &flux);
            continue;
        }

        if (pt = strstr(szLine, "windspeed=")) {
            int ret = sscanf(pt + 10, "%lf", &windSpeed);
            continue;
        }

        if (pt = strstr(szLine, "winddirection=")) {
            int ret = sscanf(pt + 14, "%lf", &windDirection);
            continue;
        }

        if (pt = strstr(szLine, "plumeheight=")) {
            int ret = sscanf(pt + 12, "%lf", &plumeHeight);
            continue;
        }

        if (pt = strstr(szLine, "windspeedsource=")) {
            if (sscanf(pt + 16, "%s", source) == 1) {
                source[0] = '\0';
                if (strstr(source, "user")) {
                    windSpeedSource = MET_USER;
                }
                else if (strstr(source, "default")) {
                    windSpeedSource = MET_DEFAULT;
                }
                else if (strstr(source, "ecmwf_forecast")) {
                    windSpeedSource = MET_ECMWF_FORECAST;
                }
                else if (strstr(source, "ecmwf_analysis")) {
                    windSpeedSource = MET_ECMWF_ANALYSIS;
                }
                else if (strstr(source, "dual_beam_measurement")) {
                    windSpeedSource = MET_DUAL_BEAM_MEASUREMENT;
                }
            }
            continue;
        }

        if (pt = strstr(szLine, "winddirectionsource=")) {
            if (sscanf(pt + 20, "%s", source) == 1) {
                source[0] = '\0';
                if (strstr(source, "user")) {
                    windDirectionSource = MET_USER;
                }
                else if (strstr(source, "default")) {
                    windDirectionSource = MET_DEFAULT;
                }
                else if (strstr(source, "ecmwf_forecast")) {
                    windDirectionSource = MET_ECMWF_FORECAST;
                }
                else if (strstr(source, "ecmwf_analysis")) {
                    windDirectionSource = MET_ECMWF_ANALYSIS;
                }
                else if (strstr(source, "triangulation")) {
                    windDirectionSource = MET_GEOMETRY_CALCULATION;
                }
            }
            continue;
        }

        if (pt = strstr(szLine, "plumeheightsource=")) {
            if (sscanf(pt + 18, "%s", source) == 1) {
                if (strstr(source, "user")) {
                    plumeHeightSource = MET_USER;
                }
                else if (strstr(source, "default")) {
                    plumeHeightSource = MET_DEFAULT;
                }
                else if (strstr(source, "ecmwf_forecast")) {
                    plumeHeightSource = MET_ECMWF_FORECAST;
                }
                else if (strstr(source, "ecmwf_analysis")) {
                    plumeHeightSource = MET_ECMWF_ANALYSIS;
                }
                else if (strstr(source, "triangulation")) {
                    plumeHeightSource = MET_GEOMETRY_CALCULATION;
                }
            }
            continue;
        }
    }
}
*/

void CScanEvaluationLogFileHandler::ResetColumns()
{
    for (int k = 0; k < MAX_N_REFERENCES_IN_EVALLOG; ++k) {
        m_tableColumnMapping.column[k] = -1;
        m_tableColumnMapping.columnError[k] = -1;
        m_tableColumnMapping.shift[k] = -1;
        m_tableColumnMapping.shiftError[k] = -1;
        m_tableColumnMapping.squeeze[k] = -1;
        m_tableColumnMapping.squeezeError[k] = -1;
    }
    m_tableColumnMapping.delta = m_tableColumnMapping.intensity = m_tableColumnMapping.position = m_tableColumnMapping.position2 = -1;
    m_tableColumnMapping.nSpecies = 0;
    m_evResult.m_referenceResult.clear();
    m_tableColumnMapping.expTime = m_tableColumnMapping.nSpec = -1;
    m_tableColumnMapping.name = -1;
}

void CScanEvaluationLogFileHandler::ResetScanInformation()
{
    m_specInfo.m_channel = 0;
    m_specInfo.m_compass = m_specInfo.m_scanAngle = 0.0;
    m_tableColumnMapping.starttime = -1; m_tableColumnMapping.stoptime = -1;
}

bool CScanEvaluationLogFileHandler::IsWindSpeedMeasurement(int /*scanNo*/)
{
    // TODO: Implement
    return false;

    /* // check so that there are some scans read, and that the scan index is ok
    if (m_scanNum < 1 || scanNo > m_scanNum || scanNo < 0)
        return false;

    return m_scan[scanNo].IsWindMeasurement(); */

}

void CScanEvaluationLogFileHandler::SortScans()
{
    std::sort(begin(m_scan), end(m_scan), [](const BasicScanEvaluationResult& result1, const BasicScanEvaluationResult& result2) {
        return result1.m_skySpecInfo.m_startTime < result2.m_skySpecInfo.m_startTime;
    });
}

void CScanEvaluationLogFileHandler::SortScanStartTimes(std::vector<CDateTime>& timestamps, std::vector<unsigned int>& sortOrder)
{
    sortOrder.resize(timestamps.size());
    std::iota(begin(sortOrder), end(sortOrder), 0);

    std::sort(begin(sortOrder), end(sortOrder), [&](unsigned int idx1, unsigned int idx2) {
        return timestamps[idx1] < timestamps[idx2];
    });
}

bool CScanEvaluationLogFileHandler::WriteEvaluationLog(const std::string fileName) {
    std::string string, specieName;
    std::string wsSrc, wdSrc, phSrc;
    CDateTime startTime;
    /*
    // 1. Test if the file already exists, if so then return false
    if (IsExistingFile(fileName))
        return false;

    // 2. Write the file
    FILE* f = fopen(fileName, "w");

    for (int scanIndex = 0; scanIndex < this->m_scanNum; ++scanIndex) {
        Evaluation::CScanResult& scan = this->m_scan[scanIndex];
        CWindField& wind = this->m_windField[scanIndex];

        scan.GetStartTime(0, startTime);

        // ----------------- Create the additional scan-information ----------------------
        string.Format("\n<scaninformation>\n");
        string.AppendFormat("\tdate=%02d.%02d.%04d\n", startTime.day, startTime.month, startTime.year);
        string.AppendFormat("\tstarttime=%02d:%02d:%02d\n", startTime.hour, startTime.minute, startTime.second);
        string.AppendFormat("\tcompass=%.1lf\n", scan.GetCompass());
        string.AppendFormat("\ttilt=%.1lf\n", scan.GetPitch());
        string.AppendFormat("\tlat=%.6lf\n", scan.GetLatitude());
        string.AppendFormat("\tlong=%.6lf\n", scan.GetLongitude());
        string.AppendFormat("\talt=%ld\n", (int)scan.GetAltitude());

        string.AppendFormat("\tvolcano=%s\n", m_specInfo.m_volcano.c_str());
        string.AppendFormat("\tsite=%s\n", m_specInfo.m_site.c_str());
        string.AppendFormat("\tobservatory=%s\n", m_specInfo.m_observatory.c_str());

        string.AppendFormat("\tserial=%s\n", scan.GetSerial().c_str());
        string.AppendFormat("\tspectrometer=%s\n", m_specInfo.m_specModelName.c_str());
        string.AppendFormat("\tchannel=%d\n", m_specInfo.m_channel);
        string.AppendFormat("\tconeangle=%.1lf\n", scan.GetConeAngle());
        string.AppendFormat("\tinterlacesteps=%d\n", m_specInfo.m_interlaceStep);
        string.AppendFormat("\tstartchannel=%d\n", m_specInfo.m_startChannel);
        string.AppendFormat("\tspectrumlength=%d\n", 2048);
        string.AppendFormat("\tflux=%.2lf\n", scan.GetFlux());
        string.AppendFormat("\tbattery=%.2f\n", scan.GetBatteryVoltage());
        string.AppendFormat("\ttemperature=%.2f\n", scan.GetTemperature());
        // The mode
        if (scan.IsDirectSunMeasurement())
            string.AppendFormat("\tmode=direct_sun\n");
        else if (scan.IsLunarMeasurement())
            string.AppendFormat("\tmode=lunar\n");
        else if (scan.IsWindMeasurement())
            string.AppendFormat("\tmode=wind\n");
        else if (scan.IsStratosphereMeasurement())
            string.AppendFormat("\tmode=stratospheric\n");
        else if (scan.IsCompositionMeasurement())
            string.AppendFormat("\tmode=composition\n");
        else
            string.AppendFormat("\tmode=plume\n");

        double maxIntensity = CSpectrometerDatabase::GetInstance().GetModel(m_specInfo.m_specModelName).maximumIntensityForSingleReadout;

        // Finally, the version of the file and the version of the program
        string.AppendFormat("\tsoftwareversion=%d.%d\n", CVersion::majorNumber, CVersion::minorNumber);
        string.AppendFormat("\tcompiledate=%s\n", __DATE__);

        string.AppendFormat("</scaninformation>\n\n");
        fprintf(f, string);

        // ----------------- Create the flux-information ----------------------
        wind.GetWindSpeedSource(wsSrc);
        wind.GetWindDirectionSource(wdSrc);
        wind.GetPlumeHeightSource(phSrc);
        string.Format("<fluxinfo>\n");
        string.AppendFormat("\tflux=%.4lf\n", scan.GetFlux());
        string.AppendFormat("\twindspeed=%.4lf\n", wind.GetWindSpeed());
        string.AppendFormat("\twinddirection=%.4lf\n", wind.GetWindDirection());
        string.AppendFormat("\tplumeheight=%.2lf\n", wind.GetPlumeHeight());
        string.AppendFormat("\twindspeedsource=%s\n", (LPCSTR)wsSrc);
        string.AppendFormat("\twinddirectionsource=%s\n", (LPCSTR)wdSrc);
        string.AppendFormat("\tplumeheightsource=%s\n", (LPCSTR)phSrc);
        //if(fabs(spectrometer.m_scanner.compass) > 360.0)
        //	string.AppendFormat("\tcompasssource=compassreading\n");
        //else
        //	string.AppendFormat("\tcompasssource=user\n");
        string.AppendFormat("</fluxinfo>\n");
        fprintf(f, string);

        // ----------------------- write the header --------------------------------
        string.Format("#scanangle\tstarttime\tstoptime\tname\tspecsaturation\tfitsaturation\tdelta\tchisquare\texposuretime\tnumspec\t");

        for (int itSpecie = 0; itSpecie < scan.GetSpecieNum(0); ++itSpecie) {
            specieName.Format("%s", scan.GetSpecieName(0, itSpecie).c_str());
            string.AppendFormat("column(%s)\tcolumnerror(%s)\t", (LPCSTR)specieName, (LPCSTR)specieName);
            string.AppendFormat("shift(%s)\tshifterror(%s)\t", (LPCSTR)specieName, (LPCSTR)specieName);
            string.AppendFormat("squeeze(%s)\tsqueezeerror(%s)\t", (LPCSTR)specieName, (LPCSTR)specieName);
        }
        string.AppendFormat("isgoodpoint\toffset\tflag");
        string.AppendFormat("\n<spectraldata>\n");

        fprintf(f, string);


        // ------------------- Then write the parameters for each spectrum ---------------------------
        for (unsigned long itSpectrum = 0; itSpectrum < scan.GetEvaluatedNum(); ++itSpectrum) {
            // 3a. Pretty print the result and the spectral info into a string
            CEvaluationResult result;
            scan.GetResult(itSpectrum, result);

            FormatEvaluationResult(&scan.GetSpectrumInfo(itSpectrum), &result, 0.0, scan.GetSpecieNum(itSpectrum), string);

            // 3b. Write it to the evaluation log file
            fprintf(f, string);
            fprintf(f, "\n");
        }
        fprintf(f, "</spectraldata>\n");

    }



    // Remember to close the file
    fclose(f);
    */

    return 0;
}

void CScanEvaluationLogFileHandler::FormatEvaluationResult(
    const CSpectrumInfo* info,
    const CEvaluationResult* result,
    double maxIntensity,
    int nSpecies,
    std::string& destination)
{
    char tempBuffer[512];

    // Scan angle
    sprintf(tempBuffer, "%.0lf\t", info->m_scanAngle);
    destination = std::string(tempBuffer);

    // The start time
    sprintf(tempBuffer, "%02d:%02d:%02d\t", info->m_startTime.hour, info->m_startTime.minute, info->m_startTime.second);
    destination += std::string(tempBuffer);

    // 4. The stop time
    sprintf(tempBuffer, "%02d:%02d:%02d\t", info->m_stopTime.hour, info->m_stopTime.minute, info->m_stopTime.second);
    destination += std::string(tempBuffer);

    // 5 The name of the spectrum
    destination += SimplifyString(info->m_name) + std::string("\t");

    // 6. The (maximum) saturation ratio of the whole spectrum,
    //			the (maximum) saturation ratio in the fit-region
    //			and the normalized maximum intensity of the whole spectrum
    if (maxIntensity > 0.0) {
        sprintf(tempBuffer, "%.2lf\t", info->m_peakIntensity / maxIntensity);
        destination += std::string(tempBuffer);
        sprintf(tempBuffer, "%.2lf\t", info->m_fitIntensity / maxIntensity);
        destination += std::string(tempBuffer);
    }
    else {
        sprintf(tempBuffer, "%.2lf\t", info->m_peakIntensity);
        destination += std::string(tempBuffer);
        sprintf(tempBuffer, "%.2lf\t", info->m_fitIntensity);
        destination += std::string(tempBuffer);
    }
    sprintf(tempBuffer, "%.2lf\t", (info->m_peakIntensity - info->m_offset) / info->m_exposureTime);
    destination += std::string(tempBuffer);

    // 7. The delta of the fit
    if (result != nullptr)
        sprintf(tempBuffer, "%.2e\t", result->m_delta);
    else
        sprintf(tempBuffer, "%.2e\t", 0.0);
    destination += std::string(tempBuffer);

    // 8. The chi-square of the fit
    if (result != nullptr)
        sprintf(tempBuffer, "%.2e\t", result->m_chiSquare);
    else
        sprintf(tempBuffer, "%.2e\t", 0.0);
    destination += std::string(tempBuffer);

    // 9. The exposure time and the number of spectra averaged
    sprintf(tempBuffer, "%ld\t%ld\t", info->m_exposureTime, info->m_numSpec);
    destination += std::string(tempBuffer);

    // 10. The column/column error for each specie
    for (int itSpecie = 0; itSpecie < nSpecies; ++itSpecie) {
        if (result != nullptr) {
            if ((fabs(result->m_referenceResult[itSpecie].m_column) > 5e-2) && (fabs(result->m_referenceResult[itSpecie].m_columnError) > 5e-2))
            {
                sprintf(tempBuffer, "%.2lf\t%.2lf\t", result->m_referenceResult[itSpecie].m_column, result->m_referenceResult[itSpecie].m_columnError);
            }
            else
            {
                sprintf(tempBuffer, "%.2e\t%.2e\t", result->m_referenceResult[itSpecie].m_column, result->m_referenceResult[itSpecie].m_columnError);
            }
            destination += std::string(tempBuffer);

            sprintf(tempBuffer, "%.2lf\t%.2lf\t", result->m_referenceResult[itSpecie].m_shift, result->m_referenceResult[itSpecie].m_shiftError);
            destination += std::string(tempBuffer);

            sprintf(tempBuffer, "%.2lf\t%.2lf\t", result->m_referenceResult[itSpecie].m_squeeze, result->m_referenceResult[itSpecie].m_squeezeError);
            destination += std::string(tempBuffer);
        }
        else {
            sprintf(tempBuffer, "0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t");
            destination += std::string(tempBuffer);
        }
    }

    // 11. The quality of the fit
    if (result != nullptr)
        sprintf(tempBuffer, "%d\t", result->IsOK());
    else
        sprintf(tempBuffer, "%d\t", 1);
    destination += std::string(tempBuffer);

    // 12. The offset
    sprintf(tempBuffer, "%.0lf\t", info->m_offset);
    destination += std::string(tempBuffer);

    // 13. The 'flag' in the spectra
    sprintf(tempBuffer, "%d\n", info->m_flag);
    destination += std::string(tempBuffer);
}
