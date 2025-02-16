#pragma once

#include <string>

#include <SpectralEvaluation/Evaluation/FitWindow.h>

namespace novac
{
    /** A CFitWindowFileHandler object is capable of reading
        and writing fit windows from a .nfw (NovacFitWindow) file. */
    class CFitWindowFileHandler
    {
    public:

        /** Reads the fit windows which are defined in the provided file.
                @param fileName The name and path of the file to read from
                @return A vector with all the fit windows defined in the file.
                @return An empty vector if the reading failed. */
        std::vector<novac::CFitWindow> ReadFitWindowFile(const std::string& fileName);

        /** Writes the supplied fit-window to a file.
                @param window - the fit window to be written to file
                @param fileName - the name and path of the file to which to write
                @param overWrite - if true the file will be overwritten, if false, the file will be appended.
                @return true if the writing succeeded successfully. */
        bool WriteFitWindow(const novac::CFitWindow& window, const std::string& fileName, bool overWrite);
    };
}