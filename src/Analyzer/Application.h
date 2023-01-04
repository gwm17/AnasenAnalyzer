#ifndef APPLICATION_H
#define APPLICATION_H

#include "Analyzer.h"

#include <filesystem>
#include <vector>
#include <memory>

namespace AnasenAnalyzer {

    class Application
    {
    public:
        Application();
        ~Application();

        bool IsValid() const { return m_isValid; }
        void LoadConfig(const std::filesystem::path& path);
        void Run();

    private:
        std::vector<std::filesystem::path> m_inputPaths;
        std::filesystem::path m_outputPath;
        std::unique_ptr<Analyzer> m_analyzerHandle;
        bool m_isValid;
    };
}

#endif