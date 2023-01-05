#include "Application.h"
#include "Track.h"

#include <fstream>
#include <string>

#include "fmt/format.h"
#include "fmt/core.h"
#include "fmt/std.h"

#include "TFile.h"
#include "TH1.h"

namespace AnasenAnalyzer {

    Application::Application() :
        m_analyzerHandle(nullptr), m_isValid(false)
    {
        TH1::AddDirectory(kFALSE);
        if(!EnforceDictionaryLinked())
        {
            fmt::print("Dictionary error!");
        }
    }

    Application::~Application() {}

    void Application::LoadConfig(const std::filesystem::path& path)
    {
        //Reset to default state (handles case where analyzer is reused)
        m_isValid = false;
        m_analyzerHandle.reset(nullptr);
        m_inputPaths.clear();
        
        std::ifstream config(path);
        if(!config.is_open())
        {
            fmt::print("Could not open config file {0} at Application::LoadConfig\n", path);
            return;
        }

        std::string junk, name;
        double minVal, maxVal;
        double density;
        CutHandler::Ref cuts = std::make_shared<CutHandler>();
        GateHandler::Ref gates = std::make_shared<GateHandler>();

        uint32_t z, a;
        int s;
        std::vector<uint32_t> zt, at;
        std::vector<int> st;

        AnalyzerFlags flags = AnalyzerFlags_None;

        fmt::print("Reading configuration from file {0}...\n", path);

        config >> junk;
        if (junk != "begin_target")
        {
            fmt::print("Configuration file not formated correctly! Must start with target information.\n");
            return;
        }
        config >> junk >> density;
        config >> junk;
        if(junk != "begin_elements")
        {
            fmt::print("Configuration file not formated correctly! Target field must contain a list of elements.\n");
            return;
        }
        config >> junk >> junk >> junk;
        while(true)
        {
            config >> junk;
            if(junk == "end_elements" || junk != "element")
                break;
            config >> z >> a >> s;
            zt.push_back(z);
            at.push_back(a);
            st.push_back(s);
        }
        config >> junk;
        config >> junk;
        if(junk != "begin_flags")
        {
            fmt::print("Configuration file not formated correctly! Must have list of analyzer flags.\n");
            return;
        }
        while(true)
        {
            config >> junk;
            if(junk == "end_flags")
                break;
            flags |= ConvertStringToFlag(junk);
        }
        config >> junk;
        if (junk != "begin_cuts")
        {
            fmt::print("Configuration file not formated correctly! Must have list of cut files!\n");
            return;
        }
        while (true)
        {
            config >> junk;
            if (junk == "end_cuts")
                break;
            config >> name;
            cuts->AddCut(junk, name);
        }

        config >> junk;
        if (junk != "begin_gates")
        {
            fmt::print("Input file not formated correctly! Must have list of gates files!\n");
            return;
        }
        while (true)
        {
            config >> junk;
            if (junk == "end_gates")
                break;
            config >> minVal >> maxVal;
            gates->AddGate({junk, minVal, maxVal});
        }

        config >> junk;
        if (junk != "begin_data")
        {
            fmt::print("Input file not formated correctly! Must end with list of data files!\n");
            return;
        }
        while (true)
        {
            config >> junk;
            if (junk == "end_data")
                break;
            m_inputPaths.push_back(junk);
        }
        config >> junk >> m_outputPath;

        m_analyzerHandle = std::make_unique<Analyzer>(cuts, gates, Target(zt, at, st, density));
        m_analyzerHandle->SetFlags(flags);

        m_isValid = true;
        fmt::print("Output path: {0}\n", m_outputPath);
        fmt::print("Parsed config in {0} successfully...\n", path);
    }

    void Application::Run()
    {
        TFile* output = TFile::Open(m_outputPath.c_str(), "RECREATE");
        if(!output || !output->IsOpen())
        {
            fmt::print("Could not open/create output file {0} at Application::Run.\n", m_outputPath);
            return;
        }

        for(auto& path : m_inputPaths)
        {
            m_analyzerHandle->Analyze(path);
        }

        output->cd();
        auto& plots = m_analyzerHandle->GetPlots();
        for(auto& plot : plots)
        {
            plot.second->Write(plot.second->GetName(), TObject::kOverwrite);
        }
        auto& cuts = m_analyzerHandle->GetCutHandler();
        for (auto& iter : cuts->GetMap())
        {
            iter.second->Write(iter.second->GetName(), TObject::kOverwrite);
        }
        output->Close();
        delete output;
    }
}