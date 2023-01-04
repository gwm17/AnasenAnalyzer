#include "CutHandler.h"

namespace AnasenAnalyzer {

    CutHandler::CutHandler() {}

    CutHandler::~CutHandler()
    {
        for (std::size_t i = 0; i < m_cutfiles.size(); i++)
        {
            m_cutfiles[i]->Close();
            delete m_cutfiles[i];
        }
    }

    bool CutHandler::AddCut(const std::string &cutfile, const std::string &cutname)
    {
        TFile *file = TFile::Open(cutfile.c_str(), "READ");
        if (!file || !file->IsOpen())
            return false;

        TCutG *cut = (TCutG *)file->Get("CUTG");
        if (!cut)
            return false;

        cut->SetName(cutname.c_str());
        m_cutfiles.push_back(file);
        m_cutMap[cutname] = cut;
        return true;
    }

    bool CutHandler::IsInside(const std::string &cutname, double valx, double valy)
    {
        auto iter = m_cutMap.find(cutname);
        if (iter != m_cutMap.end())
        {
            return iter->second->IsInside(valx, valy);
        }

        return false;
    }

}