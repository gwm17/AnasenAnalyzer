#ifndef CUT_HANDLER_H
#define CUT_HANDLER_H

#include <unordered_map>
#include <string>
#include <vector>
#include <memory>

#include "TFile.h"
#include "TCutG.h"

namespace AnasenAnalyzer {

    class CutHandler
    {
    public:
        using Ref = std::shared_ptr<CutHandler>;
        
        CutHandler();
        ~CutHandler();

        bool AddCut(const std::string &cutfile, const std::string &cutname);
        bool IsInside(const std::string &cutname, double valx, double valy);
        const std::unordered_map<std::string, TCutG*>& GetMap() const { return m_cutMap; }

    private:
        std::vector<TFile *> m_cutfiles;
        std::unordered_map<std::string, TCutG *> m_cutMap;
    };

}

#endif