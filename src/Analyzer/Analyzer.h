#ifndef ANALYZER_H
#define ANALYZER_H

#include "Dict/Track.h"
#include "Dict/SiHit.h"
#include "Dict/PCHit.h"
#include "CutHandler.h"
#include "GateHandler.h"
#include "Target.h"

#include <unordered_map>
#include <filesystem>
#include <string>
#include <memory>
#include "TObject.h"


namespace AnasenAnalyzer {

    //Cute little way to name alias bit fields for flags. Make an enum of flags (not strongly typed)
    //and then alias int to a flag type.
    enum AnalyzerFlags_
    {
        AnalyzerFlags_None = 0,
        AnalyzerFlags_OneTrack = 1 << 1,
        AnalyzerFlags_TwoTrack = 1 << 2,
        AnalyzerFlags_ThreeTrack = 1 << 3,
        AnalyzerFlags_BinnedDalitz = 1 << 4
    };
    typedef int AnalyzerFlags;

    AnalyzerFlags ConvertStringToFlag(const std::string& keyword);

    class Analyzer
    {
    public:
        using PlotMap = std::unordered_map<std::string, std::shared_ptr<TObject>>;

        struct Histo1DParams
        {
            std::string name = "";
            std::string title = "";
            int nXBins = -1;
            double xMin = 0.0;
            double xMax = 0.0;
        };

        struct Histo2DParams
        {
            std::string name = "";
            std::string title = "";
            int nXBins = -1;
            double xMin = 0.0;
            double xMax = 0.0;
            int nYBins = -1;
            double yMin = 0.0;
            double yMax = 0.0;
        };

        Analyzer(const CutHandler::Ref& cuts, const GateHandler::Ref& gates, const Target& target);
        ~Analyzer();

        const AnalyzerFlags& GetFlags() const { return m_flags; }
        PlotMap& GetPlots() { return m_plotMap; }
        void SetFlags(AnalyzerFlags flags) { m_flags = flags; }
        
        void Analyze(const std::filesystem::path& path);

        const CutHandler::Ref& GetCutHandler() const { return m_cutHandler; }

    private:
        void FillHistogram1D(const Histo1DParams& params, double xValue);
        void FillHistogram2D(const Histo2DParams& params, double xValue, double yValue);
        int GetMaxPC(double phi, const PCHit& pcHits, const std::vector<bool>& pcMask);
        static constexpr double AngularDiff(double angle1, double angle2)
        {
            double diff = std::fabs(angle1 - angle2);
            return diff > M_PI ? 2.0 * M_PI - diff : diff;
        }

        TrackEvent CreateTrackEvent(const SiHit& siHits, const PCHit& pcHits);
        void PlotCorrelations(const TrackEvent& event);
        void AnalyzeOneTrack(const TrackEvent& event);
        void AnalyzeTwoTrack(const TrackEvent& event);
        void AnalyzeThreeTrack(const TrackEvent& event);

        PlotMap m_plotMap;
        CutHandler::Ref m_cutHandler;
        GateHandler::Ref m_gateHandler;
        Target m_target;
        AnalyzerFlags m_flags;

        static constexpr double s_epsilon = 1.0e-6;
        static constexpr double s_rad2deg = 180.0/M_PI;
        static constexpr double s_anasenLength = 55.0545; //cm
        static constexpr double s_initialBeamEnergy = 17.19; // MeV
        static constexpr int s_beamZ = 4;
        static constexpr int s_beamA = 7;
    };
}

#endif