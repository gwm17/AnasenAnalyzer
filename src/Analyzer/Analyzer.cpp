#include "Analyzer.h"
#include "Reconstructor.h"
#include "Precision.h"
#include "fmt/core.h"
#include "fmt/format.h"
#include "fmt/std.h"
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

namespace AnasenAnalyzer {

    AnalyzerFlags ConvertStringToFlag(const std::string& keyword)
    {
        if(keyword == "OneTrack")
            return AnalyzerFlags_OneTrack;
        else if (keyword == "TwoTrack")
            return AnalyzerFlags_TwoTrack;
        else if (keyword == "ThreeTrack")
            return AnalyzerFlags_ThreeTrack;
        else if (keyword == "BinnedDalitz")
            return AnalyzerFlags_BinnedDalitz;
        else
            return AnalyzerFlags_None;
    }

    Analyzer::Analyzer(const CutHandler::Ref& cuts, const GateHandler::Ref& gates, const Target& target) :
        m_cutHandler(cuts), m_gateHandler(gates), m_target(target), m_flags(AnalyzerFlags_None)
    {
    }

    Analyzer::~Analyzer() {}

    void Analyzer::FillHistogram1D(const Histo1DParams& params, double xValue)
    {
        auto iter = m_plotMap.find(params.name);
        if(iter == m_plotMap.end())
        {
            std::shared_ptr<TH1F> histo = std::make_shared<TH1F>(params.name.c_str(), params.title.c_str(), params.nXBins, params.xMin, params.xMax);
            histo->Fill(xValue);
            m_plotMap[params.name] = histo;
        }
        else
            std::static_pointer_cast<TH1>(iter->second)->Fill(xValue);
    }

    void Analyzer::FillHistogram2D(const Histo2DParams& params, double xValue, double yValue)
    {
        auto iter = m_plotMap.find(params.name);
        if(iter == m_plotMap.end())
        {
            std::shared_ptr<TH2F> histo = std::make_shared<TH2F>(params.name.c_str(), params.title.c_str(), params.nXBins, params.xMin, params.xMax, params.nYBins, params.yMin, params.yMax);
            histo->Fill(xValue, yValue);
            m_plotMap[params.name] = histo;
        }
        else
            std::static_pointer_cast<TH2>(iter->second)->Fill(xValue, yValue);
    }

    int Analyzer::GetMaxPC(double phi, const PCHit& pcHits, const std::vector<bool>& pcMask)
    {
        static constexpr double phiWindow = 0.5238; //30 deg maximum angular difference

        int maxPCIndex = -1;
        int nextToMaxPCIndex = -1;
        double maxPCEnergy = -10.0;

        // phi window

        // loop over hits
        for (std::size_t k = 0; k < (*(pcHits.ReadHit)).size(); k++)
        {
            auto &hit = (*(pcHits.ReadHit))[k];
            if (pcMask[k] && AngularDiff(hit.PhiW, phi) < phiWindow) // check to see if already used and within window
            {
                if (Precision::IsFloatGreaterOrAlmostEqual(hit.Energy, maxPCEnergy, s_epsilon))
                {
                    maxPCEnergy = hit.Energy;
                    maxPCIndex = k;
                }
            }
        }
        if (nextToMaxPCIndex > 0 && (AngularDiff((*pcHits.ReadHit)[maxPCIndex].PhiW, phi) > AngularDiff((*pcHits.ReadHit)[nextToMaxPCIndex].PhiW, phi)))
            return nextToMaxPCIndex;

        return maxPCIndex;
    }

    /*
        Take in all data and create track events. Ideally this would be handled by a previous calibration stage as many things can
        go wrong and we would like to have the data pre-filtered. But here we are.

        Note that there are three types of tracks, yet only one of these types actually involves any tracking.
        Complete -- A correlated PC hit and Si hit. These are actual tracks, and are the only thing that really gets analyzed.
        SiOnly -- A single Si hit. Can be useful, but generally don't have enough info for analysis
        PCOnly -- A single PC hit. Garbage.

        Generally only the complete are actually of any use, as tracking is the only way to deterimine the interaction point...

        NOTE: The tracking method here is EXTREMELY flawed. Due to a lack of information, it is not possible to define the exact interaction point.
        Currently this is approximated by assuming that the interaction point is on axis (along z), but we know that this is not true. RESOLUT ribs typically have
        a ~2cm diameter beam spot, and straggling through the active target will only worsen this. This impacts EVERYTHING; beam energy, energy conservation, event construction, etc.
        This can only be fixed with better detector engineering. This is related to the misconception that z-position determines beam energy; the beam energy is determined by the path-length
        through the gas volume. See diagram below. (The more "on axis" the ejecta are, the more severe the difference is, which is not great)
                                                _______Si
                                                   /
                                                  /
        -----------------------------------------/-----------------------------------------PC
        Real Path-Length                        /
        -------------------------------------->* real Interaction Point
                                              /
        Assumed Path-Length                  /
        ----------------------------------->* Assumed Interaction Point--------------------Z-axis
                                           /
                                          /
                                         / Trajectory of possible solutions
                                        /
        -----------------------------------------------------------------------------------PC
    */
    TrackEvent Analyzer::CreateTrackEvent(const SiHit& siHits, const PCHit& pcHits)
    {
        TrackEvent event;

        std::vector<bool> pcMask; //Set pc used (false) or un-used(true)
        pcMask.resize(pcHits.ReadHit->size(), true);
        for (std::size_t i=0; i<siHits.ReadHit->size(); i++)
        {
            auto& sihit = siHits.ReadHit->at(i);
            //Bad detector channels
            if (sihit.Energy <= 0)
                continue;
            else if (sihit.DetID == 7 || sihit.DetID == 9)
                continue;
            else if ((sihit.DetID == 11 || sihit.DetID == 26) &&
                     (sihit.BackChannel == 1 || sihit.BackChannel == 2 || sihit.BackChannel == 3))
                continue;
            else if (sihit.DetID == 12 && sihit.FrontChannel == 3)
                continue;
            else if (sihit.DetID == 18 && sihit.BackChannel == 2)
                continue;
            else if (sihit.DetID == 19 && sihit.BackChannel == 3)
                continue;
            else if (sihit.DetID == 22 && (sihit.FrontChannel == 2 || sihit.FrontChannel == 3))
                continue;
            else if (sihit.DetID == 24 && (sihit.BackChannel == 1 || sihit.BackChannel == 2))
                continue;
            else if (sihit.DetID == 25 && sihit.FrontChannel == 2)
                continue;
            else if (sihit.DetID == 27 && (sihit.BackChannel == 2 || sihit.BackChannel == 3))
                continue;

            //Find a correlated PC; if there isn't one, store as siOnly track
            int pcIndex = GetMaxPC(sihit.PhiW, pcHits, pcMask);
            if (pcIndex < 0)
            {
                Track track;
                track.type = TrackType::SiOnly;
                track.siHit = CreateTrackedSiHit(sihit);
                event.siOnly.push_back(track);
                continue;
            }
            auto& pchit = pcHits.ReadHit->at(pcIndex);
            if(pchit.WireID == 6)
            {
                Track track;
                track.type = TrackType::SiOnly;
                track.siHit = CreateTrackedSiHit(sihit);
                event.siOnly.push_back(track);
                continue;
            }

            Track track;
            track.type = TrackType::Complete;
            track.siHit = CreateTrackedSiHit(sihit);
            track.pcHit = CreateTrackedPCHit(pchit);

            //Do the actual tracking
            double m = (track.pcHit.r - track.siHit.r) / (track.pcHit.z - track.siHit.z);
            double b = (track.pcHit.r - m * track.pcHit.z);
            track.interactionPointZ = -b / m;
            track.theta = std::atan2(track.siHit.r, track.interactionPointZ - track.siHit.z);
            if (track.theta < 0.0)
                track.theta += M_PI;
            track.phi = track.siHit.phi;
            track.pathLength = track.siHit.r / std::sin(track.theta);
            if (track.interactionPointZ > 0.0 && track.interactionPointZ < s_anasenLength)
            {
                track.beamEnergyLoss = m_target.GetEnergyLoss(s_beamZ, s_beamA, s_initialBeamEnergy, (s_anasenLength - track.interactionPointZ) / 100.0);
                track.beamEnergy = s_initialBeamEnergy - track.beamEnergyLoss;
            }

            pcMask[i] = false;

            event.complete.push_back(track);
        }

        //Store all non-paired PC data as pcOnly
        for(std::size_t i=0; i<pcHits.ReadHit->size(); i++)
        {
            if(!pcMask[i])
                continue;

            Track track;
            track.type = TrackType::PCOnly;
            track.pcHit = CreateTrackedPCHit(pcHits.ReadHit->at(i));
            event.pcOnly.push_back(track);
        }

        std::sort(event.complete.begin(), event.complete.end(), TrackEvent::SortTracksBySi);
        std::sort(event.siOnly.begin(), event.siOnly.end(), TrackEvent::SortTracksBySi);
        std::sort(event.pcOnly.begin(), event.pcOnly.end(), TrackEvent::SortTracksByPC);

        return event;
    }

    void Analyzer::PlotCorrelations(const TrackEvent& event)
    {
        for (auto& track : event.complete)
        {
            FillHistogram1D({"InteractionPointZ", "InteractionPointZ; Z(cm);", 300, -10, 56}, track.interactionPointZ);
            FillHistogram1D({"BeamEnergy", "BeamEnergy; EBeam(MeV);", 100, -10.0, 90.0}, track.beamEnergy);
            FillHistogram2D({"BeamEnergy_vs_IntPointZ", "BeamEnergy_vs_IntPointZ; Interaction Z (cm); EBeam (MeV)", 200, -10.0, 60.0, 200, -10.0, 30.0}, track.interactionPointZ, track.beamEnergy);
            FillHistogram2D({"E_de", "EdE;Si E (MeV); PC E (arb)", 200, -1, 35, 200, -0.01, 1.5}, track.siHit.energy, track.pcHit.energy);
            FillHistogram2D({"E_de_corrAll", "EdE_correctedAll;Si E (MeV); PC adj. E (arb)", 200, -1, 35, 200, -0.01, 1.5}, track.siHit.energy, track.pcHit.energy * std::sin(track.theta));
            FillHistogram2D({"E_si_vs_Theta", "E_si_vs_trackTheta; Si E(MeV); Track #theta (deg)", 500, 0.0, 200.0, 500, 0.0, 35.0}, track.theta * s_rad2deg, track.siHit.energy);
            if (track.siHit.detID < 4 && track.siHit.detID > -1)
            {
                FillHistogram1D({"InteractionPointZ_QQQ", "InteractionPointZ_QQQ; Z(cm);", 300, -10, 56}, track.interactionPointZ);
                FillHistogram1D({"BeamEnergy_QQQ", "BeamEnergy_QQQ; EBeam(MeV);", 100, -10.0, 90.0}, track.beamEnergy);
                FillHistogram2D({"BeamEnergy_vs_IntPointZ_QQQ", "BeamEnergy_vs_IntPointZ_QQQ; Interaction Z (cm); EBeam (MeV)", 200, -10.0, 60.0, 200, -10.0, 30.0}, track.interactionPointZ, track.beamEnergy);
                FillHistogram2D({"E_de_QQQ", "EdE_QQQ;Si E (MeV); PC E (arb)", 200, -1, 35, 200, -0.01, 1.5}, track.siHit.energy, track.pcHit.energy);
                FillHistogram2D({"E_de_corrAll_QQQ", "EdE_correctedAll_QQQ;Si E (MeV); PC adj. E (arb)", 200, -1, 35, 200, -0.01, 1.5}, track.siHit.energy, track.pcHit.energy * std::sin(track.theta));
                FillHistogram2D({"E_si_vs_Theta_QQQ", "E_si_vs_trackTheta_QQQ; Si E(MeV); Track #theta (deg)", 500, 0.0, 200.0, 500, 0.0, 35.0}, track.theta * s_rad2deg, track.siHit.energy);
            }
            else if (track.siHit.detID < 16 && track.siHit.detID > 3)
            {
                FillHistogram1D({"InteractionPointZ_B1", "InteractionPointZ_B1; Z(cm);", 300, -10, 56}, track.interactionPointZ);
                FillHistogram1D({"BeamEnergy_B1", "BeamEnergy_B1; EBeam(MeV);", 100, -10.0, 90.0}, track.beamEnergy);
                FillHistogram2D({"BeamEnergy_vs_IntPointZ_B1", "BeamEnergy_vs_IntPointZ_B1; Interaction Z (cm); EBeam (MeV)", 200, -10.0, 60.0, 200, -10.0, 30.0}, track.interactionPointZ, track.beamEnergy);
                FillHistogram2D({"E_de_B1", "EdE_B1;Si E (MeV); PC E (arb)", 200, -1, 35, 200, -0.01, 1.5}, track.siHit.energy, track.pcHit.energy);
                FillHistogram2D({"E_de_corrAll_B1", "EdE_correctedAll_B1;Si E (MeV); PC adj. E (arb)", 200, -1, 35, 200, -0.01, 1.5}, track.siHit.energy, track.pcHit.energy * std::sin(track.theta));
                FillHistogram2D({"E_si_vs_Theta_B1", "E_si_vs_trackTheta_B1; Si E(MeV); Track #theta (deg)", 500, 0.0, 200.0, 500, 0.0, 35.0}, track.theta * s_rad2deg, track.siHit.energy);
            }
            else if (track.siHit.detID < 28 && track.siHit.detID > 15)
            {
                FillHistogram1D({"InteractionPointZ_B2", "InteractionPointZ_B2; Z(cm);", 300, -10, 56}, track.interactionPointZ);
                FillHistogram1D({"BeamEnergy_B2", "BeamEnergy_B2; EBeam(MeV);", 100, -10.0, 90.0}, track.beamEnergy);
                FillHistogram2D({"BeamEnergy_vs_IntPointZ_B2", "BeamEnergy_vs_IntPointZ_B2; Interaction Z (cm); EBeam (MeV)", 200, -10.0, 60.0, 200, -10.0, 30.0}, track.interactionPointZ, track.beamEnergy);
                FillHistogram2D({"E_de_B2", "EdE_B2;Si E (MeV); PC E (arb)", 200, -1, 35, 200, -0.01, 1.5}, track.siHit.energy, track.pcHit.energy);
                FillHistogram2D({"E_de_corrAll_B2", "EdE_correctedAll_B2;Si E (MeV); PC adj. E (arb)", 200, -1, 35, 200, -0.01, 1.5}, track.siHit.energy, track.pcHit.energy * std::sin(track.theta));
                FillHistogram2D({"E_si_vs_Theta_B2", "E_si_vs_trackTheta_B2; Si E(MeV); Track #theta (deg)", 500, 0.0, 200.0, 500, 0.0, 35.0}, track.theta * s_rad2deg, track.siHit.energy);
            }

            if (m_cutHandler->IsInside("protonEde", track.siHit.energy, track.pcHit.energy * std::sin(track.theta)) && track.IsBeamValid(s_initialBeamEnergy))
            {
                FillHistogram1D({"InteractionPointZ_protonCut", "InteractionPointZ_protonCut; Z(cm);", 300, -10, 56}, track.interactionPointZ);
                FillHistogram2D({"IntPointZ_SiE_protonCut", "IntPointZ_SiE_protonCut;Interaction Z (cm); Si E(MeV)", 550, 0.0, 55.0, 150, 0, 15.0}, track.interactionPointZ, track.siHit.energy);
            }
            else if (m_cutHandler->IsInside("alphaEde", track.siHit.energy, track.pcHit.energy * std::sin(track.theta)) && track.IsBeamValid(s_initialBeamEnergy))
            {
                FillHistogram1D({"InteractionPointZ_alphaCut", "InteractionPointZ_alphaCut; Z(cm);", 300, -10, 56}, track.interactionPointZ);
                FillHistogram2D({"IntPointZ_SiE_alphaCut", "IntPointZ_SiE_alphaCut;Interaction Z (cm); Si E(MeV)", 550, 0.0, 55.0, 150, 0, 15.0}, track.interactionPointZ, track.siHit.energy);
            }
        }
    }

    void Analyzer::AnalyzeOneTrack(const TrackEvent& event)
    {
        static Reconstructor reconBe8(m_target, {1, 4, 1}, {2, 7, 1});
        static Reconstructor reconLi6(m_target, {1, 4, 2}, {2, 7, 3});

        ReconResult result;
        for (auto& track : event.complete)
        {
            // All histos in this section tagged with _cre (CalcRecoilE) to avoid confusion
            // in rootfile
            if (track.IsValid() && m_cutHandler->IsInside("protonEde", track.siHit.energy, track.pcHit.energy * std::sin(track.theta)))
            {
                result = reconBe8.CalcResidualTracked(track);
                FillHistogram1D({"Ex_8Be_singleTrack", "Ex_8Be_singleTrack;E_{ex} MeV;", 600, -10.0, 40.0}, result.residualExcitation);
                FillHistogram1D({"BeamKE_8Be_singleTrack", "BeamKE_8Be_singleTrack; BeamKE(MeV);", 600, -10.0, 80.0}, result.projectileVec.E() - result.projectileVec.M());
                FillHistogram1D({"BeamKE_8Be_cm_singleTrack", "Ecm_8Be_singleTrack; Ecm(MeV);", 600, -1.0, 15.0}, result.parentECM);
                FillHistogram1D({"ProtonKErxn_8Be_singleTrack", "ProtonKErxn_8Be_singleTrack; KEProton (MeV);", 600, -10.0, 20.0}, result.ejectileRxnKE);
                FillHistogram2D({"Ex8Be_vs_BeamKE_tracked_singleTrack", "Ex8Be_vs_BeamKE_singleTrack;BeamKE (MeV);Ex(MeV)", 600, -10.0, 80.0, 400, -10.0, 20.0}, result.projectileRxnKE, result.residualExcitation);
                FillHistogram2D({"ProtonKErxn_vs_Ex8Be_tracked_singleTrack", "ProtonKErxn_vs_Ex8Be_singleTrack;KEProton(MeV);Ex(MeV)", 600, -10.0, 60.0, 400, -10.0, 20.0}, result.ejectileRxnKE, result.residualExcitation);
                FillHistogram2D({"ProtonKErxn_vs_ProtonTheta_8Be_singleTrack", "ProtonKErxn_vs_ProtonTheta_8Be_singleTrack;KEProton(MeV);#theta (deg)", 600, 0.0, 180.0, 600, 0.0, 20.0}, track.theta * s_rad2deg, result.ejectileRxnKE);
            }
        }
    }

    void Analyzer::AnalyzeTwoTrack(const TrackEvent& event)
    {
        static Reconstructor reconBe8(m_target, {1, 4, 1, 2}, {2, 7, 1, 4});
        static Reconstructor reconLi5(m_target, {1, 4, 2, 1}, {2, 7, 4, 1});

        if(event.complete.size() != 2)
            return;

        ReconResult result;
        std::vector<Track> alphas;
        std::vector<Track> protons;
        for(auto& track : event.complete)
        {
            if(track.IsValid() && track.IsBeamValid(s_initialBeamEnergy))
            {
                if(m_cutHandler->IsInside("alphaEde", track.siHit.energy, track.pcHit.energy * std::sin(track.theta)))
                {
                    alphas.push_back(track);
                }
                else if (m_cutHandler->IsInside("protonEde", track.siHit.energy, track.pcHit.energy * std::sin(track.theta)))
                {
                    protons.push_back(track);
                }
            }
        }

        if(protons.size() == 1 && alphas.size() == 1)
        {
            result = reconLi5.CalcTwoParticle(protons[0], alphas[0]);
            FillHistogram1D({"Ex5Li_reconTwoParticle", "Ex5Li_reconTwoParticle;Ex(MeV);", 600, -5.0, 25.0}, result.residualExcitation);
            FillHistogram2D({"Ex5Li_beamKECM_reconTwoParticle", "Ex5Li_beamKECM_reconTwoParticle;Ex(MeV);beamKE(MeV)", 600, -5.0, 25.0, 800, -10.0, 30.0}, result.residualExcitation, result.projectileRxnKE);
            FillHistogram2D({"Ex5Li_intPoint_reconTwoParticle", "Ex5Li_intPoint_reconTwoParticle;Ex(MeV);intPointZ(cm)", 600, -5.0, 25.0, 800, -10.0, 60.0}, result.residualExcitation, result.avgInteractionZ);
        }
        else if(alphas.size() == 2)
        {
            result = reconBe8.CalcTwoParticle(alphas[0], alphas[1]);
            FillHistogram1D({"Ex8Be_reconTwoParticle", "Ex8Be_reconTwoParticle;Ex(MeV);", 600, -5.0, 25.0}, result.residualExcitation);
            FillHistogram2D({"Ex8Be_beamKECM_reconTwoParticle", "Ex8Be_beamKECM_reconTwoParticle;Ex(MeV);beamKE(MeV)", 600, -5.0, 25.0, 800, -10.0, 30.0}, result.residualExcitation, result.projectileRxnKE);
            FillHistogram2D({"Ex8Be_intPoint_reconTwoParticle", "Ex8Be_intPoint_reconTwoParticle;Ex(MeV);intPointZ(cm)", 600, -5.0, 25.0, 800, -10.0, 60.0}, result.residualExcitation, result.avgInteractionZ);
            if(m_gateHandler->IsInside("8Begs", result.residualExcitation))
            {
                FillHistogram1D({"Ecm_be8gs_reconTwoParticle", "Ecm_be8gs_reconTwoParticle;Ex(MeV);", 1600, -10.0, 30.0}, result.parentECM);
                FillHistogram1D({"Ecm_be8gs_reconTwoParticle_goodBins", "Ecm_be8gs_reconTwoParticle;Ex(MeV);", 60, 0.0, 2.7}, result.parentECM);
            }
            else if(m_gateHandler->IsInside("8Be1ex", result.residualExcitation))
            {
                FillHistogram1D({"Ecm_be81ex_reconTwoParticle", "Ecm_be81ex_reconTwoParticle;Ex(MeV);", 1600, -10.0, 30.0}, result.parentECM);
                FillHistogram1D({"Ecm_be81ex_reconTwoParticle_goodBins", "Ecm_be81ex_reconTwoParticle;Ex(MeV);", 60, 0.0, 2.7}, result.parentECM);
            }
        }
    }

    void Analyzer::AnalyzeThreeTrack(const TrackEvent& event)
    {
        static Reconstructor reconBe8(m_target, {1, 4, 1, 2}, {2, 7, 1, 4});
        static Reconstructor reconLi5(m_target, {1, 4, 2, 1}, {2, 7, 4, 1});

        if(event.complete.size() != 3)
            return;

        ReconResult result8BeTracked, result8BeQ, result5LiTracked, result5LiQ;
        std::vector<Track> alphas;
        std::vector<Track> protons;
        for(auto& track : event.complete)
        {
            if(track.IsValid() /*&& track.IsBeamValid(s_initialBeamEnergy)*/)
            {
                if(m_cutHandler->IsInside("alphaEde", track.siHit.energy, track.pcHit.energy * std::sin(track.theta)))
                {
                    alphas.push_back(track);
                }
                else if (m_cutHandler->IsInside("protonEde", track.siHit.energy, track.pcHit.energy * std::sin(track.theta)))
                {
                    protons.push_back(track);
                }
            }
        }

        if (protons.size() == 1 && alphas.size() == 2)
        {
            double deltaPhi = std::abs(alphas[0].phi - alphas[1].phi);
            double deltaTheta = std::abs(alphas[0].theta - alphas[1].theta);

            FillHistogram1D({"protonIntPointZ_3track", "protonIntPointZ_3track;IntPointZ(cm);", 550, 0.0, 55.0}, protons[0].interactionPointZ);
            FillHistogram1D({"E_cm_a1_3track", "E_cm_a1_3track;Ecm(MeV);", 180, -1, 15}, alphas[0].beamEnergy * 4 / 22);
            FillHistogram1D({"E_cm_a2_3track", "E_cm_a2_3track;Ecm(MeV);", 180, -1, 15}, alphas[1].beamEnergy * 4 / 22);
            FillHistogram2D({"proton_IntPoint_SiEnergy_3track", "protonIntPoint_SiEnergy_3track_;IntPointZ(cm);Si E(MeV)", 550, 0.0, 55.0, 150, 0.0, 15.0}, protons[0].interactionPointZ, protons[0].siHit.energy);
            FillHistogram2D({"IntPoint_a1_vs_IntPoint_a2_3track", "IntPoint_a1_vs_IntPoint_a2_3track;IntPointZ(cm);IntPointZ(cm)", 400, -10.0, 70.0, 400, -10.0, 70.0}, alphas[0].interactionPointZ, alphas[1].interactionPointZ);
            FillHistogram2D({"IntPoint_p_vs_IntPoint_a1_3track", "IntPoint_p_vs_IntPoint_a1_3track;IntPointZ(cm);IntPointZ(cm)", 400, -10.0, 70.0, 400, -10.0, 70.0}, protons[0].interactionPointZ, alphas[0].interactionPointZ);
            FillHistogram2D({"IntPoint_p_vs_IntPoint_a2_3track", "IntPoint_p_vs_IntPoint_a2_3track;IntPointZ(cm);IntPointZ(cm)", 400, -10.0, 70.0, 400, -10.0, 70.0}, protons[0].interactionPointZ, alphas[1].interactionPointZ);
            FillHistogram2D({"Theta_a1_vs_Theta_a2_3track", "Theta_a1_vs_Theta_a2_3track;#theta (deg); #theta (deg)", 400, -1.0, 190.0, 400, -1, 190.0}, alphas[0].theta * s_rad2deg, alphas[1].theta * s_rad2deg);
            FillHistogram2D({"SiPhi_a1_vs_SiPhi_a2_3track", "SiPhi_a1_vs_SiPhi_a2_3track; #phi (deg); #phi (deg)", 400, -1.0, 360.0, 400, -1.0, 360.0}, alphas[0].phi * s_rad2deg, alphas[1].phi * s_rad2deg);
            FillHistogram2D({"SiE_a1_vs_SiE_a2_3track", "SiE_a1_vs_SiE_a2_3track;Si E(MeV); Si E(MeV);", 400, -1.0, 20.0, 400, -1.0, 20.0}, alphas[0].siHit.energy, alphas[1].siHit.energy);
            FillHistogram2D({"BeamE_a1_vs_BeamE_a2_3track", "BeamE_a1_vs_BeamE_a2_3track;BeamE(MeV);BeamE(MeV);", 400, -1.0, 80.0, 400, -1.0, 80.0}, alphas[0].beamEnergy, alphas[1].beamEnergy);
            FillHistogram2D({"PCZ_a1_vs_PCZ_a2_3track", "PCZ_a1_vs_PCZ_a2_3track;PCZ(cm);PCZ(cm);", 600, -10.0, 60.0, 600, -1.0, 60.0}, alphas[0].pcHit.z, alphas[1].pcHit.z);

            result8BeTracked = reconBe8.CalcMultiParticleTracked(protons[0], alphas[0], alphas[1]);
            result8BeQ = reconBe8.CalcMultiParticleQvalue(protons[0], alphas[0], alphas[1]);
            auto temp1 = reconLi5.CalcMultiParticleQvalue(alphas[1], protons[0], alphas[0]);
            auto temp2 = reconLi5.CalcMultiParticleQvalue(alphas[0], protons[0], alphas[1]);
            if(temp1.residualExcitation < temp2.residualExcitation)
            {
                result5LiQ = temp1;
                result5LiTracked = reconLi5.CalcMultiParticleTracked(alphas[1], protons[0], alphas[0]);
            }
            else
            {
                result5LiQ = temp2;
                result5LiTracked = reconLi5.CalcMultiParticleTracked(alphas[0], protons[0], alphas[1]);
            }

            FillHistogram1D({"Ex8Be_Qval_3track", "Ex8Be_Qval_3track; Ex(MeV);", 600, -20.0, 50.0}, result8BeQ.residualExcitation);
            FillHistogram1D({"Ex8Be_Track_3track", "Ex8Be_Track_3track; Ex(MeV);", 600, -20.0, 50.0}, result8BeTracked.residualExcitation);
            FillHistogram1D({"Ex9B_Qval_3track", "Ex9B_Qval_3track;Ex(MeV);", 600, 0.0, 30.0}, result8BeQ.parentExcitation);
            FillHistogram1D({"Ex5Li_Qval_3track", "Ex5Li_Qval_3track; Ex(MeV);", 600, -20.0, 50.0}, result5LiQ.residualExcitation);
            FillHistogram1D({"Ex5Li_Track_3track", "Ex5Li_Track_3track; Ex(MeV);", 600, -20.0, 50.0}, result5LiTracked.residualExcitation);

            FillHistogram2D({"Ex8Be_pKE_Qval_3track", "Ex8Be_pKE_Qval_3track;Ex(MeV);KE(MeV)", 600, -20, 50, 600, -10, 20.0}, result8BeQ.residualExcitation, result8BeQ.ejectileRxnKE);
            FillHistogram2D({"BeamKE_Track_vs_BeamKE_Qval_3track", "BeamKE_Track_vs_BeamKE_Qval_3track;KE Track(MeV);KE Q(MeV)", 600, -10, 80, 600, -10, 80}, result8BeTracked.projectileRxnKE, result8BeQ.projectileRxnKE);
            FillHistogram2D({"Ex8Be_vs_BeamKE_Qval_3track", "Ex8Be_vs_BeamKE_Qval_3track;Ex(MeV);KE(MeV)",600, -30, 80, 600, -20, 50}, result8BeQ.projectileRxnKE, result8BeQ.residualExcitation);
            FillHistogram2D({"Ex8BeTrack_vs_Ex8BeQval_3track", "Ex8BeTrack_vs_Ex8BeQval_3track;Ex Track(MeV);Ex Q(MeV)", 600, -20, 50, 600, -20, 50}, result8BeTracked.residualExcitation, result8BeQ.residualExcitation);
            FillHistogram2D({"Beamxy_Qval_3track", "Beamxy_Qval_3track", 360, -180.0, 180.0, 360, -180.0, 180.0}, result8BeQ.projectileVec.Px(), result8BeQ.projectileVec.Py());
            FillHistogram2D({"AvgIntPoint_vs_BeamKE_8BeQval_3track", "AvgIntPoint_vs_BeamKE_8BeQval_3track;IntZ(cm);KE(MeV)",150, 0.0, 50.0, 150, 0, 15.0}, result8BeQ.avgInteractionZ, result8BeTracked.projectileRxnKE);

            FillHistogram2D({"Ex5Li_ejectKE_Qval_3track", "Ex5Li_ejectKE_Qval_3track;Ex(MeV);KE(MeV)", 600, -20.0, 50.0, 600, -10.0, 20.0}, result5LiQ.residualExcitation, result5LiQ.ejectileRxnKE);
            FillHistogram2D({"Ex5LiTracked_vs_Ex5LiQval_3track", "Ex5LiTracked_vs_Ex5LiQval_3track; ExTrack(MeV);ExQ(MeV)", 600, -20.0, 50.0, 600, -20.0, 50.0}, result5LiTracked.residualExcitation, result5LiQ.residualExcitation);
            FillHistogram2D({"Ex5Li_vs_BeamKE_Qval_3track", "Ex5Li_vs_BeamKE_Qval_3track;KE(MeV);Ex(MeV)", 600, -10.0, 80.0, 600, -20.0, 50.0}, result5LiQ.projectileRxnKE, result5LiQ.residualExcitation);
            
            double beamKE_step = 0.2;
            int beamBin = std::ceil(result8BeQ.projectileRxnKE / beamKE_step);
            int tlbin = std::ceil(std::cos(result5LiQ.ejectileThetaCM) / 0.2) + 5;

            if (result8BeQ.projectileRxnKE > 0.0)
            {
                FillHistogram2D({"Dalitz", "Dalitz;a+a (GeV^{2});p+a(GeV^{2});", 90, 55.55, 56.0, 50, 21.75, 22.0}, result8BeQ.residualMassSq, result5LiQ.residualMassSq);
                FillHistogram2D({"Dalitz_Ex", "Dalitz_Ex;a+a(MeV);p+a(MeV)", 100, -10.0, 20.0, 100, -10.0, 20.0}, result8BeQ.residualExcitation, result5LiQ.residualExcitation);
                if (m_cutHandler->IsInside("li5Dalitz", result8BeQ.residualMassSq, result5LiQ.residualMassSq))
                {
                    if (result5LiQ.parentECM > 0.267 && result5LiQ.parentECM < 0.468)
                        FillHistogram1D({"AngDist_Li5_resGated", "AngDist_Li5_resGated;cos#theta_{CM}", 10, -1.0, 1.0}, std::cos(result5LiQ.ejectileThetaCM));
                }

                if (m_cutHandler->IsInside("be8Dalitz", result8BeQ.residualMassSq, result5LiQ.residualMassSq))
                {
                    if (result8BeQ.parentECM > 0.267 && result8BeQ.parentECM < 0.468)
                        FillHistogram1D({"AngDist_Be8_resGated", "AngDist_Be8_resGated;cos#theta_{CM}", 10, -1.0, 1.0}, std::cos(result8BeQ.ejectileThetaCM));
                }

                if (m_flags & AnalyzerFlags_BinnedDalitz)
                {
                    FillHistogram2D({fmt::format("Dalitz_{0}_beamBin", beamBin), fmt::format("Dalitz_{0}_beamBin;a+a(GeV^2);a+p(GeV^2)", beamBin), 90, 55.55, 56.0, 50, 21.75, 22.0},
                                    result8BeQ.residualMassSq, result5LiQ.residualMassSq);
                    if (m_cutHandler->IsInside("li5Dalitz", result8BeQ.residualMassSq, result5LiQ.residualMassSq))
                    {

                        FillHistogram2D({fmt::format("Dalitz_{0}_beamBin_{1}_thetaBin_5LigsGated", beamBin, tlbin), fmt::format("Dalitz_{0}_beamBin_{1}_thetaBin_5LigsGated;a+a(GeV^2);a+p(GeV^2)", beamBin, tlbin),
                                         90, 55.55, 56.0, 50, 21.75, 22.0}, result8BeQ.residualMassSq, result5LiQ.residualMassSq);
                    }
                }
            }
        }
    }

    void Analyzer::Analyze(const std::filesystem::path& path)
    {
        if(!std::filesystem::exists(path))
        {
            fmt::print("Could not open file {0} at Analyzer::Analyze\n", path);
            return;
        }

        TFile* inputFile = TFile::Open(path.c_str(), "READ");
        if(!inputFile || !inputFile->IsOpen())
        {
            fmt::print("ROOT could not open file {0} at Analyzer::Analyze\n", path);
            return;
        }

        TTree* inputTree = (TTree*) inputFile->Get("MainTree");
        if(!inputTree)
        {
            fmt::print("File {0} did not contain tree MainTree\n", path);
            inputFile->Close();
            delete inputFile;
        }

        /*
			The gross ReadHit, Hit convention is based off the selective saving 
			of only parts of the class SiHit, PCHit, CsIHit, etc. to the tree rather than
			saving the whole class, even though the whole class is loaded into the dictionary!
			To make matters even worse, they then doubled down and use the class as the reader as well.
			So you do HAVE to use ReadHit to read, even though the structure makes no sense.
		*/

        SiHit siHits;
        PCHit pcHits;
        //Some wierd stuff here... set vectors to null?
        siHits.ReadDet = nullptr;
        siHits.ReadHit = nullptr;
        pcHits.ReadHit = nullptr;

		inputTree->SetBranchAddress("Si.NSiHits", &siHits.NSiHits);
		inputTree->SetBranchAddress("Si.Detector", &siHits.ReadDet);
		inputTree->SetBranchAddress("Si.Hit", &siHits.ReadHit);
		inputTree->SetBranchAddress("PC.NPCHits", &pcHits.NPCHits);
		inputTree->SetBranchAddress("PC.Hit", &pcHits.ReadHit);

        uint64_t nentries = inputTree->GetEntries();
        constexpr double flushPercent = 0.05;
        const uint64_t flushVal = nentries * flushPercent;
        uint64_t flushCount = 0;
        uint64_t count = 0;

        fmt::print("Starting processing file {0}...\n", path);
        TrackEvent event;
        for(uint64_t i=0; i<nentries; i++)
        {
            inputTree->GetEntry(i);
            count++;
            if(count == flushVal)
            {
                count = 0;
                flushCount++;
                std::cout << fmt::format("\rPercent of data processed: {0:.3}% ", flushCount * flushPercent * 100) << std::flush;
            }

            event = CreateTrackEvent(siHits, pcHits);

            PlotCorrelations(event);

            if(m_flags & AnalyzerFlags_OneTrack)
            {
                AnalyzeOneTrack(event);
            }
            if(m_flags & AnalyzerFlags_TwoTrack)
            {
                AnalyzeTwoTrack(event);
            }
            if(m_flags & AnalyzerFlags_ThreeTrack)
            {
                AnalyzeThreeTrack(event);
            }
        }
        fmt::print("\nProcessing complete.\n");

        inputFile->Close();
        delete inputFile;
    }
}