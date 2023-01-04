#include "Reconstructor.h"
#include "MassLookup.h"

#include "fmt/core.h"
#include "fmt/format.h"
#include "Math/Boost.h"

namespace AnasenAnalyzer {

    // Constructor
    Reconstructor::Reconstructor(const Target& target, const std::vector<uint32_t>& Z, const std::vector<uint32_t>& A) : 
        m_isInit(false), m_target(target)
    {
        Init(Z, A);
    }

    // Destructor
    Reconstructor::~Reconstructor() {}

    void Reconstructor::Init(const std::vector<uint32_t>& Z, const std::vector<uint32_t>& A)
    {
        if (Z.size() != A.size() || Z.size() < 3)
        {
            m_isInit = false;
            return;
        }

        MassLookup &mass = MassLookup::GetInstance();
        int zr = Z[0] + Z[1] - Z[2];
        int ar = A[0] + A[1] - A[2];
        if (zr < 0 || ar <= 0 || zr > ar)
        {
            m_isInit = false;
            return;
        }
        m_targetData = mass.FindData(Z[0], A[0]);
        m_projectileData = mass.FindData(Z[1], A[1]);
        m_ejectileData = mass.FindData(Z[2], A[2]);
        m_residualData = mass.FindData(zr, ar);
        if (MassLookup::IsInvalidMass(m_targetData.mass) || MassLookup::IsInvalidMass(m_projectileData.mass) ||
            MassLookup::IsInvalidMass(m_ejectileData.mass) || MassLookup::IsInvalidMass(m_residualData.mass))
        {
            m_isInit = false;
            return;
        }

        if (Z.size() == 4)
        {
            int zb2 = zr - Z[3];
            int ab2 = ar - A[3];
            m_breakup1Data = mass.FindData(Z[3], A[3]);
            m_breakup2Data = mass.FindData(zb2, ab2);
            if (MassLookup::IsInvalidMass(m_breakup1Data.mass) || MassLookup::IsInvalidMass(m_breakup2Data.mass))
            {
                m_isInit = false;
                return;
            }
        }

        m_isInit = true;
    }

    /*DeltaZ()
     *Caclulates uncertainties for interaction point average. In principle not really necessary for
     *7Be+d since alphas & 3He are pretty easy, but good to keep just in case
     */
    double Reconstructor::DeltaZ(double Zpc, double Zsi, double Rpc, double Rsi)
    {
        double dZsi = -(Rpc / (Rsi - Rpc));
        double dZpc = 1 - dZsi;
        double dRsi = (Rpc * (Zsi - Zpc)) / ((Rsi - Rpc) * (Rsi - Rpc));
        double dRpc = -(Rsi * (Zsi - Zpc)) / ((Rsi - Rpc) * (Rsi - Rpc));

        double dz = std::sqrt(dZpc * dZpc + dRpc * dRpc + dZsi * dZsi * 0.1 * 0.1 + dRsi * dRsi * 0.1 * 0.1);
        return dz;
    }

    /*WIntP()
     *Caclulates weighted average interaction point
     *Weights are from the uncertainty calculation DeltaZ
     */
    double Reconstructor::WeightedAverageIntP(const std::vector<Track>& tracks)
    {
        std::vector<double> dz;
        std::vector<double> intp;
        double denom = 0;
        double IntPoint = 0;
        for (auto &track : tracks)
        {
            intp.push_back(track.interactionPointZ);
            double dz_i = DeltaZ(track.pcHit.z, track.siHit.z, track.pcHit.r, track.siHit.r);
            dz.push_back(dz_i);
            denom += 1.0 / dz_i;
        }
        for (size_t i = 0; i < tracks.size(); i++)
            IntPoint += (intp[i] / dz[i]) / denom;
        return IntPoint;
    }

    /*ThetaPath()
     *Calculates angle Theta and the path length
     */
    std::pair<double, double> Reconstructor::ThetaPath(double zInt, double zSi, double rSi)
    {
        double theta, path;
        theta = std::atan2(rSi, zInt - zSi);
        if (theta < 0.0)
            theta += M_PI;
        path = std::hypot(zInt - zSi, rSi);
        return std::make_pair(theta, path);
    }

    ROOT::Math::PxPyPzEVector Reconstructor::GetLorentzVectorFromTrack(double intPoint, const Track& event, ReconParticle type)
    {
        ROOT::Math::PxPyPzEVector vec;
        std::pair<double, double> theta_path;

        theta_path = ThetaPath(intPoint, event.siHit.z, event.siHit.r);

        double KE_rxn, E_rxn, P_rxn;
        switch (type)
        {
        case ReconParticle::Ejectile:
        {
            KE_rxn = event.siHit.energy + m_target.GetReverseEnergyLoss(m_ejectileData.Z, m_ejectileData.A, event.siHit.energy, theta_path.second / 100.0);
            E_rxn = KE_rxn + m_ejectileData.mass;
            P_rxn = std::sqrt(KE_rxn * (KE_rxn + 2.0 * m_ejectileData.mass));
            break;
        }
        case ReconParticle::Breakup1:
        {
            KE_rxn = event.siHit.energy + m_target.GetReverseEnergyLoss(m_breakup1Data.Z, m_breakup1Data.A, event.siHit.energy, theta_path.second / 100.0);
            E_rxn = KE_rxn + m_breakup1Data.mass;
            P_rxn = std::sqrt(KE_rxn * (KE_rxn + 2.0 * m_breakup1Data.mass));
            break;
        }
        case ReconParticle::Breakup2:
        {
            KE_rxn = event.siHit.energy + m_target.GetReverseEnergyLoss(m_breakup2Data.Z, m_breakup2Data.A, event.siHit.energy, theta_path.second / 100.0);
            E_rxn = KE_rxn + m_breakup2Data.mass;
            P_rxn = std::sqrt(KE_rxn * (KE_rxn + 2.0 * m_breakup2Data.mass));
            break;
        }
        }
        vec.SetPxPyPzE(P_rxn * std::sin(theta_path.first) * std::cos(event.siHit.phi),
                       P_rxn * std::sin(theta_path.first) * std::sin(event.siHit.phi),
                       P_rxn * std::cos(theta_path.first),
                       E_rxn);
        return vec;
    }
    /*CalcRecoil()
     *Reconstructs a recoil nucleus using the tracking information of a single ejectile, and assumed
     *knowledge of the beam and target. Should be compared to CaclRecoil_from_Qvalue, which does not
     *rely as much on the tracking to calculate the same value.
     */
    ReconResult Reconstructor::CalcResidualTracked(const Track& ejectile)
    {
        if (!IsValid())
        {
            fmt::print("At Reconstructor::CalcRecoil() Reconstruct not initialized! Skipping.\n");
            return ReconResult();
        }

        ReconResult result;

        ROOT::Math::PxPyPzEVector ejectVec, projVec, targVec;
        ejectVec = GetLorentzVectorFromTrack(ejectile.interactionPointZ, ejectile, ReconParticle::Ejectile);

        double beamE = ejectile.beamEnergy + m_projectileData.mass;
        double beamPz = std::sqrt(beamE * beamE - m_projectileData.mass * m_projectileData.mass);

        projVec.SetPxPyPzE(0.0, 0.0, beamPz, beamE);
        targVec.SetPxPyPzE(0.0, 0.0, 0.0, m_targetData.mass);

        ROOT::Math::PxPyPzEVector parentVec = projVec + targVec;
        ROOT::Math::PxPyPzEVector residualVec = parentVec - ejectVec;

        result.avgInteractionZ = ejectile.interactionPointZ;
        result.projectileVec = projVec;
        result.residualVec = residualVec;
        result.projectileRxnKE = projVec.E() - projVec.M();
        result.ejectileRxnKE = ejectVec.E() - ejectVec.M();
        result.parentKE = parentVec.E() - parentVec.M();
        result.parentExcitation = parentVec.M() - m_targetData.mass - m_projectileData.mass;
        result.residualExcitation = residualVec.M() - m_residualData.mass;

        ROOT::Math::Boost boostVec(parentVec.BoostToCM());
        ejectVec = boostVec * ejectVec;
        parentVec = boostVec * parentVec;
        result.ejectileThetaCM = ejectVec.Theta();
        result.parentECM = parentVec.E() - m_projectileData.mass - m_targetData.mass;

        return result;
    }

    /*CalcRecoil_from_Qvalue()
     *Reconstructs a recoil nucleus from a single ejectile assuming knowledge of the beam,
     *reaction Qvalue, and target. This is mainly meant to be a comparison to CalcRecoil() to
     *ensure that the tracking isn't super screwed up.
     */
    ReconResult Reconstructor::CalcResidualFromQvalue(const Track& ejectile)
    {
        if (!IsValid())
        {
            fmt::print("At Reconstructor::CalcRecoilFromQValues() Reconstruct not initialized! Skipping.\n");
            return ReconResult();
        }

        ReconResult result;

        auto ejectVec = GetLorentzVectorFromTrack(ejectile.interactionPointZ, ejectile, ReconParticle::Ejectile);
        double ejectKE = ejectVec.E() - m_ejectileData.mass;

        double Q = m_targetData.mass + m_projectileData.mass - m_ejectileData.mass - m_residualData.mass;
        double a = (m_projectileData.mass - m_residualData.mass);
        double b = (m_residualData.mass + m_ejectileData.mass);
        double c = -Q * m_residualData.mass;
        double d = 2.0 * std::sqrt(m_projectileData.mass * m_ejectileData.mass);
        double cos2 = std::cos(ejectile.theta) * std::cos(ejectile.theta);
        double a1 = a * a;
        double a2 = (2.0 * a * (b * ejectKE + c) - d * d * ejectKE * cos2);
        double a3 = (b * ejectKE + c) * (b * ejectKE + c);
        double s1 = (-a2 + std::sqrt(a2 * a2 - 4.0 * a1 * a3)) / (2.0 * a1);
        double s2 = (-a2 - std::sqrt(a2 * a2 - 4.0 * a1 * a3)) / (2.0 * a1);

        if (s1 < s2)
            std::swap(s1, s2);
        double beamKE = s2;
        if (std::isnan(s2))
            return result;

        double beamE = beamKE + m_projectileData.mass;
        double beamPz = sqrt(beamKE * (beamKE + 2.0 * m_projectileData.mass));
        ROOT::Math::PxPyPzEVector projVec(0., 0., beamPz, beamE);
        ROOT::Math::PxPyPzEVector targVec(0., 0., 0., m_targetData.mass);
        auto parentVec = projVec + targVec;
        auto residualVec = parentVec - ejectVec;

        result.avgInteractionZ = ejectile.interactionPointZ;
        result.projectileVec = projVec;
        result.residualVec = residualVec;
        result.projectileRxnKE = projVec.E() - projVec.M();
        result.ejectileRxnKE = ejectVec.E() - ejectVec.M();
        result.parentKE = parentVec.E() - parentVec.M();
        result.parentExcitation = parentVec.M() - m_targetData.mass - m_projectileData.mass;
        result.residualExcitation = residualVec.M() - m_residualData.mass;

        ROOT::Math::Boost boostVec(parentVec.BoostToCM());
        ejectVec = boostVec * ejectVec;
        parentVec = boostVec * parentVec;
        result.ejectileThetaCM = ejectVec.Theta();
        result.parentECM = parentVec.E() - m_projectileData.mass - m_targetData.mass;

        return result;
    }

    ReconResult Reconstructor::CalcMultiParticleTracked(const Track& ejectile, const Track& breakup1, const Track& breakup2)
    {
        static double massAlpha = MassLookup::GetInstance().FindMass(2, 4);
        static double p0 = -2.17304;
        static double p1 = 0.299308;
        static double p2 = -0.0104547; 
        ReconResult result;
        if (!IsValid())
        {
            fmt::print("At Reconstructor::CalculateMultiParticle() Reconstruct not initialized! Skipping.\n");
            return result;
        }

        ROOT::Math::PxPyPzEVector targVec(0., 0., 0., m_targetData.mass);
        
        double avgIntPoint, avgBeamKE;
        if (m_ejectileData.mass == massAlpha) // 5Li
        {
            avgIntPoint = (ejectile.interactionPointZ + breakup2.interactionPointZ) / 2.0;
            avgBeamKE = (ejectile.beamEnergy + breakup2.beamEnergy) / 2.0;
        }
        else
        {
            avgIntPoint = (breakup1.interactionPointZ + breakup2.interactionPointZ) / 2.0;
            avgBeamKE = (breakup1.beamEnergy + breakup2.beamEnergy) / 2.0;
        }
        avgIntPoint = avgIntPoint - (p0 + (p1 + p2 * avgIntPoint) * avgIntPoint);

        auto ejectVec = GetLorentzVectorFromTrack(avgIntPoint, ejectile, ReconParticle::Ejectile);
        auto break1Vec = GetLorentzVectorFromTrack(avgIntPoint, breakup1, ReconParticle::Breakup1);
        auto break2Vec = GetLorentzVectorFromTrack(avgIntPoint, breakup2, ReconParticle::Breakup2);

        ROOT::Math::PxPyPzEVector sumProductsVec = ejectVec + break1Vec + break2Vec;
        double beamE = avgBeamKE + m_projectileData.mass;
        double beamPz = std::sqrt(avgBeamKE * (avgBeamKE + 2.0 * m_projectileData.mass));
        ROOT::Math::PxPyPzEVector projVec(0., 0., beamPz, beamE);
        auto parentVec = projVec + targVec;
        auto residualVec = parentVec - ejectVec;

        result.avgInteractionZ = avgIntPoint;
        result.parentKE = parentVec.E() - parentVec.M();
        result.projectileRxnKE = projVec.E() - projVec.M();
        result.ejectileRxnKE = ejectVec.E() - ejectVec.M();
        result.projectileVec = projVec;
        result.residualVec = residualVec;
        result.parentExcitation = parentVec.M() - m_targetData.mass - m_projectileData.mass;
        result.residualExcitation = residualVec.M() - m_residualData.mass;
        result.residualMassSq = residualVec.M() * residualVec.M() * s_MeV2GeV * s_MeV2GeV;

        ROOT::Math::Boost boostVec(parentVec.BoostToCM());
        parentVec = boostVec * parentVec;
        residualVec = boostVec * residualVec;
        ejectVec = boostVec * ejectVec;
        result.ejectileThetaCM = ejectVec.Theta();
        result.parentECM = parentVec.E() - m_projectileData.mass - m_targetData.mass;

        return result;
    }

    ReconResult Reconstructor::CalcMultiParticleQvalue(const Track& ejectile, const Track& breakup1, const Track& breakup2)
    {
        static double massAlpha = MassLookup::GetInstance().FindMass(2, 4);
        static double p0 = -2.17304;
        static double p1 = 0.299308;
        static double p2 = -0.0104547; 
        ReconResult result;
        if (!IsValid())
        {
            fmt::print("At Reconstructor::CalculateMultiParticle() Reconstruct not initialized! Skipping.\n");
            return result;
        }

        ROOT::Math::PxPyPzEVector targVec(0., 0., 0., m_targetData.mass);
        
        double avgIntPoint, avgBeamKE;
        if (m_ejectileData.mass == massAlpha) // 5Li
        {
            avgIntPoint = (ejectile.interactionPointZ + breakup2.interactionPointZ) / 2.0;
            avgBeamKE = (ejectile.beamEnergy + breakup2.beamEnergy) / 2.0;
        }
        else
        {
            avgIntPoint = (breakup1.interactionPointZ + breakup2.interactionPointZ) / 2.0;
            avgBeamKE = (breakup1.beamEnergy + breakup2.beamEnergy) / 2.0;
        }
        avgIntPoint = avgIntPoint - (p0 + (p1 + p2 * avgIntPoint) * avgIntPoint);

        auto ejectVec = GetLorentzVectorFromTrack(avgIntPoint, ejectile, ReconParticle::Ejectile);
        auto break1Vec = GetLorentzVectorFromTrack(avgIntPoint, breakup1, ReconParticle::Breakup1);
        auto break2Vec = GetLorentzVectorFromTrack(avgIntPoint, breakup2, ReconParticle::Breakup2);

        ROOT::Math::PxPyPzEVector sumProductsVec = ejectVec + break1Vec + break2Vec;
        auto projVec = sumProductsVec - targVec;
        auto residualVec = break1Vec + break2Vec;

        result.avgInteractionZ = avgIntPoint;
        result.parentKE = sumProductsVec.E() - sumProductsVec.M();
        result.projectileRxnKE = projVec.E() - projVec.M();
        result.ejectileRxnKE = ejectVec.E() - ejectVec.M();
        result.projectileVec = projVec;
        result.residualVec = residualVec;
        result.parentExcitation = sumProductsVec.M() - m_targetData.mass - m_projectileData.mass;
        result.residualExcitation = residualVec.M() - m_residualData.mass;
        result.residualMassSq = residualVec.M() * residualVec.M() * s_MeV2GeV * s_MeV2GeV;

        ROOT::Math::Boost boostVec(sumProductsVec.BoostToCM());
        sumProductsVec = boostVec * sumProductsVec;
        residualVec = boostVec * residualVec;
        ejectVec = boostVec * ejectVec;
        result.ejectileThetaCM = ejectVec.Theta();
        result.parentECM = sumProductsVec.E() - m_projectileData.mass - m_targetData.mass;

        return result;
    }

    double Reconstructor::BeamEnergyNoTracking(const Track &p1, const Track &p2, const Track &p3)
    {
        double product_E = p1.siHit.energy + p2.siHit.energy + p3.siHit.energy + m_ejectileData.mass + m_breakup1Data.mass + m_breakup2Data.mass;
        double proj_E = product_E - m_targetData.mass;
        return proj_E - m_projectileData.mass;
    }

    // Reconstruct recoil from daughter nuclei
    ReconResult Reconstructor::CalcTwoParticle(const Track& breakup1, const Track& breakup2)
    {
        static double massAlpha = MassLookup::GetInstance().FindMass(2, 4);
        ReconResult result;
        if (!IsValid())
        {
            fmt::print("At Reconstructor::CalcTwoParticle() Reconstruct not initialized! Skipping.\n");
            return result;
        }
        ROOT::Math::PxPyPzEVector targetVec(0.0, 0.0, 0.0, m_targetData.mass);

        // Avg alpha int points, or just take alpha in case of 5Li
        double avgIntPoint;
        double avgBeamKE;
        if (m_ejectileData.mass == massAlpha) // 5Li
        {
            avgIntPoint = breakup2.interactionPointZ;
            avgBeamKE = breakup2.beamEnergy;
        }
        else
        {
            avgIntPoint = (breakup1.interactionPointZ + breakup2.interactionPointZ) / 2.0;
            avgBeamKE = (breakup1.beamEnergy + breakup2.beamEnergy) / 2.0;
        }

        auto break1Vec = GetLorentzVectorFromTrack(avgIntPoint, breakup1, ReconParticle::Breakup1);
        auto break2Vec = GetLorentzVectorFromTrack(avgIntPoint, breakup2, ReconParticle::Breakup2);

        auto residualVec = break1Vec + break2Vec;
        double beamPz = std::sqrt(avgBeamKE * (avgBeamKE + 2.0 * m_projectileData.mass));
        ROOT::Math::PxPyPzEVector projVec(0.0, 0.0, beamPz, avgBeamKE + m_projectileData.mass);
        double energyCM = avgBeamKE * m_targetData.mass / (m_targetData.mass + m_projectileData.mass);


        result.avgInteractionZ = avgIntPoint;
        result.parentKE = 0.0;
        result.parentECM = energyCM;
        result.projectileVec = projVec;
        result.residualVec = residualVec;
        result.residualExcitation = residualVec.M() - m_residualData.mass;

        return result;
    }
}