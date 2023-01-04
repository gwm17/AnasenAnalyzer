#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#include <Math/Vector4D.h>
#include "Target.h"
#include "Track.h"
#include "CsIHit.h"
#include "SiHit.h"
#include "PCHit.h"
#include <cmath>


namespace AnasenAnalyzer {

    /*
            target(projectile, ejectile)residual->breakup1+breakup2
            For Z and A goes as 0=target, 1=projectile, 2=ejectile, 3=breakup1
            Breakup1 and 2 exist only as needed.
    */

    // Result of Reconstructor
    struct ReconResult
    {
        double avgInteractionZ = -1.0;
        double parentKE = -1.0;
        double parentECM = -1.0;
        double projectileRxnKE = -1.0;
        double ejectileRxnKE = -1.0;
        double ejectileThetaCM = -1.0;
        double parentExcitation = -1.0;
        double residualExcitation = -1.0;
        double residualMassSq = -1.0;

        ROOT::Math::PxPyPzEVector residualVec;
        ROOT::Math::PxPyPzEVector projectileVec;
    };

    enum class ReconParticle
    {
        Ejectile,
        Breakup1,
        Breakup2
    };

    class Reconstructor
    {
    public:
        Reconstructor(const Target& target, const std::vector<uint32_t>& Z, const std::vector<uint32_t>& A);
        ~Reconstructor();


        const bool IsValid() const { return m_isInit; }
        ReconResult CalcResidualTracked(const Track& ejectile);
        ReconResult CalcResidualFromQvalue(const Track& ejectile);
        ReconResult CalcMultiParticleTracked(const Track& ejectile, const Track& breakup1, const Track& breakup2);
        ReconResult CalcMultiParticleQvalue(const Track& ejectile, const Track& breakup1, const Track& breakup2);
        double BeamEnergyNoTracking(const Track& p1, const Track& p2, const Track& p3);
        ReconResult CalcTwoParticle(const Track& breakup1, const Track& breakup2);

    private:
        void Init(const std::vector<uint32_t>& Z, const std::vector<uint32_t>& A);
        double DeltaZ(double Zpc, double Zsi, double Rpc, double Rsi);
        double WeightedAverageIntP(const std::vector<Track> &tracks);
        std::pair<double, double> ThetaPath(double IntPoint, double Zsi, double Rsi);
        ROOT::Math::PxPyPzEVector GetLorentzVectorFromTrack(double intPoint, const Track& event, ReconParticle type);

        bool m_isInit;
        Target m_target;
        NucData m_projectileData, m_targetData, m_ejectileData, m_residualData, m_breakup1Data, m_breakup2Data;

        static constexpr double s_rad2deg = 180.0/M_PI;
        static constexpr double s_MeV2GeV = 1.0e-3;
    };

}

#endif