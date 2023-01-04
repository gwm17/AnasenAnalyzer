/*
    Track.h
    Defines complete event data structure. This is not legacy, however, in an ideal world, the data would already be represented as such
    by the time it reaches analysis. The fact that the analyzer both event builds and analyzes is not ideal. More appropriate as final stage of
    calibration.

    In current iteration we never write Track to disk, so can remove from dictionary

    GWM -- Dec 2022
*/
#ifndef TRACK_H
#define TRACK_H

#include <vector>
#include "SiHit.h"
#include "PCHit.h"
#include "Math/Vector4D.h"

namespace AnasenAnalyzer {

    // Where possible use enums to define state
    enum class TrackType
    {
        Complete,
        SiOnly,
        PCOnly,
        None
    };

    // New struct representing silicon/pc info for track
    // Note initializing to default invalid values. Saves a lot of headaches later.
    struct TrackedSiliconHit
    {
        int detID = -1;
        double energy = -1.0;
        double time = -1.0;
        double z = -1.0;
        double r =  -1.0;
        double phi = -1.0;
    };

    struct TrackedPCHit
    {
        int wireID = -1;
        double energy = -1.0;
        double z = -1.0;
        double r = -1.0;
        double phi = -1.0;
    };

    // Conversion functions for old data type to new data type
    static TrackedSiliconHit CreateTrackedSiHit(const SiEvent &event)
    {
        return {event.DetID, event.Energy, event.Time, event.ZW, event.RW, event.PhiW};
    }

    static TrackedPCHit CreateTrackedPCHit(const PCEvent& event)
    {
        return {event.WireID, event.Energy, event.ZW, event.RW, event.PhiW};
    }

    struct Track
    {
        TrackType type = TrackType::None;
        
        TrackedSiliconHit siHit;
        TrackedPCHit pcHit;

        // Interaction point coords
        double interactionPointZ = -10;
        // Track info
        double pathLength = -10, theta = -10, phi = -10;
        // Energetics
        double beamEnergy = -10, beamEnergyLoss = -10;

        bool IsValid() const { return siHit.detID > -1 && siHit.detID < 28 && pcHit.z > 0.0 && siHit.z >= 0.0 && interactionPointZ > 0.0 && interactionPointZ < 55.0545; } 
        bool IsBeamValid(double max) const { return beamEnergy > 0.0 && beamEnergy < max; }
    };

    // Not sure these structs are totally necessary; basically simplified versions of each
    // class' sort by hit and already in TrackEvent... Turns out it is necessary for backwards
    // compatibility with dictionaries
    struct Silicon_Event
    {
        int TrackType, DetID;
        double SiEnergy, SiTime, SiZ, SiR, SiPhi;
    };

    struct PropCounter_Event
    {
        int TrackType, WireID;
        double PCEnergy, PCZ, PCR, PCPhi, Down, Up, DownVoltage, UpVoltage;
    };

    struct CsI_Event
    {
        int TrackType, CsI_ID;
        double CsIEnergy, CsI_Phi, CsI_X, CsI_Y, CsI_R;
    };


    

    struct TrackEvent
    {
        std::vector<Track> complete;
        std::vector<Track> siOnly;
        std::vector<Track> pcOnly;

        void ClearAll()
        {
            complete.clear();
            siOnly.clear();
            pcOnly.clear();
        }
        const std::size_t GetTotalTracks() const { return complete.size() + siOnly.size() + pcOnly.size(); }
        static bool SortTracksBySi(Track& a, Track& b) { return a.siHit.energy > b.siHit.energy ? true : false; }
        static bool SortTracksByPC(Track& a, Track& b) { return a.pcHit.energy > b.pcHit.energy ? true : false; }
    };

    bool EnforceDictionaryLinked();
}

#endif