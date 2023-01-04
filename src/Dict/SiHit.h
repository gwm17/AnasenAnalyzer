/*
    Legacy Si data class. CANNOT BE MODIFIED when using old data. For new data, a new Si data structure should be developed. 

    Some attempt to make stuff more readable, but the dictionary really constrains the ability to fix this.

    Also here, event/hit terminology is swapped because legacy. (Hit should refer to single detector datum, event all data in window)

    GWM -- Dec 2022
*/
#ifndef SIHIT_H
#define SIHIT_H

#include <vector>

// The fields in both structs cannot be modified when using with old data
struct SiDetector
{
    int DetID, UpMult, DownMult, FrontMult, BackMult, HitType;
    std::vector<int> UpChNum, DownChNum, FrontChNum, BackChNum;

    // Raw Energy (units of channels)
    std::vector<double> EUp_Raw, EDown_Raw, EFront_Raw, EBack_Raw;

    // Energy after pulser cals (units of channels)
    std::vector<double> EUp_Pulser, EDown_Pulser, EFront_Pulser, EBack_Pulser;

    // Energy after relative calibration (units of channels)
    std::vector<double> EUp_Rel, EDown_Rel, EFront_Rel, EBack_Rel, SX3_ZUp, SX3_ZDown,
        EUp_RelShift, EDown_RelShift, EBack_RelShift;
    // Energy after energy cal (units of MeV)
    std::vector<double> EUp_Cal, EDown_Cal, EFront_Cal, EBack_Cal, TUp, TDown, TFront, TBack;
};

struct SiEvent
{
    int NHitsInDet, DetID, HitType;
    int FrontChannel, BackChannel, UpChannel, DownChannel;
    double EnergyBack, EnergyFront, EnergyUp, EnergyDown, Energy,
        Time, X, Y, Z, Z_linear, ZUp_Dummy, ZDown_Dummy, XW, YW, ZW, RW, PhiW;
    int TrackType;
    double RFTime, MCPTime;
};

class SiHit
{

public:
    int NSiHits; // std::vector already has .size() property, makes no sense to duplicate

    // SortByDetector/Hit is the original name (which makes no sense at all), so we just alias the much more appropriate SiDetector/Event
    struct SortByDetector : SiDetector
    {
    };
    struct SortByHit : SiEvent
    {
    };
    /*
        There is so much wrong with this. Only need vectors, no idea why there are single instances of each event type. Definitely don't repeat pointers to
        vectors. But since these were in the orig. dictionary they cannot be removed.

        Also not clear why this data structure references two different types of data (detector and event based). Only event is used by AnasenAnalyzer.
    */
    SortByDetector det_obj;
    SortByHit hit_obj;
    std::vector<SortByDetector> Detector;
    std::vector<SortByHit> Hit;
    std::vector<SortByDetector> *ReadDet;
    std::vector<SortByHit> *ReadHit;

    // Love setting pointers to 0
    SiHit()
    {
        ReadDet = 0;
        ReadHit = 0;
    };

    // Ick
    void ZeroSi_obj()
    {
        det_obj.DetID = -1;
        det_obj.UpMult = 0;
        det_obj.DownMult = 0;
        det_obj.BackMult = 0;
        det_obj.FrontMult = 0;
        det_obj.HitType = 0;

        det_obj.UpChNum.clear();
        det_obj.DownChNum.clear();
        det_obj.FrontChNum.clear();
        det_obj.BackChNum.clear();

        det_obj.EUp_Raw.clear();
        det_obj.EDown_Raw.clear();
        det_obj.EFront_Raw.clear();
        det_obj.EBack_Raw.clear();

        det_obj.EUp_Pulser.clear();
        det_obj.EDown_Pulser.clear();
        det_obj.EFront_Pulser.clear();
        det_obj.EBack_Pulser.clear();

        det_obj.EUp_Rel.clear();
        det_obj.EDown_Rel.clear();
        det_obj.EFront_Rel.clear();
        det_obj.EBack_Rel.clear();

        det_obj.EUp_RelShift.clear();
        det_obj.EDown_RelShift.clear();
        det_obj.EBack_RelShift.clear();

        det_obj.SX3_ZUp.clear();
        det_obj.SX3_ZDown.clear();

        det_obj.EUp_Cal.clear();
        det_obj.EDown_Cal.clear();
        det_obj.EFront_Cal.clear();
        det_obj.EBack_Cal.clear();

        det_obj.TUp.clear();
        det_obj.TDown.clear();
        det_obj.TFront.clear();
        det_obj.TBack.clear();

        hit_obj.NHitsInDet = 0;
        hit_obj.DetID = -1;

        hit_obj.HitType = 0;
        hit_obj.FrontChannel = -1;
        hit_obj.BackChannel = -1;
        hit_obj.UpChannel = -1;
        hit_obj.DownChannel = -1;
        hit_obj.EnergyBack = -1000;
        hit_obj.EnergyFront = -1000;
        hit_obj.EnergyUp = -1000;
        hit_obj.EnergyDown = -1000;
        hit_obj.Energy = -1000;
        hit_obj.Time = 0;
        hit_obj.X = 0;
        hit_obj.Y = 0;
        hit_obj.Z = -10;
        hit_obj.Z_linear = -10; // added 12/09/2016 for the ZPosCal_linear
        hit_obj.ZUp_Dummy = -10;
        hit_obj.ZDown_Dummy = -10;
        hit_obj.XW = 0;
        hit_obj.YW = 0;
        hit_obj.ZW = -10;
        hit_obj.RW = 0;
        hit_obj.PhiW = 0;
        hit_obj.TrackType = 0;
        // hit_obj.RFSubtract = 0;
    };

    // Gross
    void ZeroSiHit()
    {
        NSiHits = 0;
        Detector.clear();
        Hit.clear();
        ZeroSi_obj();
    };
};

#endif