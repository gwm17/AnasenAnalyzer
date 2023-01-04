/*
    Legacy PC data class. CANNOT BE MODIFIED when using old data. For new data, a new PC data structure should be developed. 

    Some attempt to make stuff more readable, but the dictionary really constrains the ability to fix this.

    Also here, event/hit terminology is swapped because legacy. (Hit should refer to single detector datum, event all data in window)

    GWM -- Dec 2022
*/
#ifndef PCHIT_H
#define PCHIT_H

#include <vector>

// These fields cannot be modified when using old data
struct PCEvent
{
    int WireID;
    double Down, Up, DownVoltage, UpVoltage, Energy,
        Z, XW, YW, ZW, RW, PhiW;
    int TrackType;
};

class PCHit
{

public:
    int NPCHits; // std::vector already has .size() property, makes no sense to duplicate
    // SortByPC is the original name (which makes no sense at all), so we just alias the much more appropriate PCEvent
    struct SortByPC : PCEvent
    {
    };
    /*
        There is so much wrong with this. Only need vectors, no idea why there are single instance of event type. Definitely don't repeat pointers to
        vectors. But since these were in the orig. dictionary they cannot be removed.
    */
    SortByPC pc_obj;
    std::vector<PCHit::SortByPC> Hit;
    std::vector<SortByPC> *ReadHit;

    PCHit(){};

    // Gross
    void ZeroPC_obj()
    {
        pc_obj.WireID = -1;
        pc_obj.Down = -10.0;
        pc_obj.Up = -10.0;
        pc_obj.DownVoltage = -10.0;
        pc_obj.UpVoltage = -10.0;
        pc_obj.Energy = -10.0;
        pc_obj.Z = -100.0;
        pc_obj.XW = 0.0;
        pc_obj.YW = 0.0;
        pc_obj.ZW = -10.0;
        pc_obj.RW = 0.0;
        pc_obj.PhiW = 0.0;
        pc_obj.TrackType = 0;
    };

    // Amazing
    void ZeroPCHit()
    {
        NPCHits = 0;
        Hit.clear();
        ZeroPC_obj();
    };
};

#endif