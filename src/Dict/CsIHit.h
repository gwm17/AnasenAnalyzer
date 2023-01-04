#ifndef CSIHIT_H
#define CSIHIT_H

#include <vector>

struct CsIEvent
{
    int CsI_ID, mADC_ID, mADC_Ch;
    double Up, Down, QQQ, Up_Shift, Down_Shift, QQQ_Shift;
    double WCsI_X, WCsI_Y, WCsI_Z, WCsI_R, WCsI_Phi;
    double Energy;
};

class CsIHit
{

public:
    std::vector<double> CsI_Energy;
    int NCsIHits;
    struct SortByCsI : CsIEvent
    {
    };
    SortByCsI csi_obj;
    std::vector<SortByCsI> Hit;
    std::vector<SortByCsI> *ReadHit;

    CsIHit(){};

    void ZeroCsI_obj()
    {
        csi_obj.CsI_ID = -1;
        csi_obj.mADC_Ch = -1;
        csi_obj.mADC_ID = -1;

        csi_obj.Up = -10;
        csi_obj.Down = -10;
        csi_obj.QQQ = -10;

        csi_obj.Up_Shift = -10;
        csi_obj.Down_Shift = -10;
        csi_obj.QQQ_Shift = -10;

        csi_obj.Energy = -10;
        csi_obj.WCsI_X = -1000;
        csi_obj.WCsI_Y = -1000;
        csi_obj.WCsI_R = -1000;
        csi_obj.WCsI_Z = -1000;
        csi_obj.WCsI_Phi = -1000;
    };

    void ZeroCsIHit()
    {
        NCsIHits = 0;
        Hit.clear();
        ZeroCsI_obj();
    };
};

#endif