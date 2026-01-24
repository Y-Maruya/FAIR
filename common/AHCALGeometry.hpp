#pragma once
#include <vector>
// #include "EventFormats/AHCALDataFragment.hpp"

namespace AHCALGeometry {
    const int channel_FEE = 73;//(36charges+36times + BCIDs )*16column+ ChipID
    const int cell_SP = 16;
    const int chip_No = 9;
    const int channel_No = 36;
    const int Layer_No = 40;
    const double _Pos_X[channel_No]={100.2411,100.2411,100.2411,59.94146,59.94146,59.94146,19.64182,19.64182,19.64182,19.64182,59.94146,100.2411,100.2411,59.94146,19.64182,100.2411,59.94146,19.64182,-20.65782,-60.95746,-101.2571,-20.65782,-60.95746,-101.2571,-101.2571,-60.95746,-20.65782,-20.65782,-20.65782,-20.65782,-60.95746,-60.95746,-60.95746,-101.2571,-101.2571,-101.2571};
    const double _Pos_Y[channel_No]={141.04874,181.34838,221.64802,141.04874,181.34838,221.64802,141.04874,181.34838,221.64802,261.94766,261.94766,261.94766,302.2473,302.2473,302.2473,342.54694,342.54694,342.54694,342.54694,342.54694,342.54694,302.2473,302.2473,302.2473,261.94766,261.94766,261.94766,221.64802,181.34838,141.04874,221.64802,181.34838,141.04874,221.64802,181.34838,141.04874};
    const int _Channel[6][6]={  { 0, 1, 2,11,12,15},
                                { 3, 4, 5,10,13,16},
                                { 6, 7, 8, 9,14,17},
                                {29,28,27,26,21,18},
                                {32,31,30,25,22,19},
                                {35,34,33,24,23,20}};
    const int _Chip[3][3]={ {2,1,0},
                            {5,4,3},
                            {8,7,6}};
                            
    const double chip_dis_X=239.3;
    const double chip_dis_Y=241.8;
    const double HBU_X=239.3;
    const double HBU_Y=725.4;
    const double x_max = 40.3*18/2;
    const double y_max = 40.3*18/2;
    const double xy_size = 40.;
    const double z_size = 3.;
    inline double Pos_X(int channel_ID,int chip_ID){
        // int HBU_ID=chip_ID/3;
        chip_ID=chip_ID%3;
        if(chip_ID!=0){
            if(channel_ID==2)channel_ID=0;
            else if(channel_ID==0)channel_ID=2;
            if(channel_ID==33)channel_ID=35;
            else if(channel_ID==35)channel_ID=33;
        }
        return (_Pos_Y[channel_ID]-chip_ID*chip_dis_Y);
    }
    inline double Pos_Y(int channel_ID,int chip_ID,int HBU_ID=0){
        HBU_ID=chip_ID/3;
        return -(-_Pos_X[channel_ID]+(HBU_ID-1)*HBU_X);
    }
    inline double Pos_Z(int layer_ID){
        return layer_ID*29.63 + 1.5; // start from first layer scintillator front face // TEMPORARY 
        // TODO : read the geometry file
    }
    inline int HBUPositionOrder[40] = {
                                39, 38, 37, 27, 14, 6, 7, 9, 12, 0,
                        2, 3, 5, 8, 10, 11, 13, 15, 16, 1, 
                        17, 18, 19, 20, 21, 22, 23, 24, 25, 4, 
                        26, 28, 29, 30, 31, 32, 33, 35, 34, 36
                        };
    inline int PosToLayerID(int position_order){
        for (int i = 0; i < 40; ++i){
            if (HBUPositionOrder[i] == position_order){
                return i;
            }
        }
        return -1; // not found
    }
}

