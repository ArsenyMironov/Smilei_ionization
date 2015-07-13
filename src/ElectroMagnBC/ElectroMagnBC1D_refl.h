
#ifndef ELECTROMAGNBC1D_REFL_H
#define ELECTROMAGNBC1D_REFL_H

#include "ElectroMagnBC.h" 

class PicParams;
class ElectroMagn;

class ElectroMagnBC1D_refl : public ElectroMagnBC {
public:
    ElectroMagnBC1D_refl( PicParams &param, LaserParams &laser_params);
    ~ElectroMagnBC1D_refl();

    virtual void apply_xmin(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi);
    virtual void apply_xmax(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi);
    virtual void apply_ymin(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi);
    virtual void apply_ymax(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi);
    
 private:
    
    //! Oversize
    unsigned int oversize_;
    
    //! Number of nodes on the primal grid
    unsigned int nx_p;

    //! Number of nodes on the dual grid
    unsigned int nx_d;
    
};

#endif

