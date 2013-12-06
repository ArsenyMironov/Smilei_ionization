#ifndef Probe0D_H
#define Probe0D_H

#include "Tools.h"


#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <hdf5.h>
#include <math.h>
#include "Species.h"
#include "Interpolator.h"

class PicParams;
class SmileiMPI;
class DiagParams;
class ElectroMagn;


class Probe0D{
    
public:
    
	Probe0D(PicParams* params, SmileiMPI* smpi,std::vector<std::vector<double> > ps_coord);
	
	inline void set_p_coor(unsigned int i, std::vector<double> values){ps_coor[i]=values;}//write it in a smarter way
	inline std::vector<std::vector<double> > get_ps_coor(){return ps_coor;}
	inline std::vector<double> get_ps_coor(unsigned int i){return ps_coor[i];}
	
	void set_proc();
	
	void set_weights(); 
    
    void run(int timestep, ElectroMagn* EMfields, Interpolator* interp);
    
    
private:
    PicParams* params_;
    SmileiMPI* smpi_;
    unsigned int n_prob;
    std::vector<std::vector<double> > ps_coor;
    std::vector<std::vector<double> > weights_primal;
    std::vector<std::vector<double> > weights_dual;
    std::vector<bool> here;
};
#endif 