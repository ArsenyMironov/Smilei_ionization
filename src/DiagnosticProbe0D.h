#ifndef DiagnosticProbe0D_H
#define DiagnosticProbe0D_H

#include "Tools.h"


#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <hdf5.h>
#include <math.h>
#include "Species.h"
#include "Interpolator.h"
#include "Particle.h"

class PicParams;
class SmileiMPI;
class DiagParams;
class ElectroMagn;


class DiagnosticProbe0D{
    
public:
    
	DiagnosticProbe0D(PicParams* params, SmileiMPI* smpi,std::vector<std::vector<double> > ps_coord);
	~DiagnosticProbe0D();
	
	void set_proc();
    void set_file_name();
    
    void run(int timestep, ElectroMagn* EMfields, Interpolator* interp);
	
	std::string probeName(int p);
    void open_file();
    void close();
 
private:
	SmileiMPI* smpi_;

    hid_t probeFileGlobal_id;
	MPI_Comm comm;
    MPI_Info info;
	std::vector<Particle*> probeParticles;
	std::vector<bool> probeId;
    std::vector<LocalFields> Eloc_fields;
    std::vector<LocalFields> Bloc_fields;
    std::vector<hid_t> hidGroup;
	int probeSize;
        
};
#endif 