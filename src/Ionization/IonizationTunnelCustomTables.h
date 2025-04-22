#ifndef IONIZATIONTUNNELCUSTOMTABLES_H
#define IONIZATIONTUNNELCUSTOMTABLES_H

#include <cmath>

#include <vector>

#include "Ionization.h"
#include "Tools.h"

class Particles;

//! calculate the particle tunnel ionization
class IonizationTunnelCustomTables : public Ionization
{

public:
    //! Constructor for IonizationTunnelCustomTables: with no input argument
    IonizationTunnelCustomTables( Params &params, Species *species );
    
    //! apply the Tunnel Ionization model to the species (with ionization current)
    void operator()( Particles *, unsigned int, unsigned int, std::vector<double> *, Patch *, Projector *, int ipart_ref = 0 ) override;
    //! method for tunnel ionization with tasks
    // void ionizationTunnelWithTasks( Particles *, unsigned int, unsigned int, std::vector<double> *, Patch *, Projector *, int, int, double *b_Jx, double *b_Jy, double *b_Jz, int ipart_ref = 0 ) override;
    
private:
    unsigned int atomic_number_;
    std::vector<double> Potential;
    std::vector<double> Azimuthal_quantum_number;
    std::vector<double> Magnetic_quantum_number;
    
    double one_third;
    std::vector<double> alpha_tunnel, beta_tunnel, gamma_tunnel;
    std::vector<double> g_factor;

    unsigned int cnl_model_;
    bool m_equal_zero_;
    std::vector<double> cnl_squared_table_;
};


#endif
