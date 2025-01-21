#ifndef IONIZATIONTUNNELCUSTOMCOEFFBSI_H
#define IONIZATIONTUNNELCUSTOMCOEFFBSI_H

#include <cmath>
#include <vector>

#include "Ionization.h"
#include "Tools.h"

class Particles;


class IonizationTunnelCustomCoeffBSI : public Ionization {
    public:
        IonizationTunnelCustomCoeffBSI(Params &params, Species *species);
        void operator()(Particles*, unsigned int, unsigned int, std::vector<double>*, Patch*, Projector*, int ipart_ref = 0) override;

    private:
        unsigned int atomic_number_;
        std::vector<double> Potential, Azimuthal_quantum_number;
        double one_third;
        std::vector<double> alpha_tunnel, beta_tunnel, gamma_tunnel; 

        std::vector<double> Magnetic_quantum_number;
        unsigned int cnl_model_;
        bool m_equal_zero_, use_g_factor_;
        std::vector<double> cnl_squared_table_;
        std::vector<double> g_factor;

};

#endif