#ifndef IONIZATIONTUNNEL_H
#define IONIZATIONTUNNEL_H

#include <cmath>
#include <functional>
#include <vector>

#include "Ionization.h"
#include "IonizationTables.h"
#include "Particles.h"
#include "Species.h"
#include "Tools.h"

class Particles;

// int Tunneling_Model : choice of the tunneling model 
//                       0 - for the ADK (l* = n*-1) model with the magnetic quantum number m set to 0 for all electrons
//                       1 - for the PPT model, in which A_nl is given by the Hartree formula and m can be non-zero for p-, d-, ... states
//                       2 - PPT model with ionization pot-s, l, m, and g quiantum numbers taken from custom tables passed in Species
//
// int BSI             : the choice of the barrier suppression ionization model
//                       0 for no barrier suppression (the tunneling formula is used)
//                       1 - Tong-Lin exponential suppression formula
//                       2 - Kostyukov-Artemenko-Golovanov model
//                       All BSI models can be combined with different tunneling formulas switched by Tunneling_Model 

template <int Tunneling_Model, int BSI>
class IonizationTunnel : public Ionization
{
   public:
    inline IonizationTunnel(Params &params, Species *species);

    inline void operator()(Particles *, unsigned int, unsigned int, std::vector<double> *, Patch *, Projector *,
                           int ipart_ref = 0) override;

   private:
    inline double ionizationRate(const int Z, const double E);

    static constexpr double one_third = 1. / 3.;
    unsigned int atomic_number_;
    std::vector<double> Potential, Azimuthal_quantum_number;
    std::vector<double> alpha_tunnel, beta_tunnel, gamma_tunnel;

    // Tong&Lin
    std::vector<double> lambda_tunnel;
};


template <int Tunneling_Model, int BSI>
IonizationTunnel<Tunneling_Model, BSI>::IonizationTunnel(Params &params, Species *species) : Ionization(params, species)
{
    DEBUG("Creating the Tunnel Ionizaton class");
    double abs_m         = 0;
    double g_factor      = 1;
    double Anl           = 4;  // the initial value is set to 4 on purpose, 
                               // as A_nl = 4*C_nl^2, where C_nl are 
                               // the Hartree coefficients; C_nl is set 
                               // to 1 for a neutral atom
    double Blm           = 1;
    double ionization_tl_parameter_;
    std::vector<double> Magnetic_quantum_number;
    std::vector<double> g_factors;

    // Ionization potential & quantum numbers (all in atomic units 1 au = 27.2116 eV)
    atomic_number_ = species->atomic_number_;
    if (Tunneling_Model == 2) {
        Potential                   = species->ionization_potentials_;
        Azimuthal_quantum_number    = species->azimuthal_quantum_numbers_;
        Magnetic_quantum_number     = species->magnetic_quantum_numbers_;
        g_factors                   = species->g_factors_;
    } else {
        Potential.resize(atomic_number_);
        Azimuthal_quantum_number.resize(atomic_number_);
    }

    alpha_tunnel.resize(atomic_number_);
    beta_tunnel.resize(atomic_number_);
    gamma_tunnel.resize(atomic_number_);

    if (BSI == 1) {
        ionization_tl_parameter_ = species->ionization_tl_parameter_;  // species->ionization_tl_parameter_ is
                                                                       // double Varies from 6 to 9. This is
                                                                       // the alpha parameter in Tong-Lin
                                                                       // exponential, see Eq. (6) in [M F
                                                                       // Ciappina and S V Popruzhenko 2020
                                                                       // Laser Phys. Lett. 17 025301 2020].
        lambda_tunnel.resize(atomic_number_);
    }

    for (unsigned int Z = 0; Z < atomic_number_; Z++) {
        DEBUG("Z : " << Z);

        if (Tunneling_Model == 1) {
            abs_m = abs(IonizationTables::magnetic_atomic_number(atomic_number_, Z));
            g_factor = IonizationTables::magnetic_degeneracy_atomic_number(atomic_number_, Z);
        } 

        if (Tunneling_Model == 2) {
            abs_m    = abs(Magnetic_quantum_number[Z]);
            g_factor = g_factors[Z];
            Potential[Z] = Potential[Z] * eV_to_au;
        } else {
            Potential[Z] = IonizationTables::ionization_energy(atomic_number_, Z) * eV_to_au;
            Azimuthal_quantum_number[Z] = IonizationTables::azimuthal_atomic_number(atomic_number_, Z);
        }

        DEBUG("Potential: " << Potential[Z] << " Az.q.num: " << Azimuthal_quantum_number[Z]);

        Blm      = ( 2.*Azimuthal_quantum_number[Z]+1.0 ) * \
                   tgamma(Azimuthal_quantum_number[Z]+abs_m+1) / \
                   ( pow( 2, abs_m )*tgamma(abs_m+1)*tgamma(Azimuthal_quantum_number[Z]-abs_m+1) );

        double cst = ((double)Z + 1.0) * sqrt(2.0 / Potential[Z]);
        if (Tunneling_Model == 0)  {
            Anl = pow( 2, cst+1.0 ) / \
                ( cst*tgamma( cst ) );
            
        } else {
            if( Z>0 ) {
                Anl = pow( 2, cst+1.0 ) / \
                                ( cst*tgamma( cst/2.0+Azimuthal_quantum_number[Z]+1 )*tgamma( cst/2.0-Azimuthal_quantum_number[Z]) );
            }
        }

        // MESSAGE(0, "Z=" << Z << ", Ip=" << Potential[Z] << ", l=" << Azimuthal_quantum_number[Z] << ", |m|=" << abs_m << ", g=" << g_factor << ", Anl=" << Anl << ", Blm=" << Blm);


        alpha_tunnel[Z] = cst - 1.0 - abs_m;
        beta_tunnel[Z] = g_factor*Anl*Blm * Potential[Z] * au_to_w0;
        gamma_tunnel[Z] = 2.0 * sqrt(2.0 * Potential[Z] * 2.0 * Potential[Z] * 2.0 * Potential[Z]);
        if (BSI == 1) {
            lambda_tunnel[Z] = ionization_tl_parameter_ * cst * cst / gamma_tunnel[Z];
        }
    }

    DEBUG("Finished Creating the Tunnel Ionizaton class");
}

template <int Tunneling_Model, int BSI>
inline void IonizationTunnel<Tunneling_Model, BSI>::operator()(Particles *particles, unsigned int ipart_min, unsigned int ipart_max,
                                                vector<double> *Epart, Patch *patch, Projector *Proj, int ipart_ref)
{
    unsigned int Z, Zp1, newZ, k_times;
    double TotalIonizPot, E, invE, factorJion, ran_p, Mult, D_sum, P_sum, Pint_tunnel;
    vector<double> IonizRate_tunnel(atomic_number_), Dnom_tunnel(atomic_number_);
    LocalFields Jion;
    double factorJion_0 = au_to_mec2 * EC_to_au * EC_to_au * invdt;

    int nparts = Epart->size() / 3;
    double *Ex = &((*Epart)[0 * nparts]);
    double *Ey = &((*Epart)[1 * nparts]);
    double *Ez = &((*Epart)[2 * nparts]);

    for (unsigned int ipart = ipart_min; ipart < ipart_max; ipart++) {
        // Current charge state of the ion
        Z = (unsigned int)(particles->charge(ipart));

        // If ion already fully ionized then skip
        if (Z == atomic_number_) {
            continue;
        }

        // Absolute value of the electric field normalized in atomic units
        E = EC_to_au * sqrt(*(Ex + ipart - ipart_ref) * *(Ex + ipart - ipart_ref) +
                            *(Ey + ipart - ipart_ref) * *(Ey + ipart - ipart_ref) +
                            *(Ez + ipart - ipart_ref) * *(Ez + ipart - ipart_ref));
        if (E < 1e-10) {
            continue;
        }

        // --------------------------------
        // Start of the Monte-Carlo routine
        // --------------------------------

        invE = 1. / E;
        factorJion = factorJion_0 * invE * invE;
        ran_p = patch->rand_->uniform();
        IonizRate_tunnel[Z] = ionizationRate(Z, E);

        // Total ionization potential (used to compute the ionization current)
        TotalIonizPot = 0.0;

        // k_times will give the nb of ionization events
        k_times = 0;
        Zp1 = Z + 1;

        if (Zp1 == atomic_number_) {
            // if ionization of the last electron: single ionization
            // -----------------------------------------------------
            if (ran_p < 1.0 - exp(-IonizRate_tunnel[Z] * dt)) {
                TotalIonizPot += Potential[Z];
                k_times = 1;
            }

        } else {
            // else : multiple ionization can occur in one time-step
            //        partial & final ionization are decoupled (see Nuter Phys.
            //        Plasmas)
            // -------------------------------------------------------------------------

            // initialization
            Mult = 1.0;
            Dnom_tunnel[0] = 1.0;
            Pint_tunnel = exp(-IonizRate_tunnel[Z] * dt);  // cummulative prob.

            // multiple ionization loop while Pint_tunnel < ran_p and still partial
            // ionization
            while ((Pint_tunnel < ran_p) and (k_times < atomic_number_ - Zp1)) {
                newZ = Zp1 + k_times;
                IonizRate_tunnel[newZ] = ionizationRate(newZ, E);
                D_sum = 0.0;
                P_sum = 0.0;
                Mult *= IonizRate_tunnel[Z + k_times];
                for (unsigned int i = 0; i < k_times + 1; i++) {
                    Dnom_tunnel[i] = Dnom_tunnel[i] / (IonizRate_tunnel[newZ] - IonizRate_tunnel[Z + i]);
                    D_sum += Dnom_tunnel[i];
                    P_sum += exp(-IonizRate_tunnel[Z + i] * dt) * Dnom_tunnel[i];
                }
                Dnom_tunnel[k_times + 1] = -D_sum;
                P_sum = P_sum + Dnom_tunnel[k_times + 1] * exp(-IonizRate_tunnel[newZ] * dt);
                Pint_tunnel = Pint_tunnel + P_sum * Mult;

                TotalIonizPot += Potential[Z + k_times];
                k_times++;
            }  // END while

            // final ionization (of last electron)
            if (((1.0 - Pint_tunnel) > ran_p) && (k_times == atomic_number_ - Zp1)) {
                TotalIonizPot += Potential[atomic_number_ - 1];
                k_times++;
            }
        }  // END Multiple ionization routine

        // Compute ionization current
        if (patch->EMfields->Jx_ != NULL) {  // For the moment ionization current is
                                             // not accounted for in AM geometry
            factorJion *= TotalIonizPot;
            Jion.x = factorJion * *(Ex + ipart);
            Jion.y = factorJion * *(Ey + ipart);
            Jion.z = factorJion * *(Ez + ipart);

            Proj->ionizationCurrents(patch->EMfields->Jx_, patch->EMfields->Jy_, patch->EMfields->Jz_, *particles, ipart, Jion);
        }

        // Creation of the new electrons
        // (variable weights are used)
        // -----------------------------

        if (k_times != 0) {
            new_electrons.createParticle();
            int idNew = new_electrons.size() - 1;
            for (unsigned int i = 0; i < new_electrons.dimension(); i++) {
                new_electrons.position(i, idNew) = particles->position(i, ipart);
            }
            for (unsigned int i = 0; i < 3; i++) {
                new_electrons.momentum(i, idNew) = particles->momentum(i, ipart) * ionized_species_invmass;
            }
            new_electrons.weight(idNew) = double(k_times) * particles->weight(ipart);
            new_electrons.charge(idNew) = -1;

            if (save_ion_charge_) {
                ion_charge_.push_back(particles->charge(ipart));
            }

            // Increase the charge of the particle
            particles->charge(ipart) += k_times;
        }

    }  // Loop on particles
}

template <int Tunneling_Model, int BSI>
inline double IonizationTunnel<Tunneling_Model, BSI>::ionizationRate(const int Z, const double E)
{
    if (BSI == 0) { // no barrier suppression, use tunneling formula
        double delta = gamma_tunnel[Z] / E;
        return beta_tunnel[Z] * exp(-delta * one_third + alpha_tunnel[Z] * log(delta));
    } else if (BSI == 1) { // Tong-Ling barrier suppression model
        const double delta = gamma_tunnel[Z] / E;
        return beta_tunnel[Z] * exp(-delta * one_third + alpha_tunnel[Z] * log(delta) - E * lambda_tunnel[Z]);
    } else if (BSI == 2) { // KAG barrier suppression model
        constexpr double IH = 13.598434005136;
        double ratio_of_IPs = IH / IonizationTables::ionization_energy(atomic_number_, Z);

        double BSI_rate_quadratic = 2.4 * (E * E) * ratio_of_IPs * ratio_of_IPs * au_to_w0;
        double BSI_rate_linear = 0.8 * E * sqrt(ratio_of_IPs) * au_to_w0;
        double delta = gamma_tunnel[Z] / E;
        double Tunnel_rate = beta_tunnel[Z] * exp(-delta / 3.0 + alpha_tunnel[Z] * log(delta));

        if (BSI_rate_quadratic >= BSI_rate_linear) {
            return BSI_rate_linear;
        } else if (std::min(Tunnel_rate, BSI_rate_quadratic) == BSI_rate_quadratic) {
            return BSI_rate_quadratic;
        } else {
            return Tunnel_rate;
        }
    }
}



#endif
