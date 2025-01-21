#include "IonizationTunnelCustomCoeffBSI.h"
#include "IonizationTables.h"
#include "continuity_tool.h"

#include <cmath>

#include "Particles.h"
#include "Species.h"
using namespace std;

IonizationTunnelCustomCoeffBSI::IonizationTunnelCustomCoeffBSI( Params &params, Species *species ) : Ionization( params, species )
{
    DEBUG( "Creating the Tunnel Ionizaton CustomCoeff + BSI class" );
    
    atomic_number_          = species->atomic_number_;
    m_equal_zero_           = species->m_equal_zero_;
    use_g_factor_           = species->use_g_factor_;
    cnl_model_              = species->cnl_model_;

    
    // Ionization potential & quantum numbers (all in atomic units 1 au = 27.2116 eV)
    Potential.resize( atomic_number_ );
    Azimuthal_quantum_number.resize( atomic_number_ );
    for( int Zstar=0; Zstar<( int )atomic_number_; Zstar++ ) {
        Potential               [Zstar] = IonizationTables::ionization_energy( atomic_number_, Zstar ) * eV_to_au;
        Azimuthal_quantum_number[Zstar] = IonizationTables::azimuthal_atomic_number( atomic_number_, Zstar );
    }

    g_factor.resize( atomic_number_ );
    std::fill(g_factor.begin(), g_factor.end(), 1);

    if( !m_equal_zero_ ) {
        Magnetic_quantum_number.resize( atomic_number_ );
        for( int Zstar=0; Zstar<( int )atomic_number_; Zstar++ ) {
            Magnetic_quantum_number[Zstar] = IonizationTables::magnetic_atomic_number( atomic_number_, Zstar );
        }
        if( use_g_factor_ ) {
            for( int Zstar=0; Zstar<( int )atomic_number_; Zstar++ ) {
                for( int i=Zstar+1; i<( int )atomic_number_; i++ ) {
                    if( (abs(Magnetic_quantum_number[Zstar]) == abs(Magnetic_quantum_number[i])) ) {
                        g_factor[Zstar] += 1;
                    } else {
                        break;
                    }
                }
            }
        }
    }

    if( cnl_model_ == 2 ) {
        cnl_squared_table_          = species->cnl_squared_table_;
    };
    
    for( unsigned int i=0; i<atomic_number_; i++ ) {
        DEBUG( "ioniz: i " << i << " potential: " << Potential[i] << " Az.q.num: " << Azimuthal_quantum_number[i] );
    }
    
    one_third = 1.0/3.0;
    
    alpha_tunnel.resize( atomic_number_ );
    beta_tunnel.resize( atomic_number_ );
    gamma_tunnel.resize( atomic_number_ );
    
    for( unsigned int Z=0 ; Z<atomic_number_ ; Z++ ) {
        DEBUG( "Z : " << Z );
        double cst      = ( ( double )Z+1.0 ) * sqrt( 2.0/Potential[Z] );
        double Cnl      = 1.;
        double Blm      = 1.;
        double abs_m    = 0.;

        if( m_equal_zero_ ) {
            Blm      = ( 2.*Azimuthal_quantum_number[Z]+1.0 );
        } else {
            abs_m    = abs(Magnetic_quantum_number[Z]);
            Blm      = ( 2.*Azimuthal_quantum_number[Z]+1.0 ) * \
                       tgamma(Azimuthal_quantum_number[Z]+abs_m+1) / ( tgamma(abs_m+1)*tgamma(Azimuthal_quantum_number[Z]-abs_m+1) );
        }

        alpha_tunnel[Z] = cst-1.0-abs_m;

        if( cnl_model_ == 0 ) {
            Cnl = pow( 2, alpha_tunnel[Z] ) / \
                        ( cst*tgamma( cst ) );
        } else if( cnl_model_ == 1 ) {
            if( Z>0 ) {
                Cnl = pow( 2, alpha_tunnel[Z] ) / \
                        ( cst*tgamma( cst/2.0+Azimuthal_quantum_number[Z]+1 )*tgamma( cst/2.0-Azimuthal_quantum_number[Z]) );
            }
        } else {
            Cnl = cnl_squared_table_[Z];
        }

        MESSAGE(0, "Z=" << Z << ", l=" << Azimuthal_quantum_number[Z] << ", |m|=" << abs_m << ", g=" << g_factor[Z] << ", Cnl=" << Cnl << ", Blm=" << Blm);

        beta_tunnel[Z] = g_factor[Z]*Cnl*4.*Blm * Potential[Z] * au_to_w0;
        gamma_tunnel[Z] = 2.0 * pow( 2.0*Potential[Z], 1.5 );
    }

    DEBUG( "Finished Creating the Tunnel Ionizaton class" );
    
}



// override the function call operator ().

void IonizationTunnelCustomCoeffBSI::operator()( Particles *particles, unsigned int ipart_min, unsigned int ipart_max, vector<double> *Epart, Patch *patch, Projector *Proj, int ipart_ref ) {
    unsigned int Z, Zp1, newZ, k_times;
    double TotalIonizPot, E, invE, factorJion, delta, ran_p, Mult, D_sum, P_sum, Pint_tunnel, Pint_BSI_linear, Pint_BSI_quadratic;
    vector<double> IonizRate_tunnel(atomic_number_), Dnom_tunnel(atomic_number_); // For MonteCarlo procedure using the tunnel ioniz. rate.
    vector<double> IonizRate_BSI_linear(atomic_number_),    Dnom_BSI_linear(atomic_number_); // For MonteCarlo procedure using the BSI ioniz. rate.
    vector<double> IonizRate_BSI_quadratic(atomic_number_), Dnom_BSI_quadratic(atomic_number_);
    LocalFields Jion; // To deposit ionization current appearing from the ionization events.
    double factorJion_0 = au_to_mec2 * EC_to_au*EC_to_au * invdt; // Will be used to calc. ionization current.
    
    int nparts = Epart->size()/3; // Over 3 because there are 3 Dims. Epart is like [part1_Ex, part2_Ex, ... partnparts_Ex, part1_Ey, part2_Ey,...partnparts_Ey, part1_Ez, part2_Ez, partnparts_Ez]}
    double *Ex = &( (*Epart)[0*nparts] ); // Pointer to beginning of the vector which holds the (3 * n_parts) values of the n_parts E-field (Physics) 3D vectors.
    double *Ey = &( (*Epart)[1*nparts] ); // Pointer to where E_y values start inside the Epart vector.
    double *Ez = &( (*Epart)[2*nparts] ); // Pointer to where E_z values start inside the Epart vector. 

    for(unsigned int ipart = ipart_min ; ipart < ipart_max; ipart++) { // Loop on particles from Particles* particles input argument
        // Current charge state of the ion.
        Z = (unsigned int) (particles->charge(ipart)) ; // current charge state of the particle number ipart of species called ''species'' with atomic no = atomic_number_
        Zp1 = Z + 1;
        // Implementation of Barrier Suppresion Ionization Critical Field E_cr.
        // need E_cr to be in atomic units. Thus need to implement E_cr as from DonnaStrickland1991, page 5, eq(5). [Laser ionization of noble gases by Coulomb-barrier suppression]
        // There the Ionization Potential appearing in formula is in AtomicUnits. 
        // Here, all values from Potential[] vector are in AtomicUnits also.
        // Thus need to implement exactly Donna's formula: E_cr = pow(Potential[atomic_number*(atomic_number_-1)/2 + Z] , 2) / (4*Z);
        // Donna's E_cr easily matches the SI units version of E_cr presented in Artmenko, Kostyukov, PhysRevA 96, 0321,06 (2017), page 2, section II.
        // Russians researchers say E_cr = atomic_unit_of_field * kappa^4 / (16*Z), where kappa^2 = Ii/IH, thus Ii is measured in units of (atomic_unit_of_energy / 2).
        // Had they measured Ii in atomic_unit_of_energy, the 16 on the denominator would become a 4. (atomic_unit_of_energy = 2*IH).
        // Thus the implementation of E_cr in atomic units, following the Russians researchers would be:
        // double ionizpots_ratio_forkappa = (Potential[ atomic_number_*(atomic_number_-1)/2 + Z ] / Potential[0]); // Potential[atomic_number_*(atomic_number_-1)/2 + Z] is IonizEnergy for current charge state of ion with atomic number equal to atomic_number_. Potential[0] is IonizPot of Hydrogen H. All the values in Pontential[] vector are in AtomicUnits. Thus the normalization dissapears when taking ratio of Ii/IH.
        // double kappa = pow(ionizpots_ratio_forkappa , 0.5);
        // double E_cr =  ( pow(kappa,4) / (16*Z) ) / atomic_unit_of_field; // E_cr in AtomicUnits, as from Artmenko, Kostyukov, PhysRevA 96, 0321,06 (2017), page 2, section II.
        
        // HOW IS E_CR USEFUL? Answer: it is useless at the moment.
        // E_cr in atomic units formula (with Potential[] vector's numerical values being in AtomicUnits):
        double E_cr = pow( Potential[Z], 2) / (4*Zp1); //  Formula from DonnaStrickland, 1991, page5, eq(5) [Laser ionization of noble gases by Coulomb-barrier suppression]
        // There's absolutely no usage of this E_cr in the code.

        // If ion already fully ionized then skip and go to next particle
        if (Z==atomic_number_) { continue; } // atom fully ionized already, no electrons remained to be stripped away.

        // Absolute value of the electric field in AtomicUnits (i.e. normalized to (in units of) E_atomic = 5.1422*10^11 V/m) 
        E = EC_to_au * sqrt( pow( *(Ex + ipart - ipart_ref), 2) // (Ex+ipart-ipart_ref) points to the location in the container Epart at which E-field_x for particle ipart sits.
                             + pow( *(Ey + ipart - ipart_ref), 2) // Similarly for y. Dereferencing it via *() means we access and use the value sitting at that location.
                             + pow( *(Ez + ipart - ipart_ref), 2) ); // EC_to_au transforms from SMILEI-normalized units to AtomicUnits.
        if(E < 1e-10) { // E is in Atomic Units
            continue;
        }

        // Check whether to follow MC routine with Tunnel rate or with BSI rate (Quadratic or Linear).
        // Prof. Tony Arber in: Arber_2015_Plasma_Phys._Control._Fusion_57_113001 gives a hint on how to do this.
        // There he explains how to make the PIC code choose, at each timestep, between Tunnel and BSI (2 piece rate).
        
        // I am doing it in a simpler way, see int continuity_tool() function implementation!
        // int continuity_tool(parameter1, ... ,) function does the job. Depending on the output of continuity_tool(,...,), SMILEI chooses which rate to use at each timestep.

        int rate_idx = continuity_tool(Zp1, E, alpha_tunnel[Z], beta_tunnel[Z], gamma_tunnel[Z], E_cr, Potential[Z], atomic_number_, au_to_w0); //
        if (rate_idx != 0 ) { // Then use one of the two BSI rates: vector<double> IonizRate_BSI_linear(atomic_number_) or vector<double> IonizRate_BSI_quadratic(atomic_number_).
            // So rate is not Tunnel. The rate_idx can only be 1 or 2.
            if (rate_idx == 1) { // then use BSI-Quadratic rate, as dictated by continuity_tool()

                double IH = 13.598434005136; // IP of atomic Hydrogen (H), in eV.
                double Ii = IonizationTables::ionization_energy(atomic_number_, Z); // returns IP of Element with atomic no = atomic_number_ and charge state Z (Z=0 is a neutral atom)
                double ratio_of_IPs = IH / Ii;
                // IonizRate_BSI_quadratic[Z] = 2.4 * (pow(E,2)) / pow(Zp1,4); // From Bauer1999: Quadratic dep. on E and scales as 1/Z^4. From the formula on this line of code, the rate emerges in AtomicUnits, as E is in Atomic Units and Zp1 is a integer.
                IonizRate_BSI_quadratic[Z] = au_to_w0 * (2.4 * (pow(E,2)) * pow(ratio_of_IPs, 2)); // 1/(Z^4) = (IH_Ii)^2
                // IonizRate_BSI_quadratic as 1 line above is now returned in SMILEI units.
                invE = 1./E;
                factorJion = factorJion_0 * invE*invE; // Used to compute the ionization current Jion.
                // Monte Carlo routine from Nuter et al. 2011
                ran_p = patch -> rand_ -> uniform();
            
                TotalIonizPot = 0.0; // Used to compute the ionization current Jion.
        
                // k_times will give the number of ionization events in this timestep.
                k_times = 0;
                Zp1 = Z + 1;

                if (Zp1 == atomic_number_) { // single-ionization, i.e. ioniz. of the last electron.
                    if ( ran_p < 1 - exp (-IonizRate_BSI_quadratic[Z]*dt) ) { // MC dictates to ionize it
                        TotalIonizPot += Potential[Z];
                        k_times += 1;
                    } // else nothing happens, no ionization of the last electron during this timestep. for loop moves to next particle, i.e. the one indexed by ipart+1.
                }
                else { // multiple ionization can happen in one timestep because we are not ionizing the last electron yet.
                    //    partial & final ionization are decoupled (see Nuter Phys. Plasmas 2011)
                    // -------------------------------------------------------------------------
                
                    // initialization
                    Mult = 1.0;
                    Dnom_BSI_quadratic[0] = 1.0;
                    Pint_BSI_quadratic = exp( -IonizRate_BSI_quadratic[Z]*dt );

                    // multiple ionization loop while Pint_BSI < ran_p and still partial ionization
                    while( (Pint_BSI_quadratic < ran_p) and (k_times < atomic_number_-Zp1) ) {
                        newZ = Zp1 + k_times;
                        double Ii_newZ = IonizationTables::ionization_energy(atomic_number_, newZ);
                        double ratio_of_IPs_newZ = IH / Ii_newZ;
                        IonizRate_BSI_quadratic[newZ] = au_to_w0 * (2.4 * (pow(E,2)) * pow(ratio_of_IPs_newZ, 2));
                        D_sum = 0.0;
                        P_sum = 0.0;
                        Mult *= IonizRate_BSI_quadratic[ Z + k_times];
                        for (unsigned int i=0; i <= k_times+1; i++) {
                            Dnom_BSI_quadratic[i] = Dnom_BSI_quadratic[i] / (IonizRate_BSI_quadratic[newZ] - IonizRate_BSI_quadratic[Z+i]) ;
                            D_sum += Dnom_BSI_quadratic[i];
                            P_sum += exp(-IonizRate_BSI_quadratic[Z+i]*dt) * Dnom_BSI_quadratic[i];
                        }
                        Dnom_BSI_quadratic[k_times+1] -= D_sum;
                        P_sum                = P_sum + Dnom_BSI_quadratic[k_times+1] * exp(-IonizRate_BSI_quadratic[newZ]*dt);
                        Pint_BSI_quadratic             = Pint_BSI_quadratic + P_sum*Mult;

                        TotalIonizPot += Potential[Z + k_times];
                        k_times++; 
                    } // END while

                    // final ionization (of last electron)
                    if ( ((1.0-Pint_BSI_quadratic) > ran_p) && (k_times == atomic_number_-Zp1) ) {
                        TotalIonizPot += Potential[atomic_number_-1];
                        k_times++;
                    } // END final ionization (of last electron)
                
                } // END multiple ionization routine.

                // Compute Ionization Currents
                if (patch->EMfields->Jx_ != NULL) {
                    factorJion *= TotalIonizPot;
                    Jion.x = factorJion *  *(Ex + ipart); 
                    Jion.y = factorJion *  *(Ey + ipart);
                    Jion.z = factorJion *  *(Ez + ipart);
                    Proj->ionizationCurrents(patch->EMfields->Jx_, patch->EMfields->Jy_, 
                                         patch->EMfields->Jz_, *particles, ipart, Jion);
                }
                // Creation of the new electrons (variable weights are used)
                if(k_times != 0) { // If ionization has taken place.
                    new_electrons.createParticle();
                    int idNew = new_electrons.size() - 1; // new_electrons.size() =  number of electrons in new_electrons. -1 from idNew because the first elem of array has index 0.
                    for(unsigned int i=0; i < new_electrons.dimension(); i++) { // new_electrons.dimension() is 1,2 or 3, depending on geometry of simulation.
                        new_electrons.position(i, idNew) = particles->position(i, ipart);
                    }
                    for(unsigned int i=0; i<3; i++) {
                        new_electrons.momentum(i, idNew) = particles->momentum(i, ipart) * ionized_species_invmass;
                    }
                    new_electrons.weight(idNew) = double(k_times) * particles -> weight(ipart);
                    new_electrons.charge(idNew) = -1;

                    // Increase charge of particle.
                    particles->charge(ipart) += k_times;
                } // END creation of new electrons.

            } // END subroutine for BSI-Quadratic rate  // END if(index == 1) {}

            else { // if (rate_idx==1) {...} loop's else statement: means that rate_idx = 2 ==> use BSI-Linear
                double IH = 13.598434005136; // IP of atomic Hydrogen (H), in eV.
                double Ii = IonizationTables::ionization_energy(atomic_number_, Z);
                double ratio_of_IPs = IH / Ii;
                IonizRate_BSI_linear[Z] = au_to_w0 * (0.8 * E * pow(ratio_of_IPs, 0.5)); // paranthesis is in AtomicUnits
                invE = 1./E;
                factorJion = factorJion_0 * invE*invE; // Used to compute the ionization current Jion.
                
                ran_p = patch -> rand_ -> uniform();

                TotalIonizPot = 0.0; // Used to compute the ionization current Jion.

                // k_times will give the number of ionization events in this timestep.
                k_times = 0;
                Zp1 = Z + 1;
                if (Zp1 == atomic_number_) { // single-ionization, i.e. ioniz. of the last electron.
                    if ( ran_p < 1 - exp (-IonizRate_BSI_linear[Z]*dt) ) { // MC dictates to ionize it
                        TotalIonizPot += Potential[Z];
                        k_times += 1;
                    } // else nothing happens, no ionization of the last electron during this timestep. for loop moves to next particle, i.e. the one indexed by ipart+1.
                }
                else { // multiple ionization can happen in one timestep because we are not ionizing the last electron yet.
                    //    partial & final ionization are decoupled (see Nuter Phys. Plasmas 2011)
                    // -------------------------------------------------------------------------
                
                    // initialization
                    Mult = 1.0;
                    Dnom_BSI_linear[0] = 1.0;
                    Pint_BSI_linear = exp( -IonizRate_BSI_linear[Z]*dt );

                    // multiple ionization loop while Pint_BSI < ran_p and still partial ionization
                    while( (Pint_BSI_linear < ran_p) and (k_times < atomic_number_-Zp1) ) {
                        newZ = Zp1 + k_times;
                        double Ii_newZ = IonizationTables::ionization_energy(atomic_number_, newZ);
                        double ratio_of_IPs_newZ = IH / Ii_newZ;
                        IonizRate_BSI_linear[newZ] = au_to_w0 * (0.8 * E * pow(ratio_of_IPs_newZ, 0.5)); // Paranthesis is in Atomic Units
                        // Linear BSI rate as 1 line above is in SMILEI units.
                        D_sum = 0.0;
                        P_sum = 0.0;
                        Mult *= IonizRate_BSI_linear[ Z + k_times];
                        for (unsigned int i=0; i <= k_times+1; i++) {
                            Dnom_BSI_linear[i] = Dnom_BSI_linear[i] / (IonizRate_BSI_linear[newZ] - IonizRate_BSI_linear[Z+i]) ;
                            D_sum += Dnom_BSI_linear[i];
                            P_sum += exp(-IonizRate_BSI_linear[Z+i]*dt) * Dnom_BSI_linear[i];
                        }
                        Dnom_BSI_linear[k_times+1] -= D_sum;
                        P_sum                = P_sum + Dnom_BSI_linear[k_times+1] * exp(-IonizRate_BSI_linear[newZ]*dt);
                        Pint_BSI_linear             = Pint_BSI_linear + P_sum*Mult;

                        TotalIonizPot += Potential[Z + k_times];
                        k_times++; 
                    } // END while

                    // final ionization (of last electron)
                    if ( ((1.0-Pint_BSI_linear) > ran_p) && (k_times == atomic_number_-Zp1) ) {
                        TotalIonizPot += Potential[atomic_number_-1];
                        k_times++;
                    } // END final ionization (of last electron)
                
                } // END multiple ionization routine.

                // Compute Ionization Currents
                if (patch->EMfields->Jx_ != NULL) {
                    factorJion *= TotalIonizPot;
                    Jion.x = factorJion * *(Ex + ipart); 
                    Jion.y = factorJion * *(Ey + ipart);
                    Jion.z = factorJion * *(Ez + ipart);
                    Proj->ionizationCurrents(patch->EMfields->Jx_, patch->EMfields->Jy_, 
                                         patch->EMfields->Jz_, *particles, ipart, Jion);
                }

                // Creation of the new electrons (variable weights are used)
                if(k_times != 0) { // If ionization has taken place.
                    new_electrons.createParticle();
                    int idNew = new_electrons.size() - 1; // new_electrons.size() =  number of electrons in new_electrons. -1 from idNew because the first elem of array has index 0.
                    for(unsigned int i=0; i < new_electrons.dimension(); i++) { // new_electrons.dimension() is 1,2 or 3, depending on geometry of simulation.
                        new_electrons.position(i, idNew) = particles->position(i, ipart);
                    }
                    for(unsigned int i=0; i<3; i++) {
                        new_electrons.momentum(i, idNew) = particles->momentum(i, ipart) * ionized_species_invmass;
                    }
                    new_electrons.weight(idNew) = double(k_times) * particles -> weight(ipart);
                    new_electrons.charge(idNew) = -1;

                    // Increase charge of particle.
                    particles->charge(ipart) += k_times;
                } // END creation of new electrons.
            } // END else {// rate_idx = 2}
        } // END if(rate_idx != 0) {}

        else { // rate_idx = 0, so use Tunnel rate. There's no other rate_idx value which can appear. It is deemed to be 0 if the code reaches this else statement.
            invE = 1./E; 
            factorJion = factorJion_0 * invE*invE; 
            delta      = gamma_tunnel[Z] * invE; 
            ran_p = patch->rand_->uniform(); // 
            IonizRate_tunnel[Z] = beta_tunnel[Z] * exp(-delta*one_third + alpha_tunnel[Z]*log(delta));
            // IonizRate_tunnel as from 1 line above is in SMILEI units (beta has at its implementation's end the multiplication by au_to_w0)
            // Total ionization potential (used to compute the ionization current)
            TotalIonizPot = 0.0;
        
            // k_times will give the number of ionization events
            k_times = 0;
            Zp1 = Z + 1;
        
            if( Zp1 == atomic_number_ ) {
            // if ionization of the last electron: single ionization
            // -----------------------------------------------------
                if( ran_p < 1.0 - exp(-IonizRate_tunnel[Z]*dt) ) { // MC dictates to ionize it.
                    TotalIonizPot += Potential[Z];
                    k_times        = 1; // ionize it
                }
            } 
            else {
                // else : MULTIPLE IONZIATION can occur in one time-step
                //        partial & final ionization are decoupled (see Nuter Phys. Plasmas 2011)
                // -------------------------------------------------------------------------
            
                // initialization
                Mult = 1.0;
                Dnom_tunnel[0] = 1.0;
                Pint_tunnel = exp(-IonizRate_tunnel[Z]*dt); // cummulative prob. (Pint_tunnel is a double)
            
                // MULTIPLE IONZIATION loop while Pint_tunnel < ran_p and still partial ionization.
                while( (Pint_tunnel < ran_p) and (k_times < atomic_number_-Zp1) ) {
                    newZ = Zp1 + k_times;
                    delta = gamma_tunnel[newZ] * invE;
                    IonizRate_tunnel[newZ] = beta_tunnel[newZ] * exp(-delta*one_third + alpha_tunnel[newZ]*log(delta));
                    D_sum = 0.0;
                    P_sum = 0.0;
                    Mult  *= IonizRate_tunnel[ Z+k_times ];
                    for(unsigned int i=0;  i < k_times+1; i++) {
                        Dnom_tunnel[i] = Dnom_tunnel[i] / (IonizRate_tunnel[newZ] - IonizRate_tunnel[Z+i]);
                        D_sum += Dnom_tunnel[i];
                        P_sum += exp( -IonizRate_tunnel[Z+i]*dt ) * Dnom_tunnel[i];
                    }
                    Dnom_tunnel[k_times+1]  = -D_sum;  //bug fix
                    P_sum                   = P_sum + Dnom_tunnel[k_times+1] * exp(-IonizRate_tunnel[newZ]*dt);
                    Pint_tunnel             = Pint_tunnel + P_sum*Mult;
                
                    TotalIonizPot += Potential[ Z+k_times ];
                    k_times++;
                } //END while
            
                // final ionization (of last electron)
                if( ((1.0-Pint_tunnel) > ran_p) && (k_times == atomic_number_-Zp1) ) {
                    TotalIonizPot += Potential[atomic_number_-1];
                    k_times++;
                }
            } // END Multiple ionization routine for Tunnel rate
        
            // Compute ionization current
            if (patch->EMfields->Jx_ != NULL) {  // For the moment ionization current is not accounted for in AM geometry
                factorJion *= TotalIonizPot;
                Jion.x = factorJion * *(Ex + ipart);
                Jion.y = factorJion * *(Ey + ipart);
                Jion.z = factorJion * *(Ez + ipart);
            
                Proj->ionizationCurrents(patch->EMfields->Jx_, patch->EMfields->Jy_, patch->EMfields->Jz_, 
                                    *particles, ipart, Jion);
            }  
        
            // Creation of the new electrons (variable weights are used)
            // -----------------------------
            if( k_times !=0 ) { // If ionization has taken place.
                new_electrons.createParticle(); // new_electrons is an instance of class Particles. it is defined in Ionization.h
                //new_electrons.initialize( new_electrons.size()+1, new_electrons.dimension() ); // ???
                int idNew = new_electrons.size() - 1; // new_electrons.size() =  number of electrons in new_electrons. -1 from idNew because the first elem of array has index 0.

                for(unsigned int i=0; i<new_electrons.dimension(); i++) { // new_electrons.dimension() is 1,2 or 3, depending on geometry of simulation.
                    new_electrons.position(i, idNew) = particles->position(i, ipart);
                }
                for(unsigned int i=0; i<3; i++) { // Momentum is a 3Dim vector. (3Velocity simulation)
                    new_electrons.momentum(i, idNew) = particles->momentum(i, ipart) * ionized_species_invmass;
                }

                new_electrons.weight(idNew) = double(k_times) * particles->weight(ipart);
                new_electrons.charge(idNew) = -1;
            
                // Increase the charge of the particle
                particles->charge(ipart) += k_times;

            } // END creation of electrons for when using Tunnel rate
        } // END whole routine for when rate_idx = 0 , i.e. when using Tunnel rate

    } // END loop on particles

} // void IonizationTunnelCustomCoeffBSI::operator()(arg1, arg2, ...) scope end.
