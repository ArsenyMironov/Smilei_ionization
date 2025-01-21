#include "IonizationTunnelCustomCoeffTL.h"
#include "IonizationTables.h"

#include <cmath>

#include "Particles.h"
#include "Species.h"

using namespace std;



IonizationTunnelCustomCoeffTL::IonizationTunnelCustomCoeffTL( Params &params, Species *species ) : Ionization( params, species )
{
    DEBUG( "Creating the Tunnel Ionizaton CustomCoeff + TL class" );
    
    atomic_number_          = species->atomic_number_;
    m_equal_zero_           = species->m_equal_zero_;
    use_g_factor_           = species->use_g_factor_;
    cnl_model_              = species->cnl_model_;
    ionization_tl_parameter_ = species->ionization_tl_parameter_; 

    
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
    lambda_tunnel.resize( atomic_number_ );
    
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
        lambda_tunnel[Z] = ionization_tl_parameter_*pow(cst, 2)/gamma_tunnel[Z];
    }

    DEBUG( "Finished Creating the Tunnel Ionizaton class" );
    
}



void IonizationTunnelCustomCoeffTL::operator()( Particles *particles, unsigned int ipart_min, unsigned int ipart_max, vector<double> *Epart, Patch *patch, Projector *Proj, int ipart_ref )
{

    unsigned int Z, Zp1, newZ, k_times;
    double TotalIonizPot, E, invE, factorJion, delta, ran_p, Mult, D_sum, P_sum, Pint_tunnel;
    vector<double> IonizRate_tunnel( atomic_number_ ), Dnom_tunnel( atomic_number_ );
    LocalFields Jion;
    double factorJion_0 = au_to_mec2 * EC_to_au*EC_to_au * invdt;
    
    int nparts = Epart->size()/3;
    double *Ex = &( ( *Epart )[0*nparts] );
    double *Ey = &( ( *Epart )[1*nparts] );
    double *Ez = &( ( *Epart )[2*nparts] );
 
    for( unsigned int ipart=ipart_min ; ipart<ipart_max; ipart++ ) {
    
        // Current charge state of the ion
        Z = ( unsigned int )( particles->charge( ipart ) );
        
        // If ion already fully ionized then skip
        if( Z==atomic_number_ ) {
            continue;
        }
        
        // Absolute value of the electric field normalized in atomic units
        E = EC_to_au * sqrt( pow( *( Ex+ipart-ipart_ref ), 2 )
                             +pow( *( Ey+ipart-ipart_ref ), 2 )
                             +pow( *( Ez+ipart-ipart_ref ), 2 ) );
        if( E<1e-10 ) {
            continue;
        }
        
        // --------------------------------
        // Start of the Monte-Carlo routine
        // --------------------------------
        
        invE = 1./E;
        factorJion = factorJion_0 * invE*invE;
        delta      = gamma_tunnel[Z]*invE;
        ran_p = patch->rand_->uniform();
        IonizRate_tunnel[Z] = beta_tunnel[Z] * exp( -delta*one_third + alpha_tunnel[Z]*log( delta ) - E*lambda_tunnel[Z] );

        // MESSAGE(0, "Z: " << Z << " I=" << Potential[Z] \
        //             << " n*=" << (alpha_tunnel[Z]+1)/2. \
        //             << " l=" << int(Azimuthal_quantum_number[Z]) \
        //             << " Gamma(n*-l)=" << tgamma( (alpha_tunnel[Z]+1)/2.-Azimuthal_quantum_number[Z] ) \
        //             << " C_n*l^2=" << pow( 2, alpha_tunnel[Z] ) * 4.0 / \
        //                   ( (alpha_tunnel[Z]+1)*tgamma( (alpha_tunnel[Z]+1)/2.+Azimuthal_quantum_number[Z]+1 )*tgamma( (alpha_tunnel[Z]+1)/2.-Azimuthal_quantum_number[Z]) ) \
        //             << " C_n*l^2[ADK]=" << pow( 2, alpha_tunnel[Z] ) * 4.0 / \
        //                   ( (alpha_tunnel[Z]+1)*tgamma(alpha_tunnel[Z]+1) ) \
        //             << " F=" << E/pow( 2.0*Potential[Z], 1.5 ) \
        //             << " W=" << IonizRate_tunnel[Z]);
        
        // Total ionization potential (used to compute the ionization current)
        TotalIonizPot = 0.0;
        
        // k_times will give the nb of ionization events
        k_times = 0;
        Zp1=Z+1;
        
        if( Zp1 == atomic_number_ ) {
            // if ionization of the last electron: single ionization
            // -----------------------------------------------------
            if( ran_p < 1.0 -exp( -IonizRate_tunnel[Z]*dt ) ) {
                TotalIonizPot += Potential[Z];
                k_times        = 1;
            }
            
        } else {
            // else : multiple ionization can occur in one time-step
            //        partial & final ionization are decoupled (see Nuter Phys. Plasmas)
            // -------------------------------------------------------------------------
            
            // initialization
            Mult = 1.0;
            Dnom_tunnel[0]=1.0;
            Pint_tunnel = exp( -IonizRate_tunnel[Z]*dt ); // cummulative prob.
            
            //multiple ionization loop while Pint_tunnel < ran_p and still partial ionization
            while( ( Pint_tunnel < ran_p ) and ( k_times < atomic_number_-Zp1 ) ) {
                newZ = Zp1+k_times;
                delta = gamma_tunnel[newZ]*invE;
                IonizRate_tunnel[newZ] = beta_tunnel[newZ]
                                            * exp( - delta*one_third + alpha_tunnel[newZ]*log( delta ) - E*lambda_tunnel[newZ] );
                D_sum = 0.0;
                P_sum = 0.0;
                Mult  *= IonizRate_tunnel[Z+k_times];
                for( unsigned int i=0; i<k_times+1; i++ ) {
                    Dnom_tunnel[i]=Dnom_tunnel[i]/( IonizRate_tunnel[newZ]-IonizRate_tunnel[Z+i] );
                    D_sum += Dnom_tunnel[i];
                    P_sum += exp( -IonizRate_tunnel[Z+i]*dt )*Dnom_tunnel[i];
                }
                Dnom_tunnel[k_times+1]  = -D_sum; // bug fix
                P_sum                   = P_sum + Dnom_tunnel[k_times+1]*exp( -IonizRate_tunnel[newZ]*dt );
                Pint_tunnel             = Pint_tunnel + P_sum*Mult;
                
                TotalIonizPot += Potential[Z+k_times];
                k_times++;
            }//END while
            
            // final ionization (of last electron)
            if( ( ( 1.0-Pint_tunnel )>ran_p ) && ( k_times==atomic_number_-Zp1 ) ) {
                TotalIonizPot += Potential[atomic_number_-1];
                k_times++;
            }
        }//END Multiple ionization routine
        
        // Compute ionization current
        if (patch->EMfields->Jx_ != NULL){  // For the moment ionization current is not accounted for in AM geometry
            factorJion *= TotalIonizPot;
            Jion.x = factorJion * *( Ex+ipart );
            Jion.y = factorJion * *( Ey+ipart );
            Jion.z = factorJion * *( Ez+ipart );
            
            Proj->ionizationCurrents( patch->EMfields->Jx_, patch->EMfields->Jy_, patch->EMfields->Jz_, *particles, ipart, Jion );
        }
        
        // Creation of the new electrons
        // (variable weights are used)
        // -----------------------------

        if( k_times !=0 ) {
            new_electrons.createParticle();
            int idNew = new_electrons.size() - 1;
            for( unsigned int i=0; i<new_electrons.dimension(); i++ ) {
                new_electrons.position( i, idNew )=particles->position( i, ipart );
            }
            for( unsigned int i=0; i<3; i++ ) {
                new_electrons.momentum( i, idNew ) = particles->momentum( i, ipart )*ionized_species_invmass;
            }
            new_electrons.weight( idNew )=double( k_times )*particles->weight( ipart );
            new_electrons.charge( idNew )=-1;
            
            if( save_ion_charge_ ) {
                ion_charge_.push_back( particles->charge( ipart ) );
            }
            
            // Increase the charge of the particle
            particles->charge( ipart ) += k_times;
        }
        
    } // Loop on particles
}


