#ifndef IonizationFactory_H
#define IonizationFactory_H

#include "Ionization.h"
#include "IonizationTunnel.h"
#include "IonizationFromRate.h"
#include "IonizationTunnelEnvelopeAveraged.h"
#include "IonizationTunnelBSI.h"
#include "IonizationTunnelTL.h"
#include "IonizationTunnelFullPPT.h"
#include "IonizationTunnelPPT.h"
#include "IonizationTunnelFullADK.h"
#include "IonizationTunnelCustomCoeff.h"
#include "IonizationTunnelCustomCoeffBSI.h"
#include "IonizationTunnelCustomCoeffTL.h"
#include "IonizationTunnelCustomCoeffG.h"

#include "Params.h"

#include "Tools.h"

#include "Species.h"

//! this class create and associate the right ionization model to species
class IonizationFactory
{
public:
    static Ionization *create( Params &params, Species *species )
    {
        Ionization *Ionize = NULL;
        std::string model=species->ionization_model_;
        
        if( model == "tunnel" ) {
            
            if( species->max_charge_ > ( int )species->atomic_number_ ) {
                ERROR( "Charge > atomic_number for species " << species->name_ );
            }

            if( params.Laser_Envelope_model ) {
                ERROR( "The ionization model for species interacting with envelope is tunnel_envelope_averaged" );
            }
            
            Ionize = new IonizationTunnel( params, species );
            
        } else if( model == "tunnel_envelope_averaged" ) {
            if( species->max_charge_ > ( int )species->atomic_number_ ) {
                ERROR( "Charge > atomic_number for species " << species->name_ );
            }
            if( species->particles->is_test ) {
                ERROR( "Cannot ionize test species " << species->name_ );
            }
            
            Ionize = new IonizationTunnelEnvelopeAveraged( params, species );
            
            if ( !params.Laser_Envelope_model ) {
                ERROR( "The ionization model tunnel_envelope_averaged needs a laser envelope");
            }

         } else if( model == "from_rate" ) {
            
            if( species->max_charge_ > ( int )species->maximum_charge_state_ ) {
                ERROR( "For species '" << species->name_ << ": charge > maximum_charge_state" );
            }
            
            Ionize = new IonizationFromRate( params, species );
  
        } else if(model == "tunnel_BSI") { // added by I. Ouatu.  Put keyword "tunnel_BSI" for the 
                                           // Tong-Lin ionization model in the species description in your namelist.
                                           // The alpha parameter is specified by "ionization_tl_parameter" in the atom species 

            if (species->max_charge_ > (int) species->atomic_number_) { // same as for simple Tunnel Ionization.
                ERROR( "Charge > atomic_number for species " << species->name_);
            }
            if( (params.Laser_Envelope_model) ) { // same as for simple Tunnel Ionziation
                ERROR( "The ionization model for species interacting with envelope is tunnel_envelope_averaged" );
            }

            Ionize = new IonizationTunnelBSI(params, species);

        } else if(model == "tunnel_TL") { // added by Arseny Mironov. Put keyword "tunnel_TL" for the Tong-Lin ionization model
                                          // in the species description in your namelist.
            if (species->max_charge_ > (int) species->atomic_number_) {
                ERROR( "Charge > atomic_number for species " << species->name_);
            }
            if( (params.Laser_Envelope_model) ) { // same as for simple Tunnel Ionziation
                ERROR( "The ionization model for species interacting with envelope is tunnel_envelope_averaged" );
            }
            Ionize = new IonizationTunnelTL(params, species);    

        } else if(model == "tunnel_full_PPT") { // added by Arseny Mironov. Put keyword "tunnel_full_PPT" for the tunneling ionization model 
                                                // with account for the magnetic quantum number in the species description in your namelist.
                                          
            if( species->max_charge_ > ( int )species->atomic_number_ ) {
                ERROR( "Charge > atomic_number for species " << species->name_ );
            }

            if( params.Laser_Envelope_model ) {
                ERROR( "The ionization model for species interacting with envelope is tunnel_envelope_averaged" );
            }

            Ionize = new IonizationTunnelFullPPT( params, species );

        } else if(model == "tunnel_PPT") { // added by Arseny Mironov. Put keyword "tunnel_PPT" for the tunneling ionization model 
                                                // with account for l instead of l* in A_nl.
                                          
            if( species->max_charge_ > ( int )species->atomic_number_ ) {
                ERROR( "Charge > atomic_number for species " << species->name_ );
            }

            if( params.Laser_Envelope_model ) {
                ERROR( "The ionization model for species interacting with envelope is tunnel_envelope_averaged" );
            }

            Ionize = new IonizationTunnelPPT( params, species );

        } else if(model == "tunnel_full_ADK") { // added by Arseny Mironov. Put keyword "tunnel_full_ADK" for the tunneling ionization model 
                                              // with account for m in the ADK formula.
                                          
            if( species->max_charge_ > ( int )species->atomic_number_ ) {
                ERROR( "Charge > atomic_number for species " << species->name_ );
            }

            if( params.Laser_Envelope_model ) {
                ERROR( "The ionization model for species interacting with envelope is tunnel_envelope_averaged" );
            }

            Ionize = new IonizationTunnelFullADK( params, species );
        } else if(model == "tunnel_custom_coefficient") { // added by Arseny Mironov. Put keyword "tunnel_custom_coefficient" for the tunneling ionization model 
                                                          // with customized Cnl coefficient and choice between m=0 or m!=0 for electron states.
                                                          // The coefficients are customized with 3 input params: cnl_model, cnl_squared_table (if needed), m_equal_zero
                                                          // Input options for Cnl coefficients:
                                                          // cnl_model = 0 for ADK (default)
                                                          // cnl_model = 1 for PPT (Hartree)
                                                          // cnl_model = 2 for custom table provided in cnl_squared_table
                                                          // cnl_squared_table is a list of Cnl^2 values, should be an array of doubles of lenght = atomic_number
                                                          // Input for m:
                                                          // m_equal_zero = True: m=0 for all electrons during ionization (default)
                                                          // m_equal_zero = False: m!=0 and is taken from IonizationTables::magnetic_atomic_number
            if( species->max_charge_ > ( int )species->atomic_number_ ) {
                ERROR( "Charge > atomic_number for species " << species->name_ );
            }

            if( params.Laser_Envelope_model ) {
                ERROR( "The ionization model for species interacting with envelope is tunnel_envelope_averaged" );
            }

            Ionize = new IonizationTunnelCustomCoeff( params, species );
        } else if(model == "tunnel_custom_coefficient_BSI") { // added by Arseny Mironov. Same as "tunnel_custom_coefficient", 
                                                              // but added BSI as in the Ouatu realization
            if( species->max_charge_ > ( int )species->atomic_number_ ) {
                ERROR( "Charge > atomic_number for species " << species->name_ );
            }

            if( params.Laser_Envelope_model ) {
                ERROR( "The ionization model for species interacting with envelope is tunnel_envelope_averaged" );
            }

            Ionize = new IonizationTunnelCustomCoeffBSI( params, species );
        } else if(model == "tunnel_custom_coefficient_TL") { // added by Arseny Mironov. Same as "tunnel_custom_coefficient", 
                                                              // but added Tong-Lin suppression factor
            if( species->max_charge_ > ( int )species->atomic_number_ ) {
                ERROR( "Charge > atomic_number for species " << species->name_ );
            }

            if( params.Laser_Envelope_model ) {
                ERROR( "The ionization model for species interacting with envelope is tunnel_envelope_averaged" );
            }

            Ionize = new IonizationTunnelCustomCoeffTL( params, species );
        } else if(model == "tunnel_custom_coefficient_G") { 
            if( species->max_charge_ > ( int )species->atomic_number_ ) {
                ERROR( "Charge > atomic_number for species " << species->name_ );
            }

            if( params.Laser_Envelope_model ) {
                ERROR( "The ionization model for species interacting with envelope is tunnel_envelope_averaged" );
            }

            Ionize = new IonizationTunnelCustomCoeffG( params, species );
        }
        
        return Ionize;
    }

};

#endif
