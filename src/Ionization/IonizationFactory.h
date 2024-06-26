#ifndef IonizationFactory_H
#define IonizationFactory_H

#include "Ionization.h"
#include "IonizationTunnel.h"
#include "IonizationFromRate.h"
#include "IonizationTunnelEnvelopeAveraged.h"
#include "IonizationTunnelBSI.h"
#include "IonizationTunnelTL.h"
#include "IonizationTunnelFullPPT.h"

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
        }
        
        return Ionize;
    }

};

#endif
