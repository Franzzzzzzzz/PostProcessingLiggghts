#ifndef MASTERH
#include "MasterHeader.h"

#else
// Classe commune à tous les type de dump
class Dump {
public : 
   // Variables
   int nbsteps ;
   gzifstream dumpin ; 
   long int inscrit ;   
   vector <Step> steps ;
   vector <struct Op> prefiltre ;
   vector <struct Op> postfiltre ;
   Multisphere multisphere ; 
   int Type ; 
   
   // Fonctions
   Dump (void) ;  
   virtual void disp(void) ; 
   int open (string chemin) ;
   void loopdat (long int *res) ;
   int write_asVTK (string chemin) {return(write_asVTKBase(chemin, 0));}
   int write_asOneVTK (string chemin) {return(write_asVTKBase(chemin, 1)) ; }
   int check_timestep (long int timestep);
   virtual int free_timestep (long int inscrit) ;
protected :
   virtual int read (unsigned char read_type, int index) = 0 ;
   int write_asVTKBase (string chemin, char type) ; 
};



#endif

