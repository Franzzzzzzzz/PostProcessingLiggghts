#ifndef MASTERH
#include "MasterHeader.h"

#else

class Surface {
public :
  // Fonctions
  int detect_surface (Dump & dump) ; 
  int detect_surface2D (Dump & dump) ;
  int write_asmatlab (string chemin) ; 
  int write_asmatlab2D (string chemin) ;
  int write_asmatlab_bottom (string chemin) ;
  
  //Variables
  vector <vector < vector <double> > > positions ;
  vector < vector <double> >  positions2D ;
  vector <vector < vector <double> > > positionsbas ;
  int nbsteps ; 
  int nbbox[2] ;
} ; 

#endif
