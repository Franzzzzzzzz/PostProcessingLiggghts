#ifndef MASTERH
#include "MasterHeader.h"

#else

// Flags
#define GP_OK 0
#define GP_OUT 1
#define GP_PBC 2
#define GP_LOST 16
#define GP_BAD 255

class Multisphere {

public :
  Multisphere ():initialized(false),ngp(0) {} 
  
  // Functions
  int init(Step & step) ; 
  int get_orientations(Step &step) ; 
  
  // Variables
  int ngp ;  // Nombre de groupes
  vector < vector <int> > gps ; 
  vector <vector <double> > data ; //Values: flag, centroid x, y, z, orientation x, y, z.
  int type ; 
  
  
private:
  vector < Vector > pts, segments ; 
  bool initialized ; 
} ;
#endif