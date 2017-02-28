#ifndef MASTERH
#include "MasterHeader.h"

#else

// Flags
#define GP_OK 0.
#define GP_OUT 1.
#define GP_PBC 2.
#define GP_LOST 16.
#define GP_BAD 255.

class Multisphere {

public :
  Multisphere ():initialized(false),ngp(0), currentstepinit(false), currentstep(-1) { for (int i=0; i<7; i++) symetrie[i]=false ; } 
  
  // Functions
  int init(Step & step) ; 
  int get_orientations(Step &step) ; 
  Matrix3d compute_K (Step &step) ; 
  double compute_dzeta (Step &step) ; 
  void compute_eigen(Step &step) ; 
  void set_current_step (int stepid) {if (stepid != currentstep) {currentstep=stepid ; currentstepinit=false ; } }
  void check() ; 
  
  // Variables
  int ngp ;  // Nombre de groupes
  vector < vector <int> > gps ; 
  vector <vector <double> > data ; //Values: flag, centroid x, y, z, orientation x, y, z.
  int type ; 
  
private:
  vector < Vector > pts, segments ; 
  bool initialized ; 
  bool symetrie[7] ; 
  int currentstep ; bool currentstepinit ;  
} ;
#endif