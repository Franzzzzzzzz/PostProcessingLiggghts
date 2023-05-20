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
  Multisphere() ;

  // Functions
  int init(Step & step) ;
  int init_flux(Step & step) ;
  int get_orientations(Step &step) {if (isflux) get_orientations_flux(step) ; else get_orientations_noflux(step) ; return 0 ;}
  int get_orientations_noflux(Step &step) ;
  int get_orientations_flux (Step & step) ;
  Matrix3d compute_K (Step &step) ;
  double compute_dzeta (Step &step) ;
  void compute_eigen(Step &step) ;
  void set_current_step (int stepid) {if (stepid != currentstep) {currentstep=stepid ; currentstepinit=false ; } }
  void check() ;
  int prepare_Writing (Step & step)  ;
  int remove_atoms (Step &step) ;
  int data2atomes (Step &step) ;
  int gp_from_atmid (int id) ;

  // Variables
  int ngp ;  // Nombre de groupes
  vector < vector <int> > gps ;
  vector <int> gps_id ;
  vector <vector <double> > data ; //Values: flag, centroid x, y, z, orientation x, y, z.
  int type ;

  double ** Wdat ;
  double * Wint ;

private:
  vector < Vector > pts, segments ;
  bool initialized ;
  bool isflux ;
  bool symetrie[7] ;
  int currentstep ; bool currentstepinit ;
} ;
#endif
