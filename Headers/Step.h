#ifndef MASTERH
#include "MasterHeader.h"

#else
#define MASK_ALWAYS_THE_SAME 128

#define UNKNOWN 255

#define TNONE 0
#define TF 1
#define TL 2
#define TCF 3
#define TCOARSE 16

class Step {
public :
  // Variables
  char Type ;    // Variable indiquant le type de step (F ou L)
  streampos posinfile ;
  unsigned char nb_idx ;
  vector<unsigned char> idx_col ;
  vector< vector<double> > datas ;
  Step * otherstep ; // lien vers le lstep pour le cfstep
  // Fonctions
  Step() ;
  int  find_idx (int id) ;
  vector<int> find_idx (vector<int>id, bool errcatch=true) ;
  int check_idx (int id) ;
  double epsilon(double valeur) {if (valeur==0.0) return 0.0 ; else return valeur ; }
  // Aiguillages
  void disp (void) {switch(Type) { case TF : Fdisp() ; break ;
                                   case TL : Ldisp() ; break ;
                                   case TCF : LCFdisp() ; break ;
                                   default : cout << "ERR : pb d'aiguillage de step" ; break ; } }
  void write_asVTK(ofstream &out) {switch(Type) { case TF : Fwrite_asVTK(out) ; break ;
                                   case TL : Lwrite_asVTK(out) ; break ;
                                   default : cout << "ERR : pb d'aiguillage de step" ; break ; } }
  void atm_rotate (Matrix3x3 & rot, int id) {switch(Type) { case TF : Fatm_rotate(rot, id) ; break ;
                                     	 	 	 	 	 case TL : Latm_rotate(rot, id) ; break ;
                                     	 	 	 	 	 case TCF : LCFatm_rotate(rot, id) ; break ;
                                     	 	 	 	 	 default : cout << "ERR : pb d'aiguillage de step" ; break ; } }
  void write_asVTK (ofstream &out, Step &step) {LCFwrite_asVTK(out, step) ;}
  void mean_forces(double fr[]) {switch (Type) {case TF :Fmean_forces(fr) ; break ;
						case TL : Lmean_forces(fr) ; break ;
						default : cout << "ERR : pb d'aiguillage de step" ; break ; }}

  Step& operator = (const Step & step) ;

  // Variables orientées L et CF
  int timestep ;
  int nb_atomes ;  // Pour TCF, ce serait plutôt nb_entries`
  int nb_atomes_supp ; // Pour TCF, chaines periodiques qui sont renvoyées à la fin du tableau.
  bool filtered ;
  bool has_periodic_chains ;

  // Variables orientées L
  double box[3][2] ; double triclinic[3] ;
  // Fonctions orientées L
  void Ldisp (void) ;
  void Lwrite_asRESTART (ofstream &out) ;
  void Lwrite_asVTK(ofstream &out) ;
  int Latm_rotate (Matrix3x3 & rot, int id) ;
  int Lmean_forces(double fr[]) ;
  bool operator == (Step &step) ;
  void del_atm (long int atm) ;
  void swap_atm (long int atm1, long int atm2) ;
  void crush_atm (long int atm1, long int atm2) ;
  void copy_atm_end (long int nb) ;
  void del_end_atms (long int nb) ;
  void write_asDUMP (ofstream & out) ;
  int wall_force(ofstream & out, double ** meanforces, int * meangrains) ;
  double gravite (int atm) ;
  int xray_projection (int dir, int w, int h, double **img, double * box) ;
  double GetVoronoi() ;
  Multisphere * multisphere ;

  // Fonctions orientées CF
  void LCFdisp (void) ;
  int LCFwrite_asVTK(ofstream &out, Step & step) ;
  int LCFatm_rotate (Matrix3x3 & rot, int id) ;
  int LCFcouple (Step & lstep, Vector &torque, Cylindre C) ;
  int grain_force(double ** meanforces, double * meangrains, ofstream & out) ;
  int LCFforce_by_particle_v2 (Step &lstep, int type) ;
  int LCFpressure_by_particle (Step &lstep);
  int LCFforce_by_tank (Step &lstep, double * resultat) ;
  map <string, double> LCFenergy (Step &lstep) ;


  // Variables orientées F
  int nb_triangles, nb_pts ;
  // Fonctions orientées F
  int Fmean_forces(double fr[]) ;
  int Fcouple(Vector & torque, Vector c) ;
  void Fwrite_asVTK(ofstream &out) ;
  void Fdisp (void) ;
  int Fatm_rotate (Matrix3x3 & rot, int id) ;
  int buildtridata(int type) ;
  //Matrix get_tri_center (int polyidx) ;
  Vector get_tri_center (int polyidx) ;
  Vector get_tri_normal (int polyidx) ;
  double get_tri_surface (int polyidx) ;

} ;

#endif
/*class FStep {
public :
 // Fonctions
 int disp(void) ;
 FStep(void) ;
 void write(ofstream &out) ;
 int mean_forces(double fr[]) ;

 // Variables
 streampos posinfile ;
 int nb_triangles, nb_pts ;
 vector<double> pressure ;
 vector<double> shearstress ;
 vector< vector<double> > forces ;
 vector< vector<double> > normales ;
 vector< vector<double> > points ;
 vector< vector<double> > polygones ;

 int alloue ;
} ;*/
