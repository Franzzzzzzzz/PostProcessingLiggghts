#ifndef MASTERH
#include "MasterHeader.h"

#else
//bool operator== (Boites b1, Boites b2) ; 
//bool operator== (Coarse c1, Coarse c2) ; 
//#include "../newmat10/newmat.h"

/*class Boites 
{
public : 
  friend Boites operator+(Boites b1, Boites b2) ; 
  void operator=(Boites b) ; 
  friend bool operator== (Boites b1, Boites b2) ; 

 Vector x0 ;
 //double w, d, h ; 
 //double volume ; 
 double phi ; // Résolu en faisant des intersections de sphères
 double nb_atomes, nb_contacts ;
 double vx, vy, vz ; 
 double mean_radius ;
 double cfbreak, cffmax ; 
 Matrix3x3 sigmatot ; 
 Matrix3x3 sigma ; 
 Matrix3x3 sigmak ;
 //Vector pressure ;
 //Matrix sigma(3,3) ;

} ; */

//----------------------------------
typedef void (*Fenetre) (double &,  Vector &,  Vector &,  Vector &) ; 
typedef void (*FenetreO) (double &,  Vector &,  Vector &, double) ;
class Coarse 
{
public :   
  //Fonctions
  Coarse () ; 
  int coarse (Step &step) ;             // Prend en charge la construction d'un coarse
  friend Coarse operator+(Coarse c1, Coarse c2) ; // Prend en charge la somme de deux coarse
  void operator=(Coarse c) ;
  friend bool operator== (Coarse c1, Coarse c2) ; 
  int write_asVTK(ofstream &out) ;
  int write_asMatlab (ofstream &out_vx, ofstream &out_phi) ;
  int write_asMatlabTCF2D (ofstream &out_sigma) ;
  int do_coarse (Step &step, Step &stepatm, int type) ; 
  
  //int do_coarse_polydisperse(Step &step) ;
  //void calc_bornes (double x, double y, double z, int DeltaBoites, int * bornes) ;
  //int calc_deltaboites (int type) ; 
  void calc_bornes_v2 (double x, double y, double z, double r, int curtype, int *bornes) ;
  void calc_bornes_v0 (double x, double y, double z, int *bornes) ;
  
  // Variables
  int nb_boites[3], nb_boites_tot ;
  Step * cstep ; 
  bool useCF ; 

  Fenetre fenetre ; FenetreO fenetreO ; 
  Vector sigmafenetre ; Vector sigmafenetreO ; 

  // Fonctions
  int linidx (int a, int b, int c){return a*nb_boites[2]*nb_boites[1]+b*nb_boites[2]+c ;} 
  Vector taille_boite (Step &step) ; 
  void addto_sigma (int idxclosest, Matrix3x3 Sig);
  void addto_sigmak (int idxclosest, Matrix3x3 Sig);
  void compute_sigmatot (void);
  Vector position_vec (int id) ;
  
  // Variables
  double borders[3][2] ;
  vector <double> VolumeO, Volume ; 

} ;

//--------------------------------
class CoarseDump : Dump
{
public : 
  // Fonctions
  int do_coarse (Dump &dump, Dump & dumpatm, int type) ;
  int do_coarse_basic (Dump &dump, string chemin)  ;
  int write(string chemin) ;
  int mean (void) ;
  //void * thread_do_coarse (void * donnees) ;

  //Variables
  vector<Coarse> coarse ;
  
  public : 
   // Variables héritées
   /*int nbsteps ;
   ifstream dumpin ; 
   long int inscrit ;   
   vector <Step> steps ;
   vector <struct Op> prefiltre ;
   vector <struct Op> postfiltre ;
   */
   // Fonctions héritées (les virtuelles sont à redéfinir ...)

   void disp(void) {} 
   int write_asVTK (string chemin)  {return 1; }
   int write_asOneVTK (string chemin)  {return 1; }
   // TODO 
   int open (string chemin) ;
   int check_timestep (long int timestep) ;
   int free_timestep (long int inscrit) ; 
protected :
   int read (unsigned char read_type, int index) {return 1 ; } 
   int write_asVTKBase (string chemin, char type) {return 1 ; } 
   
   // Cette classe permet de parcourir les boîtes des coarses dans l'ordre x le plus rapide, y, z le moins rapide
// alors que par construction les boites sont dans l'ordre z le plus rapide ... x le moins rapide
// id. idxboit=x*sy*sz+y*sz+z.
public :   
class IterStructGrid
{
private :   
  int sx, sy, sz ;
  int ix, iy, iz ; 
  bool state ; 
  void add_one (void) {ix++ ; if (ix>=sx) {ix=0 ; iy++ ; } if (iy>=sy) {iy=0 ; iz++ ; } if (iz>=sz) {iz=0 ; state=false ;}}
  void min_one (void) {ix-- ; if (ix<0) {ix=sx-1 ; iy-- ; } if (iy<0) {iy=sy-1 ; iz-- ; } if (iz<0) {iz=sz-1 ; state=false ;}}
  
public : 
  IterStructGrid () : sx(1), sy(1), sz(1), ix(0), iy(0), iz(0) {}
  IterStructGrid (int a, int b, int c) :ix(0), iy(0), iz(0) {sx=a ; sy=b ; sz=c ; }
  IterStructGrid (int *a) :ix(0), iy(0), iz(0) {sx=a[0] ; sy=a[1] ; sz=a[2] ; }
  bool operator-= (int n) {for (int i=0 ; i<n ; i++) min_one() ; return state ;}
  bool operator-- (int i) {min_one() ; return state ; } 
  bool operator+= (int n) {for (int i=0 ; i<n ; i++) add_one() ; return state ;}
  bool operator++ (int i) {add_one() ; return state ; }
  bool operator== (bool v) {return (state==v) ; }
  bool operator=  (int i)  {ix=i/(sy*sz) ; iy=(i%(sy*sz))/sz ; ix=(i%(sy*sz)%sz) ; 
			    if (ix<sx && iy<sy && iz<sz) state=true ; else state=false ; return state ; }
  operator bool () { return state ; }
  operator int () {return (ix*sy*sz + iy*sz +iz) ; }
  int idx() {return (iz*sx*sy+ iy*sx + ix) ; }
} ;

} ; 


#endif
