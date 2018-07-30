#ifndef MASTERH
#include "MasterHeader.h"

#else
// Valeurs < 128 : valeurs pour le dumps atomiques de Liggghts
// Valeurs >= 128 et <255 : valeurs pour les dumps de stress
// INFO : les composantes d'une même variable doivent toujours se suivre !!!
/*#define ID 1
#define TYPE 2
#define POSX 3
#define POSY 4
#define POSZ 5
#define VX 6
#define VY 7
#define VZ 8
#define RAYON 9
#define FX 10
#define FY 11
#define FZ 12
#define MASSE 13
#define OMEGAX 14
#define OMEGAY 15
#define OMEGAZ 16
#define SIGMAXX 17
#define SIGMAXY 18
#define SIGMAXZ 19
#define SIGMAYX 20
#define SIGMAYY 21
#define SIGMAYZ 22
#define SIGMAZX 23
#define SIGMAZY 24
#define SIGMAZZ 25
#define FORCEWALLX 26
#define FORCEWALLY 27
#define FORCEWALLZ 28

#define CELLDATAMIN 63
#define CFID1 64
#define CFID2 65
#define CFFORCEX 66
#define CFFORCEY 67
#define CFFORCEZ 68
#define CFMAG 69
#define CFR 70
#define CFTHETA 71
#define CFPHI 72
#define CFX 73
#define CFY 74
#define CFZ 75
#define CFPERIOD 76
#define CFID1X 77
#define CFID1Y 78
#define CFID1Z 79
#define CFID2X 80
#define CFID2Y 81
#define CFID2Z 82
#define CELLDATAMAX 127

#define PRESSURE 128
#define SHEARSTRESS 129
#define FORCEX 130
#define FORCEY 131
#define FORCEZ 132
#define NORMALEX 133
#define NORMALEY 134
#define NORMALEZ 135
#define POINTX 136
#define POINTY 137
#define POINTZ 138
#define POLY1 139
#define POLY2 140
#define POLY3 141
#define POLY4 142
#define STRESSX 143
#define STRESSY 144
#define STRESSZ 145
#define CENTREX 146
#define CENTREY 147
#define CENTREZ 148

#define COARPHI 192
#define COARATM 193
#define COARCONTACTS 194
#define COARRAD 195
#define COARSIGKXX 196
#define COARSIGKXY 197
#define COARSIGKXZ 198
#define COARSIGKYX 199
#define COARSIGKYY 200
#define COARSIGKYZ 201
#define COARSIGKZX 202
#define COARSIGKZY 203
#define COARSIGKZZ 204
#define COARSIGTOTXX 205
#define COARSIGTOTXY 206
#define COARSIGTOTXZ 207
#define COARSIGTOTYX 208
#define COARSIGTOTYY 209
#define COARSIGTOTYZ 210
#define COARSIGTOTZX 211
#define COARSIGTOTZY 212
#define COARSIGTOTZZ 213
*/
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
