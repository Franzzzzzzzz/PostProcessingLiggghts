#ifndef MASTERH
#include "MasterHeader.h"

#else

struct Op{
Op() : TypeDump(TNONE) {} ;   
string operation ; 
int type ; 
int idx1, idx2 ; 
double valeur ;
int TypeDump ; 
Dump * alter_dump ; 
Dump * wall_dump ;
} ;

// Modification de Filter
class Filter ; 
typedef int(Filter::*pf)(Step &, struct Op) ; 
class Filter {
public : 
  Filter() ;

  // Initialisations : 
  map <string, pf> Operations ;   // Données des opérations disponibles
  int add_operation (const char nom[], pf fonction) ; 

  // Lecture des arguments, réalisation du filtre, affichage
  vector <struct Op> parse_arg (const char arg[], int TypeDump) ; 
  int do_filter (Step &step, vector <struct Op> lst_op) ; 
  static void disp (vector <struct Op> lst_op) ; 

  static int filtrer (Step &step, const char * commande, int TypeDump) ; 
  
private : 
  struct Op parse_operation (char com[500], int length, int TypeDump) ; 

  // fonctions de traitement directement appelées par do_filter via la map Operations
  int do_nothing(Step & step, struct Op op) ; 
  int fill (Step & step, struct Op op) ; 
  int quicksort_base (Step &step, struct Op op) ; 
  int op_compare (Step &step, struct Op op) ; 
  int op_collapse (Step &step, struct Op op) ; 
  int op_collapseby (Step &step, struct Op op) ; 
  int op_operation (Step &step, struct Op op) ;
  int swappbc (Step &step, struct Op op) ;
  int zrotate (Step &step, struct Op op) ; 
  int cf_postproc (Step &step, struct Op op) ; 
  int wallforceatm (Step &step, struct Op op) ; 
  int pourcentage (Step & step, struct Op op) ; 
  int moveto_xplusyplus (Step & step, struct Op op) ;
  int copypartiallink (Step & step, struct Op op) ;
  int no_periodic_chain (Step & step, struct Op op) ;
  int lieedwardswap (Step & step, struct Op op) ;
  int moveref (Step & step, struct Op op) ;
  int multisphere (Step &step, struct Op op) ;  
  int multisphere_remove_atm (Step &step, struct Op op) ; 

  int translate_ref (Step &step, struct Op op)  {return (translate_refgen(step, op, 7)) ; } 
  int translate_refx (Step &step, struct Op op) {return (translate_refgen(step, op, 1)) ; } 
  int translate_refy (Step &step, struct Op op) {return (translate_refgen(step, op, 2)) ; }
  int translate_refz (Step &step, struct Op op) {return (translate_refgen(step, op, 4)) ; } 
  int translate_refgen (Step &step, struct Op op, int dir) ; 
  int make2D3D (Step & step, struct Op op) ; 

  //int reorder (vector<vector<double> > &tableau, int index) ;

  // Sous-fonctions du quicksort
  int quicksort (vector<vector<double> > &tableau, int premier, int dernier, int idx, int nbidx, int recurseleft) ;
  int selectsort (vector<vector<double> > &tableau, int premier, int dernier, int idx, int nbidx) ;
  int partitionner (vector<vector<double> > &tableau, int premier, int dernier, int pivot, int idx, int nbidx) ;
  void echanger (vector< vector<double> > &tableau, int idx1, int idx2, int nbidx) ;
  int heapsort(vector <vector<double> > &tableau, int deb, int fin, int idx, int nbidx) ;
  int tamiser(vector <vector<double> > &tableau, int noeud, int max, int deb, int idx, int nbidx) ;
} ; 


/*class Filter {
public : 
  // Fonctions de traitement de l'entrée de filtrage
  Filter() ;
  vector<struct Op> lst_op ; 
  enum Operations {
       null, gt, lt, gte, lte, eq, neq, gt_id, lt_id, gte_id, lte_id, eq_id, neq_id, sort, zrot, cfdumpdo
       } ; 
  int parse_arg (char arg[]) ; 
  int do_filter (Step &step) ; 
//  int do_filter (Dump &dump) ; 
  void disp (void) ; 

private : 
  int quicksort (vector<vector<double> > &tableau, int premier, int dernier, int idx, int nbidx) ;
  int partitionner (vector<vector<double> > &tableau, int premier, int dernier, int pivot, int idx, int nbidx) ;
  void echanger (vector< vector<double> > &tableau, int idx1, int idx2, int nbidx) ;
  int parse_operation (char com[500], int length) ; 
  
  int reorder (vector<vector<double> > &tableau, int index) ;
  int op_compare (Step &step, int idx, double valeur, int op) ; 
  int op_compare_id (Step &step, int idx, int idx2, int op) ; 
  int zrotate (Step &step, double angle) ; 
} ;*/


#endif
