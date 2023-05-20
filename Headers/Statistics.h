#include <float.h>

#ifndef MASTERH
#include "MasterHeader.h"

#else

#define MAX_DENOM 16

class Stats {
public : 
  // Variables
 
  //vector <unsigned char> idx_col ;
  vector <double> min, max ; 
  vector <double> mean ; // TODO
  vector <bool> same ; 
  vector <bool> canbefloat ; 
  vector <unsigned short int> denom ; 
  vector <vector <double> > denombrables ;
  vector <unsigned char> idx_col ;
  
  // Fonctions
  int compute_global (Dump &dump) ; 
  int compute_step (Dump & dump, int ts) ; 
  int compute (Step &step) ; 
  double minmaxdst(Dump &dump, int ts );
  int disp () ; 
  
  bool check_index (Dump &dump) ; 
  bool check_box   (Dump &dump) ; // TODO
} ;

#endif
