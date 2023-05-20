#ifndef MASTERH
#include "MasterHeader.h"

#else
class Donnees {
public :
  int longueur ;
  char* datas;  
  char C ; short int SI ; int I ; long int LI; float F ; double D ; 
} ; 


class Compresser {
public : 
  int compress (Dump &dump, string nom, FILE * out) ;
  int uncompress(FILE * in) ; // Deprecated
  // Variables
  Stats champs ;
  vector <unsigned char> formats ;
  vector <bool> same ; 
  unsigned char binary_swap (unsigned char octet, unsigned int bits) ; 
  char format (int i) ; 

  enum Formats {
    CHAR=1,UCHAR,SINT,USINT,INT,UINT,LINT,ULINT,FLOAT,DOUBLE,
    CHAR_CST=17,UCHAR_CST,SINT_CST,USINT_CST,INT_CST,UINT_CST,LINT_CST,ULINT_CST,FLOAT_CST,DOUBLE_CST,
    CHAR_DENOM_2=33,UCHAR_DENOM_2,SINT_DENOM_2,USINT_DENOM_2,INT_DENOM_2,UINT_DENOM_2,LINT_DENOM_2,ULINT_DENOM_2,FLOAT_DENOM_2,DOUBLE_DENOM_2,
    CHAR_DENOM_4=49,UCHAR_DENOM_4,SINT_DENOM_4,USINT_DENOM_4,INT_DENOM_4,UINT_DENOM_4,LINT_DENOM_4,ULINT_DENOM_4,FLOAT_DENOM_4,DOUBLE_DENOM_4 ,
    CHAR_DENOM_16=65,UCHAR_DENOM_16,SINT_DENOM_16,USINT_DENOM_16,INT_DENOM_16,UINT_DENOM_16,LINT_DENOM_16,ULINT_DENOM_16,FLOAT_DENOM_16,DOUBLE_DENOM_16
  } ;  
  
private :
  int nonebynow ; 
  bool check_all (LucDump & dump) ; 
  Donnees ecrirefmt (double valeur, int fmt, FILE * out ) ; 
  int create_datas_denom (vector <unsigned char> & datasout, unsigned char fmt, int idx, vector <double> &datasin) ;
} ;
#endif
