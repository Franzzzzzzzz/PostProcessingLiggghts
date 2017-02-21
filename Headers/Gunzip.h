#ifndef MASTERH
#include "MasterHeader.h"

#else

#define CLUSTER 32767 // taille des blocs lus

// wrapper for ifstream to be able to decompress while reading
class gzifstream {
public :
  gzifstream() : isgz(false), nomoredata(false), cur(0), qttavail(0), curtot(0) {} ;   
  
  // Wrapped functions
  void open (const char* filename) ;
  bool operator!() {return (!in);}
  bool eof() {if (isgz) return(gzeof()) ; else return (in.eof()); } 
  void getline (string & ligne) ; 
  gzifstream& operator>> (int& val) ;
  gzifstream& operator>> (float& val) ;
  gzifstream& operator>> (double& val) ;
  unsigned long int tellg() {if (isgz) return (curtot) ; else return (in.tellg()) ;  }
  gzifstream& seekg (unsigned long int pos) ;
  gzifstream& seekg (long int off, ios_base::seekdir way) ; 
  void getline (char** s, streamsize n ) {if (isgz) gzgetline((unsigned char**)s,n) ; else in.getline(*s, n) ;  }
  void clear() {in.clear() ; }
  unsigned long int size() ; 
  //bool isgz() {return (isgz) ; }
  //bool fail() {if (isgz) printf("GZ") ; else return(in.fail());  }
  //ios_base::iostate rdstate() {if (isgz) printf("GZ") ; else return(in.rdstate()) ;  }
  
  unsigned char mybuffer [CLUSTER] ; 
  
  unsigned char getnextbyte() ;
  void getmot(char * mot) ;
  bool isgz ; 
  
private:
  bool nomoredata ; 
  char name[500] ; 
  ifstream in ;
  unsigned long int curtot, cur, qttavail ; 
  //int cur, qttavail, curtot ; 
  FILE * gzin ; 
  z_stream strm;
  unsigned char inbuffer [CLUSTER] ; 
  
  // Pointeurs de fonctions (pour les fonctions précédentes qui sont appelées souvent, accelère?) : 
  bool (*peof) () ; 
  void (*pgetlines) (string & ligne) ; 
  istream & (*pgtgti) (int & val) ; 
  istream & (*pgtgtf) (float & val) ; 
  istream & (*pgtgtd) (double & val) ; 
  void (*pgetlineC) (char* s, streamsize n ) ; 
  
  // Fonctions de dezippage
  int gzopen(const char * filename) ; 
  bool gzeof() {if (nomoredata && cur>=qttavail) return true ; else return false ; }  ; 
  void gzgetline (string & ligne) ; 
  void gzgetline (unsigned char** s, streamsize n ) ;
  int getmoredata() ;
  int unzip(unsigned long int qttwanted) ;
  void gzmove (unsigned long int pos) ; 
} ;



#endif