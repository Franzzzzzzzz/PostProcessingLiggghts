#ifndef MASTERH
#include "MasterHeader.h"

#else

// Classe pour les dump atomiques de base de liggghts
class LDump : virtual public Dump {
public :
  // Variables
  // Fonctions
  int write_asDUMP(string chemin) ;
  int write_asRESTART(string chemin) ;
  bool operator == (LDump &dmp) ;
  // Fonctions masquantes
  virtual void disp(void) ;
};


// Classe pour les dump atomiques non compressés
class LucDump : public LDump
{
public :
  bool operator == (Dump &dmp) ;
  int write_forcestot(string chemin) ;
  int write_wallforce (string chemin) ;
  int write_xray (string chemin) ; 
  int write_multisphere_dumbell(string chemin) ; 
  int write_Voronoi (string chemin) ;

private :
  // Fonctions
  int isitem (string ligne) ;
  int sparselabels(Step &step , string ligne) ;
  virtual int read (unsigned char read_type, int index);
};

#include "Statistics.h"
#include "Compress.h"

// Classe pour les dump compressés
class LcpDump : public LDump, public FDump
{
public :
  LcpDump() : isfirst(true),dumpinc(NULL) {}
  virtual int free_timestep (long int inscrit) ;
  virtual void disp(void) ;
  int uncompress () ;
  Compresser cp_dat ;
  bool isfirst ;
  FILE * dumpinc ;
  char nomoriginal[500] ;
private :
  virtual int read (unsigned char read_type, int index);
  bool checkfile (void) ;
};

#endif

