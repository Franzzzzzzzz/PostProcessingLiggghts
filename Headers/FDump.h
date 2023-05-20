#ifndef MASTERH
#include "MasterHeader.h"

#else
// Classe pour les pas de force (lecture de VTK
class FDump : virtual public Dump {
public :
  // Fonctions
  virtual void disp(void) ;
  int write_forcestot(string chemin) ;
  int write_coupletot(string chemin) ;
  virtual int read (unsigned char read_type, int index) ;
  int mean_stress(string chemin) ;

private :
  // Fonctions
  int identify_ligne_vtk (string ligne) ;
};

#endif
