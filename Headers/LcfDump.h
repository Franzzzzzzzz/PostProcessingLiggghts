#ifndef MASTERH
#include "MasterHeader.h"

#else

// Classe pour les dump par pair de contact non compressé (chain forces)
class LcfDump : public LDump {
public :
  int write_asVTK (string chemin, LucDump &dump) ; // remplace la fonction write_asVTK de la classe mère.
  int write_grainforce (string chemin, LucDump & ldump) ;
  int write_grainforcetot (string chemin) ;
  int write_radiuscontact (void) ;
  int coupletot (string chemin, LucDump & ldump) ;
  int forcetank(string chemin, LucDump & ldump) ;

private :
  // Fonctions
  int isitem (string ligne) ;
  int sparselabels(Step &step , string ligne) ;
  virtual int read (unsigned char read_type, int index);
} ;

#endif

