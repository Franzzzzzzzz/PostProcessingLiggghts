#ifndef MASTERH
#include "MasterHeader.h"

#else

/* Espaces de définition :
 * LDUMP : de 0 à 63 (éviter 63 : plutôt 0 à 62)
 * CFDUMP : de 64 à 127
 * VTK : de 128 à 191
 * COARSE : de 192 à 254
 * RESERVED : 255
 */

#define IDSMAX 256

#define MASK_ALWAYS_THE_SAME 128

#define TNONE 0
#define TF 1
#define TL 2
#define TCF 3
#define TCOARSE 16

class IDSCLASS {
private : 
  vector <string> idname ; 
  vector < vector <string> > idalias ;
  
  int find (string nom) {for (int i=0 ; i<idname.size() ; i++) {if (idname[i]==nom) return i ; }
    for (int i=0 ; i<idalias.size() ; i++) {
      for (int j=0 ; j<idalias[i].size() ; j++) 
	if (idalias[i][j]==nom) return i ; }
    return -1 ; 
  }
  void default_assign() ;
  
public :
  IDSCLASS() {default_assign() ; }
  //int new_id(char *name, int type) {return 1 ; }
  //int new_id(string name, int type) ; 
  //void default_assign() ;
  //vector <string> alias(int id) {vector <string> none ; return none ;}
  //int operator() (string name) {return 1 ;} 
  int new_id(string name, int type) ; 
  int new_id(const char *name, int type) {string nom ; nom=name ; return (new_id(nom, type)) ;} 
  int new_alias(const char *name, const char *alias) {int id=find(name) ; if (id!=-1) {idalias[id].push_back(alias) ; return (idalias[id].size()-1) ; }
									  else return -1 ; } 
  int display_all() ; 
  int operator() (string name) {return find(name) ;}  
  string operator() (int id) {return idname[id] ;}  
  vector <string> alias(int id) ; //{vector <string> retour ; for (int i=0 ; i<idalias[id].size() ; i++) retour.push_back(idalias[id][i]) ; return retour ; /*retour=idalias[id] ; return retour ; */}
  //vector <string> alias(string nom) {int u ; u=find(nom) ; if (u>=0) return idalias[u] ; else {vector <string> nullstr ; return nullstr ;}} 

  bool isTL (string nom) {if (find(nom)>=0 && find(nom)<=63) return true ; else return false ;} 
  bool isTCF (string nom) {if (find(nom)>=64 && find(nom)<=127) return true ; else return false ;} 
  bool isTF (string nom) {if (find(nom)>=128 && find(nom)<=191) return true ; else return false ;} 
  bool isTCOARSE (string nom) {if (find(nom)>=192 && find(nom)<=254) return true ; else return false ;} 
  bool exist (string nom) {if (find(nom)!=-1) return true ; else return false ; }
  bool isTL (int nom) {if (nom>=0 && nom<=63) return true ; else return false ;} 
  bool isTCF (int nom) {if ((nom)>=64 && (nom)<=127) return true ; else return false ;} 
  bool isTF (int nom) {if ((nom)>=128 && (nom)<=191) return true ; else return false ;} 
  bool isTCOARSE (int nom) {if ((nom)>=192 && (nom)<=254) return true ; else return false ;} 
  
} ;
#endif