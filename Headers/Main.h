#ifndef MASTERH
#include "MasterHeader.h"

#else

class Arguments {
public : 
  int id ; 
  bool set ; 
  string commande ; 
  string description ; 
  vector <int> dependances ; 
  vector <double> params ; 
  vector <string> p_names ;

  Arguments() {set=false ; id=-1 ; } 
  double operator[] (string name) ;
  int manual_set (string name, double value) ;

} ; 


class Actions {
public :
  // Fonctions
  Actions() ; 
  void aide(void) ; 

  // Fonctions de gestion de l'affichage de la progression
  void set_progress (long int *lp) ; 
  void disp_progress (void) ; 
  void clean(string chemin) ;  
  friend void * progression (void * rien) ; 
  void write_history (int todo, int argc, char * argv[]) ;
  void initialize (void) ;
  void regenerate_restart2vtk_script() ; // Special function to generate restart2vtk script (marre de la réécrire de 0 à chaque fois.
  
  // Fonctions de gestion des arguments
  vector <Arguments> arguments ; 
  // Avec surcharge ...
  int new_arg (string commande, string description, int nb_args, int nb_dep) ; 
  int new_arg (string commande, string description, int nb_args, char * noms_args[], int nb_dep) ;
  int new_arg (string commande, string description, int nb_args, int nb_dep, int id_dep[]) ;
  int new_arg (string commande, string description, int nb_args, char * noms_args[], int nb_dep, int id_dep[]) ;
  int parse_args (int argc, char * argv[]) ; 
  int parse_args_fromfile (void) ;
  int parse_alias (char com[], int argc, char * argv[]) ;
  int check_dependences (void) ; 
  void genere_bash_complete (void) ; 
  Arguments & operator[] (string name) ; 
  int copyarg (string from, string to) ;

  //Variables de timer
  double valeur, total ; 
  bool using_loop ; long int loop[3] ; 
  
  //Variable filtre global
  vector <struct Op> filtre_global ;
  
  // Chaines contenant les dumps
  int numdump ; 
  vector <string> dumpnames ; 

  // Gestion des constantes écrites en dur : 
  //  One () { Cst.add("hello", 129) ; Cst.add("Bye", 72) ;Cst.add("testons", 19) ;Cst.add("Bla", 35.625) ;}
  class Constantes 
  {
  public : 
    //Constantes () : count(0) {}
    Constantes (Actions & a) : act(a) { }
    void add(string n, double v) {char *arg[1] ; arg[0]=(char *) "valeur" ; count++ ; valeur.push_back(v); nom.push_back(n) ; info.push_back(false) ; 
      act.new_arg ("set"+n,"Régler la valeur d'une constante",1, arg,0) ;
    }
    void set_constantes () {
    for (int i=0 ; i<count ; i++)
      if (act["set"+nom[i]].set) valeur[i]=act["set"+nom[i]]["valeur"] ; 
    }
    double operator[] (string n) ; 
  private :
    int count ; 
    Actions & act ;
    vector <double> valeur ; 
    vector <string> nom ; 
    vector <bool> info ; 
  } ; 
  Constantes Cst ; 
  
  // Variable d'écriture
  Writing ecrire ; 
} ;

// Variable globale (là c'est vraiment plus simple ! )
extern Actions actions ;
extern IDSCLASS IDS ; 
void * progression (void * rien) ; 

#endif
