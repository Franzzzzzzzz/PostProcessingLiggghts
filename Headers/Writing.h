#ifndef MASTERH
#include "MasterHeader.h"

#else

#define READY 1
#define ASKING 2
#define WAITING 3
 
enum WSignal {WHATAREYOU=1, FIRSTTS, NEXTTS, 
  ASKING1D, ASKING2D, ASKING3D, ASKINGND, ASKINGGRID, ASKINGTS, 
  LDUMP, CFDUMP, COARSEDUMP, VTKDUMP, OTHERDUMP,
  CHOOSELDUMP, CHOOSECFDUMP,
  OK, FINISH, FINI, ASKINGMULTISPHERE, UNABLE
} ; 
#define WFORMAT_NUMBER 5
enum WFormat {VTK=1, MAT, ASCII, BINARY, NONE} ; 
enum WVTKType {POLYDATA, STRUCTUREDGRID, UNSTRUCTUREDGRID} ; 

class Writing {

private : 
string WSignal_name[22] ; 
  
// Gestion de threads 
public : 
WSignal Signal ; 
pthread_t th_ecrire ;
pthread_mutex_t mutex;
pthread_cond_t sigin, sigout;
void * mfp() ; 

int *i ; double **d ; 
string *chemin ; 
int tsinfo ; 
bool init, mkdir ; 

private:
   static void * InternalThreadEntryFunc(void * This) {((Writing *)This)->start(NULL); return NULL ; }

// Variables  
public:
  //map <string, int> vardef;
  map <string, int> dim ; 
  map <string, WFormat> formats ; 
  
private:  
  int format_writing ; 
  vector < unsigned char * > varset ;   
  bool nodef, use_matlib ; 
  int char1D, chardim ; //charbysubsp ;
  vector <vector <string> > unk_var ; 
    
// Fonctions
private: 
  WFormat get_nextwritingformat() ; 
  int set_var(string nom, string var, bool val) ; 
  int set_varset(string nom, string var) ;
  int set_print(string nom, bool val) ;
  int set_def()  ; 
  int check_unregistered () ; 
  int makedir() ; 
  
  int  function_call(string function, string space, string var) ;
  void disp() ;
  void sendout (WSignal sig) ;
  
public :
  Writing() ;
  ~Writing() ;
  int startall(void) ;
  int parse (string commande) ; 
  bool get_varset(string space, string var) ;  
  int  get_nextfield (string space) ;
  int  get_nextfield (string space, string sp) ;
  string get_nextname (string var) ;
  int get_nextvardef (string var) ;
  int get_nextvardef (string var, string sp) ; 
  int get_numelem (string var) ; 
  int get_numelem (int var) ; 
  
  void disp_Signal(void) {printf("Current signal %d -> %s\n", Signal, WSignal_name[Signal].c_str()) ; fflush(stdout) ; }
  void disp_Signal(WSignal s) {printf("Given signal %d -> %s\n", s, WSignal_name[s].c_str()) ; fflush(stdout) ; }
  
  void sendin (WSignal sig) ;
  void * start (void * rien) ;
  int ldump () ;
  int lcfdump () ;
  int coarsedump () ;
  int VTKFile (WVTKType f, string nom) ; 
  
  int VTKPolyData (ofstream & out, double **datas, int n) ;
  int VTKStructuredGrid(ofstream &out, double **datas, int n, int * sizes) ;
  int VTKUnstructuredGrid (ofstream & out, double **datas, int n) ;
  int VTKCellLines (ofstream & out, double **datas, int nb) ; 
  int VTKField (ofstream & out, string name, double **datas, int dim, int n);  
  int VTKStartPointData (ofstream & out, int n) ;
  int VTKStartCellData (ofstream & out, int nb) ;
  
  int ASCIIFichier (string nom, int type) ;
  int ASCIIwrite (FILE * out, double *datas, int lig, int col) ;
  int ASCIIwrite (FILE * out, double **datas, int ndim, int n) ;
  int ASCIIwrite (ofstream & out, double **datas, int ndim, int n) ;
  int BINARYwrite (FILE * out, double *datas, int lig, int col) ;
  int BINARYwrite (FILE * out, double **datas, int ndim, int n) ;
  
  int MATFichier (string nom) ;
#ifdef MATLAB
  int MATMatrix (MATFile * out, string name, double *datas, int ndim, ...) ;
#endif
  int MATHeader_perso (FILE *out) ;
  int MATMatrix_perso (FILE *out, string name, double *datas, int ndim, ...) ; 
  
#ifdef USETIFF  
  static int TIFF_bw_write (string nom, int w, int h, double **datas) ;
#endif  
  
} ;
#endif
