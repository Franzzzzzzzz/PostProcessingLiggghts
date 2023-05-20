#include "Headers/Writing.h"

Writing::Writing() : init(false), mkdir(false)
{
  int i, j, nb ;

  WSignal_name[0]="NotASignal" ; WSignal_name[1]="WHATAREYOU"; WSignal_name[2]="FIRSTTS"; WSignal_name[3]="NEXTTS";
  WSignal_name[4]="ASKING1D"; WSignal_name[5]="ASKING2D"; WSignal_name[6]="ASKING3D"; WSignal_name[7]="ASKINGND";
  WSignal_name[8]="ASKINGGRID"; WSignal_name[9]="ASKINGTS"; WSignal_name[10]="LDUMP"; WSignal_name[11]="CFDUMP";
  WSignal_name[12]="COARSEDUMP"; WSignal_name[13]="VTKDUMP"; WSignal_name[14]="OTHERDUMP"; WSignal_name[15]="CHOOSELDUMP";
  WSignal_name[16]="CHOOSECFDUMP"; WSignal_name[17]="OK"; WSignal_name[18]="FINISH"; WSignal_name[19]="FINI";
  WSignal_name[20]="ASKINGMULTISPHERE" ; WSignal_name[21]="UNABLE" ;
  if (UNABLE !=21) DISP_Err("Writing: erreur dans le programme. WSignal_name ne correspond pas aux WSignaux") ;

 //dim=vardef ;
 dim["general1D"]=1 ;
 dim["pos"]=dim["v"]=dim["f"]=dim["omega"]=2 ;
 dim["sigma"]=3 ; dim["cff"]=2 ; dim["cfforce"]=2 ;
 dim["sigmak"]=3 ;
 dim["sigmatot"]=3 ;
 dim["orient"]=3;
 dim["forcewall"]=2 ;

 formats["vtk"]=VTK ;
 formats["mat"]=MAT ;
 formats["ascii"]=ASCII ;
 formats["binary"]=BINARY ;

 varset.resize(formats.size() , NULL) ;
 char1D=ceil(IDSMAX/8.0) ;
 chardim=ceil(dim.size()/8.0) ;
 nb=char1D+chardim ;
 for (i=0 ; i<formats.size() ; i++)
 {varset[i]=(unsigned char *) malloc(nb) ;
  for (j=0 ; j<nb ; j++) varset[i][j]=0 ;
 }

 vector <string> empty ;
 unk_var.resize(formats.size(), empty) ;

 format_writing=0 ;
 set_print("vtk", true) ;
 nodef=false ; use_matlib=false ;
}

Writing::~Writing()
{
  for (int i=0 ; i<formats.size() ; i++) free(varset[i]) ;
  varset.clear() ;
}

//==============================================
void Writing::sendout (WSignal sig)
{
//pthread_mutex_lock(&mutex);
Signal=sig ; //printf("->%d %d) ", Signal, WHATAREYOU) ; fflush(stdout) ;
pthread_cond_signal(&sigout);
//pthread_mutex_unlock(&mutex);
}

void Writing::sendin (WSignal sig)
{
//pthread_mutex_lock(&mutex);
Signal=sig ; //printf("<-%d) ", Signal) ; fflush(stdout) ;
pthread_cond_signal(&sigin);
//pthread_mutex_unlock(&mutex);
}
//======================================================
int Writing::startall(void)
{pthread_create(&th_ecrire, NULL, InternalThreadEntryFunc, this); return 1 ; }

void * Writing::start (void * rien)
{
pthread_mutex_lock(&mutex);
sendout(WHATAREYOU) ;
pthread_cond_wait(&sigin, &mutex);

if (nodef==false) {set_def() ; }
if (mkdir) makedir() ;

switch(Signal)
{
  case LDUMP : ldump() ; break ;
  case CFDUMP : lcfdump() ; break ;
  case COARSEDUMP : coarsedump() ; break ;
  case VTKDUMP : break ;
  case OTHERDUMP : DISP_Warn(_("Writing : Je ne sais pas traiter ce dump.\n")) ; break ;
  default : DISP_Err(_("Writing : signal inattendu.\n")) ;
}
sendout(FINI) ; pthread_mutex_unlock(&mutex);
return NULL ;
}

int Writing::makedir()
{
 // 1 : get last part of name
 std::size_t found ;
 string dirs ; dirs=*chemin ;
 int res ;

 found = dirs.find("dump.force");
 if (found!=std::string::npos) {dirs.erase(0, 10) ; }
 found = dirs.find("dump.test");
 if (found!=std::string::npos) {dirs.erase(0, 9) ; }
 found = dirs.find("dump.forcewall");
 if (found!=std::string::npos) {dirs.erase(0, 14) ; }

 if (dirs.length()<=0) {DISP_Warn (_("Writing: impossible d'utiliser un dossier au nom vide. L'écriture se fera dans le dossier actuel.\n")) ; return 2 ; }

 // Bad chars are suppressed
 dirs.erase(std::remove(dirs.begin(), dirs.end(), ','), dirs.end());
 dirs.erase(std::remove(dirs.begin(), dirs.end(), ' '), dirs.end());
 dirs.erase(std::remove(dirs.begin(), dirs.end(), '/'), dirs.end());
 dirs.erase(std::remove(dirs.begin(), dirs.end(), '\\'), dirs.end());
 dirs.erase(std::remove(dirs.begin(), dirs.end(), '?'), dirs.end());
 dirs.erase(std::remove(dirs.begin(), dirs.end(), '%'), dirs.end());
 dirs.erase(std::remove(dirs.begin(), dirs.end(), '*'), dirs.end());
 dirs.erase(std::remove(dirs.begin(), dirs.end(), ':'), dirs.end());
 dirs.erase(std::remove(dirs.begin(), dirs.end(), '|'), dirs.end());
 dirs.erase(std::remove(dirs.begin(), dirs.end(), '<'), dirs.end());
 dirs.erase(std::remove(dirs.begin(), dirs.end(), '>'), dirs.end());
 dirs.erase(std::remove(dirs.begin(), dirs.end(), '.'), dirs.end());

 res=std::system(("mkdir "+dirs+ " >/dev/null 2>&1 ").c_str()) ;
 if (res==0) {DISP_Info(_("Writing: nouveau répertoire : ")) ; printf("%s\n", dirs.c_str()) ; }
 else if (res==1 || res==256) {DISP_Info(_("Writing: utilisation d'un répertoire existant : ")) ; printf("%s\n", dirs.c_str()) ;  }
 else {DISP_Warn(_("Writing: Une erreur s'est produite en essayant de créer un répertoire\n")) ; printf("%s\n", dirs.c_str()) ; return 3 ; }

 *chemin=(dirs+"/")+*chemin ;
 return 1 ;
}
//=========================================
int Writing::ldump ()
{
WFormat format ; string chemin2 ;
char strnum[10] ;
 while ((format=get_nextwritingformat())!=NONE)
 {
  switch(format)
  {
    case VTK :
       sendout(FIRSTTS) ; pthread_cond_wait(&sigin, &mutex);
       check_unregistered() ;
       while (Signal==OK)
       {
        sprintf(strnum, "%d", *i) ;
        try{(*chemin).erase((*chemin).rfind(".gz")) ; } catch(...){}
	chemin2=(*chemin)+"-"+strnum+".vtk" ;
        VTKFile(POLYDATA, chemin2) ;
	sendout(NEXTTS) ; pthread_cond_wait(&sigin, &mutex);
       }
      break ;
    case MAT :
      DISP_Info(_("Attention : pour l'écrire matlab de LDump, les données sont d'abord chargées en mémoire, ce qui peut prendre beaucoup de temps et risque de dépasser la capacité mémoire.\n")) ;
      try {(*chemin).erase((*chemin).rfind(".gz")) ; } catch (...) {}
      chemin2=(*chemin) + ".mat" ;
      sendout(FIRSTTS) ; pthread_cond_wait(&sigin, &mutex);
      check_unregistered() ;
      MATFichier(chemin2) ;
      break ;
    case ASCII :
      DISP_Info(_("Attention : pour l'écriture ASCII de LDump, les données sont d'abord chargées en mémoire, ce qui peut prendre beaucoup de temps et risque de dépasser la capacité mémoire.\n")) ;
      try {(*chemin).erase((*chemin).rfind(".gz")) ;} catch (...) {}
      chemin2=(*chemin) ;
      sendout(FIRSTTS) ; pthread_cond_wait(&sigin, &mutex);
      check_unregistered() ;
      ASCIIFichier(chemin2, 1) ;
      break ;
    case BINARY :
      DISP_Info(_("Attention : pour l'écriture ASCII de LDump, les données sont d'abord chargées en mémoire, ce qui peut prendre beaucoup de temps et risque de dépasser la capacité mémoire.\n")) ;
      try {(*chemin).erase((*chemin).rfind(".gz")) ; } catch (...) {}
      chemin2=(*chemin) ;
      sendout(FIRSTTS) ; pthread_cond_wait(&sigin, &mutex);
      check_unregistered() ;
      ASCIIFichier(chemin2, 2) ;
      break ;
    default : DISP_Warn(_("Writing : Unsupported format\n")) ;
  }
 }
 return 1 ;
}
//=========================================
int Writing::lcfdump ()
{
WFormat format ; string chemin2 ;
char strnum[10] ;
 while ((format=get_nextwritingformat())!=NONE)
 {
  switch(format)
  {
    case VTK :
       sendout(FIRSTTS) ; pthread_cond_wait(&sigin, &mutex);
       check_unregistered() ;
       while (Signal==OK)
       {
        sprintf(strnum, "%d", *i) ;
	try {(*chemin).erase((*chemin).rfind(".gz")) ; } catch (...) {}
	chemin2=(*chemin)+"-"+strnum+".vtk" ;
        VTKFile(UNSTRUCTUREDGRID, chemin2) ;
	sendout(NEXTTS) ; pthread_cond_wait(&sigin, &mutex);
       }
      break ;
    case MAT :
      DISP_Info(_("Attention : pour l'écrire matlab de LDump, les données sont d'abord chargées en mémoire, ce qui peut prendre beaucoup de temps et risque de dépasser la capacité mémoire.\n")) ;
      try {(*chemin).erase((*chemin).rfind(".gz")) ; } catch (...) {}
      chemin2=(*chemin) + ".mat" ;
      sendout(FIRSTTS) ; pthread_cond_wait(&sigin, &mutex);
      check_unregistered() ;
      MATFichier(chemin2) ;
      break ;
    case ASCII :
      DISP_Info(_("Attention : pour l'écriture ASCII de LDump, les données sont d'abord chargées en mémoire, ce qui peut prendre beaucoup de temps et risque de dépasser la capacité mémoire.\n")) ;
      try {(*chemin).erase((*chemin).rfind(".gz")) ; } catch (...) {}
      chemin2=(*chemin) ;
      sendout(FIRSTTS) ; pthread_cond_wait(&sigin, &mutex);
      check_unregistered() ;
      ASCIIFichier(chemin2,1) ;
      break ;
    case BINARY :
      DISP_Info(_("Attention : pour l'écriture ASCII de LDump, les données sont d'abord chargées en mémoire, ce qui peut prendre beaucoup de temps et risque de dépasser la capacité mémoire.\n")) ;
      try {(*chemin).erase((*chemin).rfind(".gz")) ; } catch (...) {}
      chemin2=(*chemin) ;
      sendout(FIRSTTS) ; pthread_cond_wait(&sigin, &mutex);
      check_unregistered() ;
      ASCIIFichier(chemin2,2) ;
      break ;
    default : DISP_Warn(_("Writing : Unsupported format\n")) ;
  }
 }
 return 1 ;
}
//=========================================
int Writing::coarsedump ()
{
WFormat format ; string chemin2 ;
char strnum[10] ;
 while ((format=get_nextwritingformat())!=NONE)
 {
  switch(format)
  {
    case VTK :
       sendout(FIRSTTS) ; pthread_cond_wait(&sigin, &mutex);
       check_unregistered() ;
       while (Signal==OK)
       {
        sprintf(strnum, "%d", *i) ;
	try {(*chemin).erase((*chemin).rfind(".gz")) ; } catch (...) {}
	chemin2=(*chemin)+"-coarse-"+strnum+".vtk" ;
        VTKFile(STRUCTUREDGRID, chemin2) ;
	sendout(NEXTTS) ; pthread_cond_wait(&sigin, &mutex);
       }
      break ;
    case MAT :
      DISP_Info(_("Attention : pour l'écrire matlab de CoarseDump, les données sont d'abord chargées en mémoire, ce qui peut prendre beaucoup de temps et risque de dépasser la capacité mémoire.\n")) ;
      try {(*chemin).erase((*chemin).rfind(".gz")) ; } catch (...) {}
      chemin2=(*chemin) + "-coarse.mat" ;
      sendout(FIRSTTS) ; pthread_cond_wait(&sigin, &mutex);
      check_unregistered() ;
      MATFichier(chemin2) ;
      break ;
    case ASCII :
      DISP_Info(_("Attention : pour l'écriture ASCII de CoarseDump, les données sont d'abord chargées en mémoire, ce qui peut prendre beaucoup de temps et risque de dépasser la capacité mémoire.\n")) ;
      try {(*chemin).erase((*chemin).rfind(".gz")) ;} catch (...) {}
      chemin2=(*chemin)+"-coarse";
      sendout(FIRSTTS) ; pthread_cond_wait(&sigin, &mutex);
      check_unregistered() ;
      ASCIIFichier(chemin2,1) ;
      break ;
    case BINARY :
      DISP_Info(_("Attention : pour l'écriture ASCII de CoarseDump, les données sont d'abord chargées en mémoire, ce qui peut prendre beaucoup de temps et risque de dépasser la capacité mémoire.\n")) ;
      try {(*chemin).erase((*chemin).rfind(".gz")) ; } catch (...) {}
      chemin2=(*chemin)+"-coarse";
      sendout(FIRSTTS) ; pthread_cond_wait(&sigin, &mutex);
      check_unregistered() ;
      ASCIIFichier(chemin2,2) ;
      break ;
    default : DISP_Warn(_("Writing : Unsupported format\n")) ;
  }
 }
 return 1 ;
}

//=============================================
int Writing::VTKFile (WVTKType f, string nom)
{
ofstream out ;
int idx[9] ;
int j, res, res2, nbpt, nbcell ;
map<string,int>::iterator it ;
bool firstcell, iscell ;

d=(double **) malloc(9*sizeof(double)) ;

sendout(CHOOSELDUMP) ; pthread_cond_wait(&sigin, &mutex);
sendout(ASKINGTS) ; pthread_cond_wait(&sigin, &mutex);
if (Signal==UNABLE) tsinfo=0 ;
else tsinfo=*i ;
out.open(nom.c_str(), ios::out) ;
idx[0]=IDS("POSX") ; idx[1]=IDS("POSY") ; idx[2]=IDS("POSZ") ;
i=idx ;
sendout(ASKING2D) ; pthread_cond_wait(&sigin, &mutex);

if (Signal==UNABLE) {DISP_Err(_("Writing: Les données de position sont indispensables")) ; }
if (f==POLYDATA) {VTKPolyData(out, d, *i) ; nbpt=*i ; }
else if (f==UNSTRUCTUREDGRID) {VTKUnstructuredGrid(out, d, *i) ; nbpt=*i ; }
else if (f==STRUCTUREDGRID)
{ nbpt=*i ;
  sendout(ASKINGGRID) ; pthread_cond_wait(&sigin, &mutex);
  VTKStructuredGrid(out, d, nbpt, i) ;
}
else DISP_Err(_("Writing : not implemented")) ;

sendout(CHOOSECFDUMP) ; pthread_cond_wait(&sigin, &mutex);
if (Signal==OK)
{
idx[0]=IDS("CFID1") ; idx[1]=IDS("CFID2") ; idx[2]=-1 ;
i=idx ;
sendout(ASKINGND) ; pthread_cond_wait(&sigin, &mutex);
if (Signal==UNABLE) {DISP_Err(_("Writing: Les données de CFID sont indispensables en CFDUMP")) ; }
VTKCellLines(out, d, *i) ; nbcell=*i ;
}

// Point datas

firstcell=false ; iscell=false ;

for (int dos=0 ; dos<2 ; dos++)
{
  if (dos==0)
  {
    sendout(CHOOSELDUMP) ; pthread_cond_wait(&sigin, &mutex);
    VTKStartPointData(out, nbpt) ;
  }
  else
  {
    sendout(CHOOSECFDUMP) ; pthread_cond_wait(&sigin, &mutex);
    if (Signal==UNABLE) { break ; }
    VTKStartCellData(out, nbcell) ;
  }
  while ( (res=get_nextfield("vtk"))  != -1)
  {
    // Construction des données
    if (res>=IDSMAX)
    {
      it=dim.begin() ; std::advance(it,res-IDSMAX);
      for (j=0 ; (res2=get_nextvardef(it->first)) != -1 ; j++ )
	idx[j]=res2 ;
    }
    else {idx[0]=res ; it=dim.find("general1D") ;}

    if (f!=STRUCTUREDGRID)
    {
     if (dos == 0 && IDS.isTCF(idx[0]))  continue ;
     if (dos == 1 && !IDS.isTCF(idx[0]))  continue ;
    }

    i=idx ;
    switch(it->second)
    {
      case 1 : sendout(ASKING1D) ; break ;
      case 2 : sendout(ASKING2D) ; break ;
      case 3 : sendout(ASKING3D) ; break ;
      default : DISP_Warn(_("Writing: Dimension inconnue")) ;
    }
    pthread_cond_wait(&sigin, &mutex);
    if (Signal==UNABLE) continue ;

    if (res>=IDSMAX) VTKField(out, it->first, d,it->second, *i) ;
    else VTKField (out, (IDS.alias(res))[0], d, 1, *i) ;
  }
}
out.close() ;
free(d) ;
return 1 ;
}


//--------------------------------------------
int Writing::parse (string commande)
{
int i, j ;
int ptoc ;
bool newoc, newvar ;
string tok, space, var, function ;
int type ;
// On enlève tous les blancs :
commande.erase(remove_if( commande.begin(), commande.end(), ::isspace ), commande.end() );
commande.push_back(';') ;  // Just to be sure ...

newoc=false ; newvar=false ;

ptoc=0 ;
for (i=1 ; i<commande.length() ; i++)
{
 if (commande[i]=='[' || commande[i]==']' ||
     commande[i]=='{' || commande[i]=='}' ||
     commande[i]==',' || commande[i]==';') // The end of a commande, the beginning of another
 {
   if (i-ptoc-1>0) tok=commande.substr(ptoc+1, i-ptoc-1) ;
   if      (commande[i]=='}') {type=1 ; function=tok ; } // Functions
   else if (commande[i]==']') {type=2 ; space=tok ;} // spaces
   else
   {
     var=tok ; // variables
     if (type==1)
       function_call (function, space, var) ;
     else if (type==2)
       set_varset(space,var) ;
     else
       DISP_Warn(_("Writing : pas d'espace ou de fonction definie")) ;
   }
   ptoc=i ;
 }
}
//disp() ;
return 0 ;
}

int Writing::set_print(string nom, bool val)
{
 map<string,WFormat>::iterator it2 ;
 it2=formats.find(nom) ;

 if (it2==formats.end()) {DISP_Warn(_("Writing : chaine de print non trouvée")) ; printf("%s/", nom.c_str()) ; }

 if (val==true) format_writing=(format_writing | (1<<(std::distance(formats.begin(), it2)))) ;
 else format_writing=(format_writing & (~(1<<(std::distance(formats.begin(), it2))))) ;
 return 1 ;
}


int Writing::set_varset (string nom, string var)
{
bool set ; string nom2;
 if (nom[0]=='-')      {nom2=nom.substr(1, string::npos) ; set_print(nom2, false) ;}
 else if (nom[0]=='+') {nom2=nom.substr(1, string::npos) ; set_print(nom2, true) ; }
 else nom2=nom ;

 if (!var.compare("nodef")) {nodef=true ; return 1 ; }
 if (!var.compare("def")) {nodef=false ; return 1; }
 if (var[0]=='-') set_var (nom2, var.substr(1, string::npos), false) ;
 else set_var(nom2, var, true) ;
 return 1 ;
}

int Writing::set_var(string nom, string var, bool val)
{
 int form, v ;
 map<string,int>::iterator it ;
 form=std::distance(formats.begin(), formats.find(nom)) ;

 // 1 : on cherche dans toutes les variables de base (1D)
 v=IDS(var) ;
 if (v>=0) // trouvé
 {
   if (val) {varset[form][v/8]=(varset[form][v/8] | (1<<v%8)) ; }
   else     varset[form][v/8]=(varset[form][v/8] & (~(1<<v%8))) ;

   return v ;
 }

 // 2: on cherche dans les variables dim (2D et 3D).
 it=dim.find(var) ;
 if (it!=dim.end()) // trouvé
 {
  v=std::distance(dim.begin(),it) ;
  v+=IDSMAX ;
  if (val) varset[form][v/8]=(varset[form][v/8] | (1<<v%8)) ;
  else     varset[form][v/8]=(varset[form][v/8] & (~(1<<v%8))) ;
  return v ;
 }

 // 3: rien trouvé, on enregistre la variable, elle existera peut-être plus tard.
 if (val) {unk_var[form].push_back(var) ; return -1 ; }
 return -2 ;
}

int Writing::set_def()
{
 switch (Signal) {
   case LDUMP : /*set_print("vtk", true) ;*/ set_varset("mat", "pos") ; break ;
   case CFDUMP: /*set_print("vtk", true) ;*/ break ;
   case COARSEDUMP : /*parse("[+vtk]phi,v,sigma[mat]pos,phi,v,sigma") ; */ break ;
   case VTKDUMP : break ;
   case OTHERDUMP : DISP_Info(_("Writing : pas de mise à défaut pour OTHERDUMP\n")) ; break ;
   default : DISP_Warn(_("Writing : signal incohérent pour setdef\n")) ;
 }
 return 1 ;
}

//---------------------------------------
int Writing::function_call(string function, string space, string var)
{
int val ;
string var2 ;
if (function=="expand")
{
  do
  {
    var2=get_nextname(var) ;
    if (var2!="false") set_varset(space,var2) ;
  } while (var2!= "false") ;
}
else if (function=="matlib")
{
 if (space=="mat")
#ifdef MATLAB
   use_matlib=true ;
#else
   {use_matlib=false ; DISP_Warn(_("Impossible d'utiliser la librairie Matlab, PostProc n'a pas été compilé avec le support Matlab.\n")) ;}
#endif
}
else if (function=="nomatlib")
{
 if (space=="mat") use_matlib=false ;
}
else
  DISP_Warn(_("Unknown writing function")) ;

return 1 ;
}
/*
if (!function.compare("print"))
 {
  if (!var.compare("all")) {typewrite=0 ; typewrite=(~typewrite) ; }
  else if (!var.compare("all:vtk")) {parse("{print}ldump:vtk,cfdump:vtk,coarse:vtk") ; }
  else if (!var.compare("all:mat")) {parse("{noprint}ldump:mat,cfdump:mat,coarse:mat") ; }
  else
  {
    try {
      val=types.at(var) ;
      typewrite=(typewrite|(1<<val)) ;
    }
    catch (...) { DISP_Warn("Writing: Variable de fonction inconnue \n") ; }
  }
 }

 else if (!function.compare("noprint"))
 {
  if (!var.compare("all")) {typewrite=0 ; }
  else if (!var.compare("all:vtk")) {parse("{noprint}ldump:vtk,cfdump:vtk,coarse:vtk") ; }
  else if (!var.compare("all:mat")) {parse("{noprint}ldump:mat,cfdump:mat,coarse:mat") ; }
  else
  {
    try {
      val=types.at(var) ;
      typewrite=(typewrite&(~(1<<val))) ;
    }
    catch (...) { DISP_Warn("Writing: Variable de fonction inconnue \n") ; }
  }
 }

 else if (!function.compare("default"))
 {
 parse("{print} all:vtk [ldump:vtk] id [ldump:mat] id, pos [cfdump:vtk] f [cfdump:mat] f [coarse:vtk] v, sigma, phi [coarse:mat] pos, v, sigma, phi") ;
 }

 else
   DISP_Warn("Writing: Fonction inconnue") ;
 */


void Writing::disp()
{
  map <string, WFormat>::iterator itfrmt=formats.begin() ;
  printf(_("-----> Selected formats 0x%X\n"), format_writing) ;
  for (vector<unsigned char *>::iterator it=varset.begin() ; it<varset.end() ; it++)
  { printf("%-8s>", itfrmt->first.c_str()) ; itfrmt++ ;
    for (int i=0 ; i<char1D+chardim ; i++)
      printf("%02x|", (*it)[i]) ;
  printf("\n") ;
  }
}

int Writing::get_nextfield (string format)
{
 static int numero=-1 ;
 int form ;
 form=std::distance(formats.begin(), formats.find(format)) ;

 if (numero==-1) {check_unregistered() ; numero=0 ;}
 while (numero<(char1D+chardim)*8)
 {
  if ((varset[form][numero/8] & (1<<(numero%8)))>0)
   {
    numero++ ;
    return (numero-1) ;
   }
  numero++ ;
 }
 numero=-1 ;
 return numero ;
}

WFormat Writing::get_nextwritingformat()
{
 static int i=-1 ;
 map<string,WFormat>::iterator it ;
 if (i==-1) i=0 ;

 it=formats.begin() ;
 for (i=i ; i<formats.size() ; i++)
 {if ((format_writing & (1<<i))>0) {std::advance(it, i) ; i++ ; return it->second ; } }
 i=-1 ;
 return NONE ;
}

/*int Writing::get_nextfield (string space, string sp) // Using an alternative subspace than the linked one
{
 static int numero=-1 ;
 if (numero==-1) numero=0 ;
 while (numero<charbysubsp*8)
 {
  if ((varset[types[space]][subspaces[sp]*charbysubsp+numero/8] & (1<<(numero%8)))>0)
   {
    //printf("[%d %d %d]", types[space], lntypessubsp[space]*charbysubsp+numero/8, numero%8) ;
    numero++ ;
    return (numero-1) ;
   }
  numero++ ;
 }
 numero=-1 ;
 return numero ;
}*/


int Writing::get_numelem (string var)
{
 map<string,int>::iterator it ;

 it=dim.find(var) ;
 if (it==dim.end()) return 1 ;
 else if (it->second==2) return 3 ;
 else if (it->second==3) return 9 ;
 else {DISP_Warn (_("Writing : dimension inconnue")) ; return 0 ;}
}
int Writing::get_numelem (int var)
{
 if (var==1) return 1 ;
 else if (var==2) return 3 ;
 else if (var==3) return 9 ;
 else {DISP_Warn (_("Writing : dimension inconnue")) ; return 0 ;}
}

string Writing::get_nextname (string var)
{
 static int num=-1 ;
 map<string,int>::iterator it ;
 it=dim.find(var) ;
 if (num==-1) num=0 ;
 if (it==dim.end())
 {
   if (num==0) {num++ ; return (var) ; }
   else {num=-1 ; return "false" ; }
 }
 else if (it->second==2)
 {
   switch(num)
   {
     case 0 : num++ ; return  (var+"x") ;
     case 1 : num++ ; return  (var+"y") ;
     case 2 : num++ ; return  (var+"z") ;
     case 3 : num=-1 ;return  ("false") ;
     default : DISP_Warn(_("Writing : should not happen\n")) ; break ;
  }
 }
 else if (it->second==3)
 {
   switch(num)
   {
     case 0 : num++ ; return  (var+"xx") ;
     case 1 : num++ ; return  (var+"xy") ;
     case 2 : num++ ; return  (var+"xz") ;
     case 3 : num++ ; return  (var+"yx") ;
     case 4 : num++ ; return  (var+"yy") ;
     case 5 : num++ ; return  (var+"yz") ;
     case 6 : num++ ; return  (var+"zx") ;
     case 7 : num++ ; return  (var+"zy") ;
     case 8 : num++ ; return  (var+"zz") ;
     case 9 : num=-1 ;return  ("false") ;
     default : DISP_Warn(_("Writing : should not happen\n")) ; break ;
  }
 }
 else
   DISP_Warn(_("Writing : dimension inconnue \n")) ;
num=-1 ;
return "false" ;
}

int Writing::get_nextvardef (string var)
{
return(get_nextvardef(var, "atm")) ;
}

int Writing::get_nextvardef (string var, string sp)
{
  string a=get_nextname (var) ;
  if (!a.compare("false")) return -1 ;
  if (!sp.compare("cf")) a="cf"+a ;
  return (IDS(a)) ;
}

int Writing::check_unregistered ()
{
 int i, j, v, cpt=0 ;
 map<string,WFormat>::iterator it ;
for (i=0 ; i<unk_var.size() ; i++)
{
 it=formats.begin() ; std::advance(it, i) ;
 for (j=0 ; j<unk_var[i].size() ; j++)
 {
  v=set_var(it->first, unk_var[i][j], true) ;
  if (v==-1) unk_var[i].pop_back() ;
  else cpt++ ;
 }
}
return cpt ;
}

//=========================================================
//-------------- ACTUAL VTK WRITING FUNCTIONS -----------------
//=========================================================
int Writing::VTKPolyData (ofstream & out, double **datas, int n)
{
  int ii, j ; double dst ;
  out << "# vtk DataFile Version 2.0\n" << "TS=" << tsinfo << ". Made by PostProcessing (CC-BY-NC)\n" << "ASCII \n" ;
  out << "DATASET POLYDATA\n" ;
  out << "POINTS " << n << " double\n" ;

  for (ii=0 ; ii<n ; ii++)
  {
   for (j=0 ; j<3 ; j++)
     out << datas[j][ii] << " " ;
   out << "\n" ;
  }

  out << "VERTICES " <<n << " "<< n*2 << "\n" ;
  for (ii=0 ; ii<n ; ii++)
  {out << "1 " << ii <<"\n" ;}

  if (actions["multisphere"].set || actions["multisphereflux"].set)
  {
    double ** dtmp ; int * itmp ;
    dtmp=d ; itmp=i ;
    sendout(ASKINGMULTISPHERE) ; pthread_cond_wait(&sigin, &mutex);
    // Recast le ** double en *<vector <vector <int>>
    Multisphere * ptr ;
    int ngp, ngp_real=0, nnumbers=0 ; int deb, end, delta ;
    ptr=(Multisphere *)(d) ;
    ngp=*i ;
    if (actions["multisphere"].set) {deb=1 ; end=ngp+1 ; delta=1 ; }
    else {deb=0 ; end=ngp ; delta=0 ; }
    for (ii=deb ; ii<end ; ii++) if (ptr->gps[ii][0]>0 && ptr->data[0][ii]==GP_OK) {ngp_real++ ; nnumbers+=ptr->gps[ii][0]+1 ;}
    out << "LINES " << ngp_real << " "<< nnumbers << "\n" ;
    for (ii=deb ; ii<end ; ii++)
    {
      if (ptr->gps[ii][0]>0 && ptr->data[0][ii]==GP_OK)
      {
        out << ptr->gps[ii][0] << " " ;
        for (j=1 ; j<ptr->gps[ii][0]+1 ; j++)
            out << ptr->gps[ii][j]-delta << " " ;
        out << "\n" ;
      }
    }
    i=itmp ; d=dtmp ;
  }

  return 1 ;
}

int Writing::VTKUnstructuredGrid (ofstream & out, double **datas, int n)
{
  int i, j ;
  out << "# vtk DataFile Version 2.0\n" << "TS=" << tsinfo << ". Made by PostProcessing (CC-BY-NC)\n" << "ASCII \n" ;
  out << "DATASET UNSTRUCTURED_GRID\n" ;
  out << "POINTS " << n << " double\n" ;

  for (i=0 ; i<n ; i++)
  {
   for (j=0 ; j<3 ; j++)
     out << datas[j][i] << " " ;
   out << "\n" ;
  }
  return 1 ;
}

int Writing::VTKStructuredGrid(ofstream &out, double **datas, int n, int * sizes)
{
int i,j,k ;
out << "# vtk DataFile Version 2.0\n" << "TS=" << tsinfo << ". Made by PostProcessing (CC-BY-NC)\nASCII\n" << "DATASET STRUCTURED_GRID\n" ;
out << "DIMENSIONS " <<sizes[0] << " " << sizes[1] << " " << sizes[2] << "\n" ;
out << "POINTS " << n << " double\n" ;
for (i=0 ; i<n; i++)
  {
    for (j=0 ; j<3 ; j++)
      out << datas[j][i] << " " ;
    out << "\n" ;
  }
return 1 ;
}

int Writing::VTKCellLines (ofstream & out, double **datas, int nb)
{
  int i, j ;
  out << "CELLS " << nb << " " << (nb*3) << "\n" ;
  for (i=0 ; i<nb ; i++)
  {out << "2 " << (int)datas[0][i]-1 << " " << (int)datas[1][i]-1 << "\n" ; }
  out << "CELL_TYPES " << nb << "\n" ;
  for (i=0 ; i<nb ; i++) out << "3\n" ;
  out << "\n" ;
return 1 ;
}

int Writing::VTKField (ofstream & out, string name, double **datas, int dim, int n)
{
 int i, j ;
 int nbcomp ;

 nbcomp=get_numelem(dim);

 if (dim==1)  out << "SCALARS " ;
 else if (dim==2) out << "VECTORS " ;
 else if (dim==3) out << "TENSORS " ;
 else DISP_Warn(_("Writing: dimension annormale pour les données à écrire VTK")) ;

 out << name <<" double " ;
 if (dim==1)  out << "1 \nLOOKUP_TABLE default\n" ;

 for (i=0 ; i<n ; i++)
 {
  for (j=0 ; j<nbcomp ; j++)
  {
   //if (datas[j][i]<=std::numeric_limits<double>::min()) datas[j][i]=0. ;
   out<< datas[j][i] << " ";
  }
  out << "\n" ;
 }
 out << "\n" ;
 return 1 ;
}

int Writing::VTKStartPointData (ofstream & out, int n)
{
out << "POINT_DATA " << n << "\n" ;
return 1 ;
}

int Writing::VTKStartCellData (ofstream & out, int nb)
{
out << "CELL_DATA " << nb << "\n" ;
return 1 ;
}

//=========================================================
//-------------- ACTUAL WRITING MATLAB FUNCTIONS ----------
//=========================================================
int Writing::MATFichier (string nom)
{
  double ** datas, *boite ;
  int *idx, *dims, j, k, l, cpt, cptvar, ts, natm, ncell, res, res2, cfiddeb ;
  int natm_orig, ncell_orig ;
  int **nbelem ; char ** noms ;

  boite=NULL ;
  idx=(int *)malloc(1*sizeof(int)) ; dims=(int *)malloc(1*sizeof(int)) ;
  nbelem=(int **)malloc(1*sizeof(int*)) ; noms =(char **) malloc(1*sizeof(char *)) ;
  if (nbelem==NULL || idx==NULL || dims==NULL || noms==NULL) DISP_Err(_("Writing : l'allocation a échouée")) ;
  cpt=0 ; cptvar=0 ;
  while ( (res=get_nextfield("mat")) != -1 )
  {
    cptvar++ ;
    nbelem=(int**)realloc(nbelem, cptvar*sizeof(int*)) ;
    noms=(char **)realloc(noms, cptvar*sizeof(char *)) ;
    dims=(int*)realloc(dims, cptvar*sizeof(int)) ;

    if (res>=IDSMAX)
    {
     map<string,int>::iterator it ;
     it=dim.begin() ; std::advance(it,res-IDSMAX) ;
     noms[cptvar-1]=(char*) malloc (sizeof(char)*(it->first.length()+1)) ;
     strcpy(noms[cptvar-1], it->first.c_str() ) ;
     dims[cptvar-1]=get_numelem(it->second) ;
     while ((res2=get_nextvardef(it->first)) != -1)
      { cpt++ ;
      idx=(int*)realloc(idx, cpt*sizeof(int)) ;
      idx[cpt-1]=res2 ;
      if (IDS.isTCF(res2)) nbelem[cptvar-1]=&ncell ;
      else nbelem[cptvar-1]=&natm ;
      }
    }
    else
    {
     noms[cptvar-1]=(char*) malloc (sizeof(char)*((IDS.alias(res))[0].length()+1)) ;
     strcpy(noms[cptvar-1], (IDS.alias(res))[0].c_str() ) ;
     dims[cptvar-1]=1 ;
     cpt++ ;
     idx=(int*)realloc(idx, cpt*sizeof(int)) ;
     idx[cpt-1]=res ;
     if (IDS.isTCF(res)) nbelem[cptvar-1]=&ncell ;
     else nbelem[cptvar-1]=&natm ;
    }
  }
  idx=(int*)realloc(idx, (cpt+1)*sizeof(int)) ; idx[cpt]=-1 ;


  // reordering par ordre croissant (les différentes composantes d'une variable de dimension 2 ou 3 se suivent toujours !!!
  int deplacements, tmp ; int * tmpptr ; char * tmpstr ;
  do
  {
    deplacements=0 ; k=1 ;
    for (j=dims[0] ; j<cpt ; j+=dims[k-1])
    {
     if (idx[j]<idx[j-1]) // il faut permuter, et ça c'est pas drôle
      { deplacements++ ;
	for (l=0 ; l<dims[k] ; l++)
	{
	  for (int m=0 ; m<dims[k-1] ; m++)
	  {
	   tmp=idx[j+l-m] ; idx[j+l-m]=idx[j+l-m-1] ; idx[j+l-m-1]=tmp ;
	  }
	}
	j-=dims[k-1] ;
	tmp=dims[k] ; dims[k]=dims[k-1] ; dims[k-1]=tmp ;
	j+=dims[k-1] ;
	tmpptr=nbelem[k] ; nbelem[k]=nbelem[k-1] ; nbelem[k-1]=tmpptr ;
	tmpstr=noms[k] ; noms[k]=noms[k-1] ; noms[k-1]=tmpstr ;
      }
      k++ ;
    }
  } while (deplacements>0) ;
  cfiddeb=-1 ;
  for (j=0 ; j<cpt ; j++) {if (IDS.isTCF(idx[j])) {cfiddeb=j ; break ; }}

  datas=(double **)malloc(cptvar*sizeof(double *)) ;
  d=(double **) malloc(cpt*sizeof(double*)) ;
  if (datas==NULL || d==NULL) DISP_Err(_("Writing:allocation impossible\n")) ;

  sendout(FIRSTTS) ; pthread_cond_wait(&sigin, &mutex);
  sendout(ASKINGGRID) ; pthread_cond_wait(&sigin, &mutex);
  if (Signal != UNABLE)
   {
#ifdef MATLAB
     if (use_matlib)
	boite=(double *)mxCalloc(3, sizeof(double)) ;
     else
        boite=(double *)calloc(3, sizeof(double)) ;
#else
     boite=(double *)calloc(3, sizeof(double)) ;
#endif
   boite[0]=(double)(*i) ; boite[1]=(double)(*(i+1)) ; boite[2]=(double)(*(i+2)) ;
   }

  ts=0 ;
  sendout(FIRSTTS) ; pthread_cond_wait(&sigin, &mutex);
  while (Signal==OK)
  {
    int tmp, val ;
    ncell=-1 ; natm=-1 ;
    for (int dos=0 ; dos<2 ; dos++)
    {
    if (cfiddeb==-1) {if(dos==1) continue ; sendout(CHOOSELDUMP) ; pthread_cond_wait(&sigin, &mutex); i=idx ;}
    else if (cfiddeb== 0) {if(dos==1) continue ; sendout(CHOOSECFDUMP) ; pthread_cond_wait(&sigin, &mutex); i=idx ;}
    else {
        if (dos==0) {sendout(CHOOSECFDUMP) ; pthread_cond_wait(&sigin, &mutex); i=idx+cfiddeb ; }
        else {for (k=cpt-cfiddeb-1 ; k>=0 ; k--) d[cfiddeb+k]=d[k] ;
	  ncell=val ;
	  sendout(CHOOSELDUMP) ; pthread_cond_wait(&sigin, &mutex); i=idx ; tmp=idx[cfiddeb] ; idx[cfiddeb]=-1 ; }
      }
    sendout(ASKINGND) ; pthread_cond_wait(&sigin, &mutex);
    val=*i ;
    }
    if (cfiddeb>0) {idx[cfiddeb]=tmp ; natm=val ; }
    if (cfiddeb==-1) {natm=val ;}
    if (cfiddeb==0) {ncell=val ;}
    if (ts==0)
    {
      ncell_orig=ncell ; natm_orig=natm ;
      if (ncell_orig==0 || natm_orig==0) DISP_Err(_("Writing : des données ont été demandées mais sont vides pour ce ts.")) ;
    }
    else
    {
      if (ncell!=ncell_orig || natm!=natm_orig)
	DISP_Err(_("Writing : pas le même nombre d'atomes ou de cellules entre ts. Impossible d'écrire Matlab")) ;
    }

    cpt=0 ;
    for (j=0 ; j<cptvar ; j++)
    {
#ifdef MATLAB
      if (use_matlib)
      {
        if (ts==0) {datas[j]=(double *)mxCalloc((*nbelem[j])*(ts+1)*dims[j], sizeof(double)) ; }
        else 	  datas[j]=(double *)mxRealloc(datas[j], (*nbelem[j])*(ts+1)*dims[j]*sizeof(double)) ;
      }
      else
      {
	if (ts==0) {datas[j]=(double *)calloc((*nbelem[j])*(ts+1)*dims[j], sizeof(double)) ; }
        else 	  datas[j]=(double *)realloc(datas[j], (*nbelem[j])*(ts+1)*dims[j]*sizeof(double)) ;
      }
#else
      if (ts==0) {datas[j]=(double *)calloc((*nbelem[j])*(ts+1)*dims[j], sizeof(double)) ; }
      else 	  datas[j]=(double *)realloc(datas[j], (*nbelem[j])*(ts+1)*dims[j]*sizeof(double)) ;
#endif
      if (datas[j]==NULL) {DISP_Err(_("Writing: pas assez de mémoire à allouer")) ; fflush(stdout) ; }

      for (k=0 ; k<dims[j] ; k++)
      {
	for (l=0 ; l<(*nbelem[j]) ; l++)
	{datas[j][ts* (*nbelem[j]) *dims[j]+k*(*nbelem[j])+l]=d[cpt][l] ;}
        cpt++ ;
      }
    }
    ts++ ;
    sendout(NEXTTS) ; pthread_cond_wait(&sigin, &mutex);
  }

#ifdef MATLAB
 if (use_matlib)
 {
   MATFile * mout ;
   mout=matOpen(nom.c_str(), "w");
   if (boite!=NULL) {MATMatrix (mout, "Grille", boite, 2, 3,1) ; mxFree(boite) ;}
   for (j=0 ; j<cptvar ; j++)
   {
     MATMatrix (mout, noms[j], datas[j], 3, (*nbelem[j]), dims[j], ts) ;
   }
   matClose(mout) ;
 }
 else {
#endif
 FILE *out ;
 out=fopen(nom.c_str(), "wb") ;
 MATHeader_perso(out) ;
 if (boite!=NULL) {MATMatrix_perso (out, "Grille", boite, 2, 3,1) ; free(boite) ;}
 for (j=0 ; j<cptvar ; j++)
 {
   MATMatrix_perso (out, noms[j], datas[j], 3, (*nbelem[j]), dims[j], ts) ;
 }
 fclose(out) ;

#ifdef MATLAB
      }
#endif
 return 1 ;
}
//--------------------------------------------
#ifdef MATLAB
int Writing::MATMatrix (MATFile * out, string name, double *datas, int ndim, ...)
{
mxArray *pm ;
mwSize n, *dims ;
va_list dimsizes;
int i, tmp ;

n=ndim ;
dims=(mwSize *) malloc(ndim*sizeof(mwSize)) ;
va_start(dimsizes, ndim);
for (i=0 ; i<ndim ; i++)
{
  tmp=va_arg(dimsizes,int);
  dims[i]=tmp ;
}
va_end(dimsizes) ;

pm=mxCreateNumericArray(n, dims, mxDOUBLE_CLASS, mxREAL);

mxSetData (pm, datas) ;
matPutVariable(out, name.c_str() , pm);

return 1 ;
//out=matOpen("Test.mat", "w");
//const unsigned long int one=1 ;
//pm=mxCreateStructArray(1, &one, 3, field_names);

//matClose(out) ;
}
#endif

#define miINT8 1
#define miUINT8 2
#define miINT16 3
#define miUINT16 4
#define miINT32 5
#define miUINT32 6
#define miDOUBLE 9
#define miMATRIX 14
#define mxDOUBLE_CLASS 6
int Writing::MATHeader_perso (FILE *out)
{
  unsigned int i ; char description[116] ;
  time_t t ;
  // Constants ...
  unsigned int zero=0 ; short int version=0x0100 ;
  unsigned char endianstr[2] ; endianstr[0]='I' ; endianstr[1]='M' ;

  for (i=0 ; i<116 ; i++) description[i]=0 ;
  sprintf(description, "MATLAB 5.0 MAT-file, Platform: xxx, Created on: %s by POSTPROCESSING", ctime(&t)) ;
  fwrite(description, 1, 116, out) ;
  fwrite(&zero, 4, 1, out) ; fwrite(&zero, 4, 1, out) ;
  fwrite(&version, 2, 1, out) ; fwrite(&endianstr, 1, 2, out) ;

  return 1 ;
}
//-------------------------------------------------------------
int Writing::MATMatrix_perso (FILE *out, string name, double *datas, int ndim, ...)
{
unsigned int *dims, tmp, i, tot ;
va_list dimsizes;

dims=(unsigned int *) malloc(ndim*sizeof(unsigned int)) ;
va_start(dimsizes, ndim);
for (i=0 ; i<ndim ; i++)
{
  tmp=va_arg(dimsizes,int);
  dims[i]=tmp ;
}
va_end(dimsizes) ;

tot=1 ;
for (i=0 ; i<ndim ; i++) tot*=dims[i] ;

unsigned int headdat[20] ; int n=0 ;
unsigned char czero=0 ;
headdat[n++]=miMATRIX ; headdat[n++]=4*(6+ndim+ndim%2+2+ceil(name.length()/8.)*2+2)+tot*8 ; // total size all tableau
// Flag subelement
headdat[n++]=miUINT32 ; headdat[n++]=8 ; headdat[n++]=mxDOUBLE_CLASS ; headdat[n++]=0 ;
// Dimensions subelement
headdat[n++]=miINT32 ; headdat[n++]=ndim*4 ; // Number of dimension
for (i=0 ; i<ndim ; i++) headdat[n++]=dims[i] ;
if (ndim%2==1) headdat[n++]=0 ; // padding
// Array name subelement
headdat[n++]=miINT8 ; headdat[n++]=name.length() ;
fwrite(headdat, 4, n, out) ;
fwrite(name.c_str(), 1, name.length(),out) ;
for (i=0 ; i<(8-name.length()%8)%8 ; i++) fwrite(&czero, 1, 1,out) ;
// Datas
headdat[0]=miDOUBLE ; headdat[1]=tot*8 ; fwrite(headdat, 4, 2,out) ;
fwrite(datas, 8, tot, out) ; // Pas de padding car les doubles tombent sur les frontières de 64 bits ...
free(dims) ;
return 1 ;
}

//=========================================================
//-------------- ACTUAL WRITING ASCII FUNCTION ----------
//=========================================================
int Writing::ASCIIFichier (string nom, int type)
{
  FILE * out ; char fmtstr[20] ;
  double ** datas, *boite ;
  int *idx, *dims, j, k, l, cpt, cptvar, ts, natm, ncell, res, res2, cfiddeb ;
  int natm_orig, ncell_orig ;
  int **nbelem ; char ** noms ;
  string chemin ;

  boite=NULL ;
  idx=(int *)malloc(1*sizeof(int)) ; dims=(int *)malloc(1*sizeof(int)) ;
  nbelem=(int **)malloc(1*sizeof(int)) ; noms =(char **) malloc(1*sizeof(char *)) ;
  if (nbelem==NULL || idx==NULL || dims==NULL || noms==NULL) DISP_Err(_("Writing : l'allocation a échouée")) ;
  cpt=0 ; cptvar=0 ;
  if (type==1) strcpy(fmtstr, "ascii") ;
  else strcpy (fmtstr, "binary") ;
  while ( (res=get_nextfield(fmtstr)) != -1 )
  {
    cptvar++ ;
    nbelem=(int**)realloc(nbelem, cptvar*sizeof(int*)) ;
    noms=(char **)realloc(noms, cptvar*sizeof(char *)) ;
    dims=(int*)realloc(dims, cptvar*sizeof(int)) ;
    if (res>=IDSMAX)
    {
     DISP_Info(_("Pour l'écriture ASCII, les variables de dimensions 2 et 3 ne peuvent être écrites que moyennées en temps / avec 1 seul ts / avec 1 seul atome.\n")) ;
     map<string,int>::iterator it ;
     it=dim.begin() ; std::advance(it,res-IDSMAX) ;
     noms[cptvar-1]=(char*) malloc (sizeof(char)*(it->first.length()+1)) ;
     strcpy(noms[cptvar-1], it->first.c_str() ) ;
     dims[cptvar-1]=get_numelem(it->second) ;
     while ((res2=get_nextvardef(it->first)) != -1)
      { cpt++ ;
      idx=(int*)realloc(idx, cpt*sizeof(int)) ;
      idx[cpt-1]=res2 ;
      if (IDS.isTCF(res2)) nbelem[cptvar-1]=&ncell ;
      else nbelem[cptvar-1]=&natm ;
      }
    }
    else
    {
     noms[cptvar-1]=(char*) malloc (sizeof(char)*((IDS.alias(res))[0].length()+1)) ;
     strcpy(noms[cptvar-1], (IDS.alias(res))[0].c_str() ) ;
     dims[cptvar-1]=1 ;
     cpt++ ;
     idx=(int*)realloc(idx, cpt*sizeof(int)) ;
     idx[cpt-1]=res ;
     if (IDS.isTCF(res)) nbelem[cptvar-1]=&ncell ;
     else nbelem[cptvar-1]=&natm ;
    }
  }
  idx=(int*)realloc(idx, (cpt+1)*sizeof(int)) ; idx[cpt]=-1 ;


  // reordering par ordre croissant (les différentes composantes d'une variable de dimension 2 ou 3 se suivent toujours !!!
  int deplacements, tmp ; int * tmpptr ; char * tmpstr ;
  do
  {
    deplacements=0 ; k=1 ;
    for (j=dims[0] ; j<cpt ; j+=dims[k-1])
    {
     if (idx[j]<idx[j-1]) // il faut permuter, et ça c'est pas drôle
      { deplacements++ ;
	for (l=0 ; l<dims[k] ; l++)
	{
	  for (int m=0 ; m<dims[k-1] ; m++)
	  {
	   tmp=idx[j+l-m] ; idx[j+l-m]=idx[j+l-m-1] ; idx[j+l-m-1]=tmp ;
	  }
	}
	j-=dims[k-1] ;
	tmp=dims[k] ; dims[k]=dims[k-1] ; dims[k-1]=tmp ;
	j+=dims[k-1] ;
	tmpptr=nbelem[k] ; nbelem[k]=nbelem[k-1] ; nbelem[k-1]=tmpptr ;
	tmpstr=noms[k] ; noms[k]=noms[k-1] ; noms[k-1]=tmpstr ;
      }
      k++ ;
    }
  } while (deplacements>0) ;
  cfiddeb=-1 ;
  for (j=0 ; j<cpt ; j++) {if (IDS.isTCF(idx[j])) {cfiddeb=j ; break ; }}
  datas=(double **)malloc(cptvar*sizeof(double *)) ;
  d=(double **) malloc(cpt*sizeof(double*)) ;

  sendout(FIRSTTS) ; pthread_cond_wait(&sigin, &mutex);
  sendout(ASKINGGRID) ; pthread_cond_wait(&sigin, &mutex);
  if (Signal != UNABLE)
   { boite=(double *)calloc(3, sizeof(double)) ;
   boite[0]=(double)(*i) ; boite[1]=(double)(*(i+1)) ; boite[2]=(double)(*(i+2)) ; }

  sendout(FIRSTTS) ; pthread_cond_wait(&sigin, &mutex);
  ts=0 ;
  while (Signal==OK)
  {
    int tmp, val ;
    ncell=-1 ; natm=-1 ;
    for (int dos=0 ; dos<2 ; dos++)
    {
    if (cfiddeb==-1) {if(dos==1) continue ; sendout(CHOOSELDUMP) ; pthread_cond_wait(&sigin, &mutex); i=idx ;}
    else if (cfiddeb== 0) {if(dos==1) continue ; sendout(CHOOSECFDUMP) ; pthread_cond_wait(&sigin, &mutex); i=idx ;}
    else {
        if (dos==0) {sendout(CHOOSECFDUMP) ; pthread_cond_wait(&sigin, &mutex); i=idx+cfiddeb ; }
        else {for (k=cpt-cfiddeb-1 ; k>=0 ; k--) d[cfiddeb+k]=d[k] ;
	  ncell=val ;
	  sendout(CHOOSELDUMP) ; pthread_cond_wait(&sigin, &mutex); i=idx ; tmp=idx[cfiddeb] ; idx[cfiddeb]=-1 ; }
      }
    sendout(ASKINGND) ; pthread_cond_wait(&sigin, &mutex);
    val=*i ;
    }
    if (cfiddeb>0) {idx[cfiddeb]=tmp ; natm=val ; }
    if (cfiddeb==-1) {natm=val ;}
    if (cfiddeb==0) {ncell=val ;}

    if (ts==0)
    {
      ncell_orig=ncell ; natm_orig=natm ;
      if (ncell_orig==0 || natm_orig==0) DISP_Err(_("Writing : des données ont été demandées mais sont vides pour ce ts.")) ;
    }
    else
    {
      if (ncell!=ncell_orig || natm!=natm_orig)
	DISP_Err(_("Writing : pas le même nombre d'atomes ou de cellules entre ts. Impossible d'écrire ASCII")) ;
    }

    cpt=0 ;
    for (j=0 ; j<cptvar ; j++)
    {
      if (ts==0) {datas[j]=(double *)calloc((*nbelem[j])*(ts+1)*dims[j], sizeof(double)) ; }
      else 	  datas[j]=(double *)realloc(datas[j], (*nbelem[j])*(ts+1)*dims[j]*sizeof(double)) ;
      if (datas[j]==NULL) {DISP_Err(_("Writing: pas assez de mémoire à allouer")) ; }
      for (k=0 ; k<dims[j] ; k++)
      {
	for (l=0 ; l<(*nbelem[j]) ; l++)
	{datas[j][ts* (*nbelem[j]) *dims[j]+k*(*nbelem[j])+l]=d[cpt][l] ; }
        cpt++ ;
      }
    }
    ts++ ;
    sendout(NEXTTS) ; pthread_cond_wait(&sigin, &mutex);
  }

// Let's find the best way to write
bool multifile =false ;
if (ts > 1) multifile=true ;
else
{
  for (j=1 ; j<cptvar ; j++)
  {
    if (*nbelem[j]!=*nbelem[j-1] ) {multifile=true ; break ;}
    if (dims[j]>1) {multifile=true ; break ;}
  }
}

if (multifile)
{
 for (j=0 ; j<cptvar ; j++)
 {
   if (dims[j]>1 && ts>1 && *nbelem[j]>1) {DISP_Warn(_("Writing : impossible d'écrire une variable de cette dimension en ascii (trop de ts / dim >1).\n")) ; continue ; }
   if (ts==1)
   {
    if (type==1)
    {
      chemin=nom+"-"+noms[j]+".txt" ;
      out=fopen (chemin.c_str(), "w") ;
      ASCIIwrite(out, datas[j], (*nbelem[j]), dims[j]) ;
    }
    else
    {
      chemin=nom+"-"+noms[j]+".dat" ;
      out=fopen (chemin.c_str(), "wb") ;
      BINARYwrite(out, datas[j], (*nbelem[j]), dims[j]) ;
      printf(_("\n\nDimensions du fichier binaire %s : %dx%d\n"), chemin.c_str(), dims[j],(*nbelem[j])) ;
    }
    fclose(out) ;
   }
   else if (dims[j]==1)
   {
    if (type==1)
    {
      chemin=nom+"-"+noms[j]+".txt" ;
      out=fopen (chemin.c_str(), "w") ;
      ASCIIwrite(out, datas[j], (*nbelem[j]), ts) ;
    }
    else
    {
      chemin=nom+"-"+noms[j]+".dat" ;
      out=fopen (chemin.c_str(), "wb") ;
      BINARYwrite(out, datas[j], (*nbelem[j]), ts) ;
      printf(_("\n\nDimensions du fichier binaire %s : %dx%d\n"), chemin.c_str(), ts,(*nbelem[j])) ;
    }
    fclose(out) ;
   }
   else if (*nbelem[j]==1)
   {
    if (type==1)
    {
      chemin=nom+"-"+noms[j]+".txt" ;
      out=fopen (chemin.c_str(), "w") ;
      ASCIIwrite(out, datas[j], dims[j], ts) ;
    }
    else
    {
      chemin=nom+"-"+noms[j]+".dat" ;
      out=fopen (chemin.c_str(), "wb") ;
      BINARYwrite(out, datas[j], dims[j], ts) ;
      printf(_("\n\nDimensions du fichier binaire %s : %dx%d\n"), chemin.c_str(), ts,dims[j]) ;
    }
    fclose(out) ;
   }
 }
}
else
{
 if (type==1)
 {
   chemin=nom+"-"+"ascii"+".txt" ;
   out=fopen (chemin.c_str(), "w") ;
   fprintf(out, "%%# ") ;
   for (j=0 ; j<cptvar ; j++) fprintf(out, "%s ", noms[j]) ;
   fprintf(out, "\n") ;
   ASCIIwrite(out, datas, cptvar, (*nbelem[0])) ;
 }
 else
 {
   chemin=nom+"-"+"binary"+".dat" ;
   out=fopen (chemin.c_str(), "wb") ;
   printf(_("\n\n%%#Liste des variables du fichier binaire : \n%%#")) ;
   for (j=0 ; j<cptvar ; j++) printf("%s ", noms[j]) ;
   printf("\n") ;
   BINARYwrite(out, datas, cptvar, (*nbelem[0])) ;
 }
 fclose(out) ;
}

// Cleaning
for (j=0 ; j<cptvar ; j++) {free(datas[j]) ; free(noms[j]) ; }
free(datas) ; free(boite) ; free(idx) ; free (dims) ; free(nbelem) ; free(noms) ;
return 1 ;
}
//---------------------------------------------------------------------
int Writing::ASCIIwrite (FILE * out, double *datas, int lig, int col)
{
  int i, j ;
  for (i=0 ; i<col ; i++)
  {for (j=0 ; j<lig ; j++)
    {
     fprintf(out, "%g", datas[i*lig+j]) ;
     if (j<lig-1) fprintf(out, ",") ;
    }
    fprintf(out, "\n") ; }
return 1 ;
}

int Writing::ASCIIwrite (FILE * out, double **datas, int ndim, int n)
{
  int i, j ;
  for (i=0 ; i<n ; i++)
  {for (j=0 ; j<ndim ; j++)
    {fprintf(out, "%6g ", datas[j][i]) ; }
    fprintf(out,"\n") ; }
return 1 ;
}

int Writing::ASCIIwrite (ofstream & out, double **datas, int ndim, int n)
{
  int i, j ;
  for (i=0 ; i<n ; i++)
  {for (j=0 ; j<ndim ; j++)
    {out << datas[j][i] << " " ; }
    out << "\n" ; }
return 1 ;
}

int Writing::BINARYwrite (FILE * out, double *datas, int lig, int col)
{
  int i, j ;
  for (i=0 ; i<col ; i++)
  {for (j=0 ; j<lig ; j++)
     fwrite(datas+i*lig+j,sizeof(double),1,out) ; //out << datas[j][i] << " " ;
  }
return 1 ;
}

int Writing::BINARYwrite (FILE * out, double **datas, int ndim, int n)
{
  int i, j ;
  for (i=0 ; i<n ; i++)
  {for (j=0 ; j<ndim ; j++)
    {fwrite(&(datas[j][i]),sizeof(double),1,out) ; }
  }
return 1 ;
}


//=========================== TIFF WRITING FUNCTION =======================
#ifdef USETIFF
int Writing::TIFF_bw_write (string nom, int w, int h, double **datas)
{
TIFF *out= TIFFOpen(nom.c_str(), "w");
int sampleperpixel = 1;
int i, j, k ;
unsigned char *image=new unsigned char [w*h*sampleperpixel];

TIFFSetField (out, TIFFTAG_IMAGEWIDTH, w);  // set the width of the image
TIFFSetField(out, TIFFTAG_IMAGELENGTH, h);    // set the height of the image
TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, sampleperpixel);   // set number of channels per pixel
TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);    // set the size of the channels
TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    // set the origin of the image.
//   Some other essential fields to set that you do not have to understand for now.
TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG); //?????
TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);       // ??????

tsize_t linebytes = sampleperpixel * w;     // length in memory of one row of pixel in the image.

for (i=0 ; i<w ; i++)
  for (j=0 ; j<h ; j++)
   for (k=0 ; k<sampleperpixel ; k++)
      image[i*sampleperpixel+j*w*sampleperpixel+k]= (unsigned char)(datas[i][j]*255) ;

unsigned char *buf = NULL;        // buffer used to store the row of pixel information for writing to file
//    Allocating memory to store the pixels of current row
if (TIFFScanlineSize(out)<linebytes)
    buf =(unsigned char *)_TIFFmalloc(linebytes);
else
    buf = (unsigned char *)_TIFFmalloc(TIFFScanlineSize(out));

// We set the strip size of the file to be size of one row of pixels
TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(out, w*sampleperpixel));

//Now writing image to the file one strip at a time
for (i = 0; i < h; i++)
{
    memcpy(buf, &image[(h-i-1)*linebytes], linebytes);    // check the index here, and figure out why not using h*linebytes
    if (TIFFWriteScanline(out, buf, i, 0) < 0)
    break;
}
(void) TIFFClose(out);
if (buf) _TIFFfree(buf);

free(image) ;
return 0 ;
}

#endif
