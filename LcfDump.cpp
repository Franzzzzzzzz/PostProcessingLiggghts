#include "Headers/LcfDump.h"

// =====================================================
// Fonctions de la sous classe LcfDump =================
//======================================================
int LcfDump::read (unsigned char read_type, int index)
{
// read_type==1 : sauver les données lues
 //* read_type==2 : ne lire qu'un dump
 //* read_type==4 : creer la strcuture step et incrémenter nb_step

  string ligne;
  Step null_step ;
  int i, j, type, stop=0, idx ;
  static bool first=true ; 
  char *ligneC, *actu, *suiv ; static bool first_naninf=true ; 
  
  Type=TCF ; 
  ligneC=(char*)malloc(5000) ; 
  
  if (read_type & 4) {nbsteps=0; printf("Lecture initiale du dump, peut prendre du temps ...\n") ; }
  if (read_type & 2) idx=index ;

  while (!dumpin.eof() && stop<2 )
  {
   dumpin.getline(ligne) ;
   type=isitem(ligne) ;
   switch (type)
   {
     case 1 : if (read_type & 4)
		{
		 nbsteps++ ; // Create new step and read 1 line to have timestep number.
		 steps.resize(nbsteps, null_step) ;
		 steps[nbsteps-1].nb_atomes=0 ; steps[nbsteps-1].Type=TCF ;
		 steps[nbsteps-1].posinfile=dumpin.tellg()-(streampos)1 -(streampos)ligne.length() ;
                 if (index==-1) idx=nbsteps-1 ;
                 actions.valeur=dumpin.tellg() ;
		}
	      if (read_type & 2) {stop++ ; if (stop>=2) break ; }
	      dumpin >> steps[idx].timestep ;
	      break ;

     case 2 : if (read_type & 1) {dumpin >> steps[idx].nb_atomes ; }
     //printf("       [%d|%d]",idx, steps[idx].nb_atomes ) ; fflush(stdout) ;
	      break ;

     case 3 : if (read_type & 1)
	      {
	       dumpin >> steps[idx].box[0][0] >> steps[idx].box[0][1] ;
	       dumpin >> steps[idx].box[1][0] >> steps[idx].box[1][1] ;
	       dumpin >> steps[idx].box[2][0] >> steps[idx].box[2][1] ;
	      }
	      break ;
	      
     case 4 : if (read_type & 1)
	      {
               //cout << "o" ; fflush(stdout) ;
	       sparselabels(steps[idx], ligne) ; // Nouvelles données atomes. On continue la lecture de la ligne et on écrit la liste d'index
	      
	       // On alloue le tableau de données nécessaire
	       steps[idx].datas.resize(steps[idx].nb_idx, vector<double>(1,0)) ;
	       for (i=0 ; i<steps[idx].nb_idx; i++)
	       {steps[idx].datas[i].resize(steps[idx].nb_atomes, 0) ;}
	       // On lit l'ensemble des données que l'on range
	       
	       // TEST optimize version
	       for (i=0 ; i<steps[idx].nb_atomes ; i++)
	        {
		     dumpin.getline(&ligneC,5000) ;
	             actu=ligneC ; 
		     for (j=0 ; j<steps[idx].nb_idx ; j++)
		       {
			steps[idx].datas[j][i]=strtod(actu, &suiv) ;
			if (actu==suiv) 
			{ 
			  DISP_Err("Conversion to double failed ! Something really bad happened. Going to next field...\n") ; 
			  while (*actu!=' ' && *actu!=0) actu++ ; 
			}
			else 
			  actu=suiv ; 
			if ( !isfinite(steps[idx].datas[j][i]))
			{
			 if (actions["nofileerror"].set && first)
			 { DISP_Warn("Une valeur lue est inf ou nan \n") ; 
			   first=false ; }
			 else if (!actions["nofileerror"].set && first_naninf)
			 {
                           
			  DISP_Warn("\n\n!!! Une valeur lue est inf ou nan.  ") ; 
			  DISP_Warn("Elle n'a PAS été supprimée.\n") ;
                          first_naninf=false ; 
			  //steps[idx].datas[j][i]=0 ;
			 }
			}
		       }
		}
	       
	       
	       /*for (i=0 ; i<steps[idx].nb_atomes ; i++)
	        {
		     for (j=0 ; j<steps[idx].nb_idx ; j++)
		       {
			int st=dumpin.rdstate() ;  
			dumpin >> steps[idx].datas[j][i];
			if (dumpin.fail())
			{
			 dumpin.clear() ; 
			 if (actions["nofileerror"].set && first)
			 { DISP_Warn("Une erreur de lecture est apparue (LCF). nofileerror étant activé, les erreurs suivantes seront silencieuses.\n") ; 
			   first=false ; }
			 else if (!actions["nofileerror"].set)
			 {
			  DISP_Warn("Une erreur de lecture (LCF) est apparue à l'octet : ") ; 
			  printf("%d (error number %d, ts %d). ", (int)dumpin.tellg(), st, index) ;
			  DISP_Warn("Elle a été supprimée et la valeur est mise à 0.\n") ;
			  steps[idx].datas[j][i]=0 ; exit (1) ; 
			 }
			}
		       }
		}*/
	       } break ;
   }   
  }
  
dumpin.clear() ;
return 0 ;
}
//-------------------------------------------------
int LcfDump::isitem (string ligne)
{
size_t position ;

position=ligne.find("ITEM: ") ;
if (position!=string::npos)
  {
  position=ligne.find("TIMESTEP") ;
  if (position!=string::npos)
    return 1 ;
  position=ligne.find("NUMBER OF ENTRIES") ;
  if (position!=string::npos)
    return 2 ;
  position=ligne.find("BOX BOUNDS") ;
  if (position!=string::npos)
    return 3 ;
  position=ligne.find("ENTRIES") ;
  if (position!=string::npos)
   return 4 ;

  cerr << "WARN1 : item inconnu " << ligne << ".\n" ;
  return -1 ;
  }
else
   return 0 ;
}
//-------------------------------------------------------------
int LcfDump::sparselabels(Step &step , string ligne)
{
  int compteur=0, res ;
  size_t espace ;
  string ligne2, item ;

  ligne2=ligne.substr(14) ;
  while (1)
  {
    espace=ligne2.find_first_of(' ') ;
    if (espace==ligne2.npos) break ;
    item=ligne2.substr(0, espace) ;
    ligne2=ligne2.substr(espace+1) ;

    step.idx_col.resize(compteur+1) ;
    res=IDS(item) ;
    
    if (res==-1) 
    {DISP_Warn ("Un type de donné n'est pas référencé, IDS va l'ajouter\n") ; 
     res=IDS.new_id(item, TCF) ; 
     if (res!=-1) step.idx_col[compteur]=res ;
     else 
     {
       DISP_Err("Le type inconnu n'a pas pu être ajouté !!\n"); 
       step.idx_col[compteur]=IDS("UNKNOWN") ; 
     }
    }
    else
     step.idx_col[compteur]=res ;

    compteur++ ;
  }
  step.nb_idx=compteur ;
return 1 ;   
}
/*    if (actions["noperiodicity"].set)
     {
     if (item=="c_cout[1]")            step.idx_col[compteur]=CFID1 ;
     else if (item == "c_cout[2]")      step.idx_col[compteur]=CFID2 ;
     else if (item == "c_cout[3]")      step.idx_col[compteur]=CFFORCEX ;
     else if (item == "c_cout[4]")      step.idx_col[compteur]=CFFORCEY ;
     else if (item == "c_cout[5]")     step.idx_col[compteur]=CFFORCEZ ;
     else {cerr << "WARN2cf : le type de données" << item <<"est inconnu\n" ; step.idx_col[compteur]=UNKNOWN ; }
     }
    else
     {
     if (item=="c_cout[1]")            step.idx_col[compteur]=CFID1 ;
     else if (item == "c_cout[2]")      step.idx_col[compteur]=CFID2 ;
     else if (item == "c_cout[3]")      step.idx_col[compteur]=CFPERIOD ;
     else if (item == "c_cout[4]")      step.idx_col[compteur]=CFFORCEX ;
     else if (item == "c_cout[5]")      step.idx_col[compteur]=CFFORCEY ;
     else if (item == "c_cout[6]")     step.idx_col[compteur]=CFFORCEZ ;
     else if (item == "c_contacts[7]")  step.idx_col[compteur]=CFID1 ;
     else if (item == "c_contacts[8]")      step.idx_col[compteur]=CFID2 ;
     else if (item == "c_contacts[9]")      step.idx_col[compteur]=CFPERIOD ;
     else if (item == "c_contacts[1]")      step.idx_col[compteur]=CFID1X ;
     else if (item == "c_contacts[2]")      step.idx_col[compteur]=CFID1Y ;
     else if (item == "c_contacts[3]")      step.idx_col[compteur]=CFID1Z ;
     else if (item == "c_contacts[4]")      step.idx_col[compteur]=CFID2X ;
     else if (item == "c_contacts[5]")      step.idx_col[compteur]=CFID2Y ;
     else if (item == "c_contacts[6]")      step.idx_col[compteur]=CFID2Z ;
     else if (item == "c_contacts[10]")      step.idx_col[compteur]=CFFORCEX ;
     else if (item == "c_contacts[11]")      step.idx_col[compteur]=CFFORCEY ;
     else if (item == "c_contacts[12]")     step.idx_col[compteur]=CFFORCEZ ;
     else {cerr << "WARN2cf : le type de données" << item <<"est inconnu\n" ; step.idx_col[compteur]=UNKNOWN ; }
     }*/

//------------------------------------------------------
int LcfDump::write_asVTK (string chemin, LucDump &dump)
{
long int loop[3] ;
int ts, i ; 
vector <int> idx; 
Writing * W ;
bool stop ; 
Dump *dmp ; 

W=&(actions.ecrire) ; 
W->chemin=&(chemin) ; 

pthread_mutex_init(&W->mutex, NULL);
pthread_cond_init (&W->sigin, NULL);
pthread_cond_init (&W->sigout, NULL);
pthread_mutex_lock(&W->mutex);
W->startall() ; 

stop=false ; 
loopdat(loop) ;
cout << "LcfDump::write       " ;

actions.total=loop[2]-loop[0] ; actions.disp_progress() ;
while (stop==false) 
{
  pthread_cond_wait(&W->sigout, &W->mutex);
  switch (W->Signal)
  {
    case WHATAREYOU : W->sendin(CFDUMP) ; break ; 
    case FIRSTTS : ts=loop[0] ; W->i=&ts ; fflush(stdout) ; check_timestep(ts) ; W->sendin(OK) ; fflush(stdout) ; break ;
    case CHOOSELDUMP :  dmp=&(dump) ; W->sendin(OK) ; break ; //Always the case ...
    case CHOOSECFDUMP : dmp=this ; W->sendin(OK) ;  break ; 
    case ASKINGTS : W->i=&(dmp->steps[ts].timestep) ; W->sendin(OK) ; break ; 
    case ASKINGGRID : W->sendin(UNABLE) ; break ; 
    case NEXTTS : 
      ts+=loop[1] ; actions.valeur=ts-loop[0] ; 
      if (ts>=loop[2]) 
	{ts=-1 ; W->i=&ts ; W->sendin(FINISH) ; } 
      else 
	{check_timestep(ts) ; W->i=&ts ; W->sendin(OK) ; }
      break ; 
    case ASKING1D : 
      idx.clear() ;
      idx.push_back(dmp->steps[ts].find_idx(*(W->i))) ;
      if (idx[0]==-1) {W->sendin(UNABLE) ; break ; }
      W->d[0]=&(dmp->steps[ts].datas[idx[0]][0]) ; W->i=&(dmp->steps[ts].nb_atomes) ; W->sendin(OK) ; break ;  
    case ASKING2D : 
      idx.clear() ; 
      for (i=0 ; i<3 ; i++) 
      { idx.push_back(dmp->steps[ts].find_idx(*(W->i+i))) ; 
	if (idx[i]==-1) {W->sendin(UNABLE) ; break ; }
      }
      for (i=0 ; i<3 ; i++) 
        W->d[i]=&(dmp->steps[ts].datas[idx[i]][0]) ; 
      if (i==3) {W->i=&(dmp->steps[ts].nb_atomes) ; W->sendin(OK) ;} break ;  
    case ASKING3D : 
      idx.clear() ; 
      for (i=0 ; i<9 ; i++) 
      {
	idx.push_back(dmp->steps[ts].find_idx(*(W->i+i))) ; 
	if (idx[i]==-1) {W->sendin(UNABLE) ; break ; }
      }
      for (i=0 ; i<9 ; i++)
	W->d[i]=&(dmp->steps[ts].datas[idx[i]][0]) ; 
      if (i==9) {W->i=&(dmp->steps[ts].nb_atomes) ; W->sendin(OK) ;} break ;  
    case ASKINGND : 
      idx.clear() ; 
      for (i=0 ; *(W->i+i)!=-1 ; i++)
      { idx.push_back(dmp->steps[ts].find_idx(*(W->i+i))) ;
	if (idx[i]==-1) {W->sendin(UNABLE) ; break ; }
      }
      for (i=0 ; *(W->i+i)!=-1 ; i++)
         W->d[i]=&(dmp->steps[ts].datas[idx[i]][0]) ; 
      if (*(W->i+i)==-1) {W->i=&(dmp->steps[ts].nb_atomes) ; W->sendin(OK) ;} break ;    
    case FINI : stop=true ; W->sendin(OK) ; break ; 
    default : DISP_Warn("Signal d'écriture inconnu"); W->disp_Signal() ;  break ; 
  }
}
return 1;    
}
//------------------------------------------------------
int LcfDump::write_grainforce (string chemin, LucDump & ldump)
{
string chem, chem2 ;
char infos[500] ;
long int i, j, loop[3], idx[7], sgn, idgrain;
double forces[3] ; double rayon ; 
double ** meanforces ; double *meangrains ; int denom=0 ;
ofstream out, out2 ;
bool byangle, bytot, byduration ; 
map<int,int> contact_duration;
map<int, int> contact_histogram ; int totalcf=0 ; 

//double forces[3] ;

//chem=chemin ; chem.append("-ForceByAngle.txt") ;
//out.open(chem.c_str(), ios::out) ;

loopdat(loop) ;

cout << "\nLcfDump::write_grainforce         " ;
actions.total=loop[2] ; actions.disp_progress() ;
DISP_Info("On travail sur du 2D en xy\n") ;

chem2=chemin ; try{chem2.erase(chem2.rfind(".gz"));} catch(...){} chem2.append("-GrainForce.txt") ;
out2.open(chem2.c_str(), ios::out) ;

// Création des tableaux de données moyennes si nécessaires
if (actions["grainforce-by-angle"].set && actions["mean"].set)
	  {
	  meanforces=(double **)malloc((int)(actions["grainforce-by-angle"]["nbbox_theta"])*sizeof(double *)) ;
	  meangrains=(double *)malloc((int)(actions["grainforce-by-angle"]["nbbox_theta"])*sizeof(double)) ;
	  for (i=0 ; i<actions["grainforce-by-angle"]["nbbox_theta"] ; i++)
	     {meanforces[i]=(double *)malloc(7*sizeof(double)) ; meanforces[i][0]=meanforces[i][1]=meanforces[i][2]=0 ; meangrains[i]=0 ;}
	  }

// Boucle sur les ts
byangle=actions["grainforce-by-angle"].set ; 
bytot=actions["grainforce"].set ; 
byduration=actions["grainforce-duration"].set ; 
if (bytot) idgrain=actions["grainforce"]["id"] ; 
//if (byangle) idgrain=actions["grainforce-by-angle"]["id"] ; // TODO would make sense to do it as well ...
for (i=loop[0] ; i<loop[2] ; i+=loop[1])
{
  actions.valeur=i ;
  check_timestep(i) ;
  
  if (byduration)
    for (multimap<int,int>::iterator it=contact_duration.begin(); it!=contact_duration.end(); it++)
      (*it).second*=-1 ; // Everyone put to negative. The lost contact remain negative after going through the chainforce.
  
  if (bytot)
  {
    idx[0]=steps[i].find_idx(IDS("CFFORCEX")) ; idx[1]=steps[i].find_idx(IDS("CFFORCEY")) ; idx[2]=steps[i].find_idx(IDS("CFFORCEZ")) ;
    idx[3]=steps[i].find_idx(IDS("CFID1")) ; idx[4]=steps[i].find_idx(IDS("CFID2")) ; 
    //idx[5]=ldump.steps[i].find_idx(IDS("RAYON")) ; idx[6]=ldump.steps[i].find_idx(IDS("ID")) ;
    forces[0]=forces[1]=forces[2]=0 ; rayon=0 ; 
    for (j=0 ; j<steps[i].nb_atomes+steps[i].nb_atomes_supp ; j++)
    {
      //if (datas[idx[5]][i]==1) return 0 ;
      if (steps[i].datas[idx[4]][j]==idgrain) sgn=-1;
      else if (steps[i].datas[idx[3]][j]==idgrain) sgn=1 ;  //Tout va bien
      else {DISP_Warn("Il n'y a pas le bon ID dans la chaines !") ; continue ; }
    
      forces[0]+=(steps[i].datas[idx[0]][j]*sgn) ;
      forces[1]+=(steps[i].datas[idx[1]][j]*sgn) ;
      forces[2]+=(steps[i].datas[idx[2]][j]*sgn) ;
      
      // Mesure de la durée des contacts 
      if (byduration)
      {
        if (sgn==1)
        {
          contact_duration[steps[i].datas[idx[4]][j]]*=-1; 
          contact_duration[steps[i].datas[idx[4]][j]]+=1 ; totalcf++ ; 
        }
        else
        {
          contact_duration[steps[i].datas[idx[3]][j]]*=-1 ; 
          contact_duration[steps[i].datas[idx[3]][j]]+=1 ; totalcf++ ; 
        }
      }
    }
    out2 << steps[i].has_periodic_chains << " " << forces[0] <<" " << forces[1] << " " << forces[2] << "\n" ;
  }
  
  if (byangle)
    denom=denom+steps[i].grain_force(meanforces, meangrains, out) ;
  
  if (byduration)
  {
    bool changes=true ;
    while (changes)
    {
      changes=false ; 
      for (multimap<int,int>::iterator it=contact_duration.begin(); it!=contact_duration.end(); it++)
      {
        if ((*it).second<0)
        {
          contact_histogram[-(*it).second]++ ; 
          contact_duration.erase(it) ; changes=true ; 
        }
      }
    }
  }
  
}

if (bytot) out2.close() ; 
//out.close() ;

//DISP_Info("Nombre de timestep otés en raison de traversée de pbc : ") ;
//printf("%d %d\n", (int) ((loop[2]-loop[0])/(double)loop[1] - denom), denom) ;
if (byangle)
{
  sprintf(infos, "La largeur en angle des boites est de %f degrés (%f radians).\n", 360.0/actions["grainforce-by-angle"]["nbbox_theta"], 2*M_PI/actions["grainforce-by-angle"]["nbbox_theta"]) ;    ;
  DISP_Info(infos) ;

  //sprintf(infos, "La commande matlab à utiliser est reshape(X, %d, %d, %d).\n", 5, (int)actions["grainforce-by-angle"]["nbbox_theta"], (int)denom) ;
  //DISP_Info(infos) ;

  // Création du fichier avec moyennage si besoin
  if (actions["mean"].set)
  {
    chem=chemin ; try{chem.erase(chem.rfind(".gz")) ;} catch(...){}  chem.append("-GrainForceByAngle-mean.txt") ;
    out.open(chem.c_str(), ios::out) ;
    for (i=0; i<actions["grainforce-by-angle"]["nbbox_theta"] ; i++ )
    {
	out << (2*i*M_PI)/actions["grainforce-by-angle"]["nbbox_theta"]+M_PI/actions["grainforce-by-angle"]["nbbox_theta"] << " "  << meangrains[i]/denom << " " ;
	out << meanforces[i][0]/denom << " " << meanforces[i][1]/denom << " " << meanforces[i][2]/denom << "\n" ;
    }
   out.close() ;
   }
}

if (actions["grainforce-by-angle"].set && actions["mean"].set)
{
  for (i=0 ; i<actions["grainforce-by-angle"]["nbbox_theta"] ; i++) free(meanforces[i]) ; 
  free(meanforces) ; free(meangrains) ; 
}

if (byduration)
{
FILE * outhist ; 
DISP_Info("Nombre total de chainforce: ") ; printf("%d\n", totalcf) ; 
chem=chemin ; try{chem.erase(chem.rfind(".gz")) ;} catch(...){} chem.append("-DurationHistogram.txt") ;
out.open(chem.c_str(), ios::out) ;
out << "DureeTS Occurence\n" ;  
for (multimap<int,int>::iterator it=contact_histogram.begin(); it!=contact_histogram.end(); ++it)
  out << (*it).first << " " << (*it).second << "\n" ; 
out.close() ; 
}

return 0 ;
}

//------------------------------------------------------
int LcfDump::write_grainforcetot (string chemin)
{
  
DISP_Err("Outdated function   LcfDump::write_grainforcetot\n") ;
  
/*  
string chem ;
char infos[500] ;
long int i, j, loop[3], idx[6], sgn;
double forces[3] ;

ofstream out ;
//double forces[3] ;

//chem=chemin ; chem.append("-ForceByAngle.txt") ;
//out.open(chem.c_str(), ios::out) ;

loopdat(loop) ;

cout << "\nLcfDump::write_grainforcetot          " ;
actions.total=loop[2] ; actions.disp_progress() ;
DISP_Info("On travail sur du 2D en xy\n") ;

// Boucle sur les ts
chem=chemin ; chem.erase(chem.rfind(".gz")) ; chem.append("-GrainForce.txt") ;
out.open(chem.c_str(), ios::out) ;
int nbstepok=0 ;
for (i=loop[0] ; i<loop[2] ; i+=loop[1])
	{
	actions.valeur=i ;
	forces[0]=forces[1]=forces[2]=0 ;
	check_timestep(i) ;
	//if (steps[i].has_periodic_chains==true) continue ;

	idx[0]=steps[i].find_idx(IDS("CFFORCEX")) ; idx[1]=steps[i].find_idx(IDS("CFFORCEY")) ; idx[2]=steps[i].find_idx(IDS("CFFORCEZ")) ;
	idx[3]=steps[i].find_idx(IDS("CFID1")) ; idx[4]=steps[i].find_idx(IDS("CFID2")) ; //idx[5]=find_idx(IDS("CFPERIOD")) ;
	for (j=0 ; j<steps[i].nb_atomes+steps[i].nb_atomes_supp ; j++)
	 {
		//if (datas[idx[5]][i]==1) return 0 ;
	    if (steps[i].datas[idx[4]][j]==actions["grainforce"]["id"]) sgn=-1;
		else if (steps[i].datas[idx[3]][j]==actions["grainforce"]["id"]) sgn=1 ;  //Tout va bien
		else {DISP_Warn("Il n'y a pas le bon ID dans la chaines !") ; continue ; }

		forces[0]+=(steps[i].datas[idx[0]][j]*sgn) ;
		forces[1]+=(steps[i].datas[idx[1]][j]*sgn) ;
		forces[2]+=(steps[i].datas[idx[2]][j]*sgn) ;
	 }

	nbstepok++ ;
	out << steps[i].has_periodic_chains << " " << forces[0] <<" " << forces[1] << " " << forces[2] ;
	out << "\n" ;
	}
out.close() ;

char chaine[500] ;
sprintf(chaine, "Nombre de ts utilises : %d sur %ld.", nbstepok, (loop[2]-loop[0])/loop[1]) ;
DISP_Info(chaine) ;
*/
return 0 ;
}
//==============================================
int LcfDump::write_radiuscontact (void)
{

char infos[500] ;
long int i, j, loop[3], idx[2], idxOK, sgn;
double radius=0, *rayons; int *ids ;
FILE *in ;

//double forces[3] ;

//chem=chemin ; chem.append("-ForceByAngle.txt") ;
//out.open(chem.c_str(), ios::out) ;

loopdat(loop) ;

cout << "\nLcfDump::write_radiuscontact          " ;
actions.total=loop[2] ; actions.disp_progress() ;
DISP_Info("On travail sur du 2D en xy\n") ;

// Chargement du tableau des rayons (compliqué sinon...)
DISP_Warn("On utilise un fichier externe de rayons Radius.txt. Ce n'est certes pas joli, mais beaucoup plus simple !!\n") ;
in=fopen("Radius.txt", "r") ;
double tmpray ; int id, nbids=0, tmpid, idold=-1 ; int res ; 
ids=(int *)malloc(sizeof(int)) ; rayons=(double *)malloc(sizeof(double)) ;
while (!feof(in))
 {
 res=fscanf(in, "%d %lf\n", &id, &tmpray) ;
 nbids++;
 ids=(int *)realloc(ids, nbids*sizeof(int)) ;
 rayons=(double *)realloc(rayons, nbids*sizeof(double)) ;
 ids[nbids-1]=id ;
 rayons[nbids-1]=tmpray ;
 }

for(i=0 ; i<nbids ; i++)
 {
 for (j=0 ; j<nbids ; j++)
	if (ids[j]==i) break ;
 if (j==nbids)
    {nbids++ ;
    ids	  = (int *)realloc(ids, nbids*sizeof(int)) ;
    rayons= (double *)realloc(rayons, nbids*sizeof(double)) ;
    ids[nbids-1]=ids[i] ; ids[i]=i ;
    rayons[nbids-1]=rayons[i] ; rayons[i]=0 ;
    }
 else
 	{
	tmpray=rayons[j] ; rayons[j]=rayons[i] ; rayons[i]=tmpray ;
 	tmpid = ids[j] ;   ids[j]   = ids[i]   ; ids[i]   = tmpid ;
 	}
 }
fclose(in) ;


// Boucle sur les ts

int nbstepok=0 ; int nbgrains=0 ;
for (i=loop[0] ; i<loop[2] ; i+=loop[1])
	{
	actions.valeur=i ;
	check_timestep(i) ;

	if (steps[i].has_periodic_chains==true) continue ;

	idx[0]=steps[i].find_idx(IDS("CFID1")) ; idx[1]=steps[i].find_idx(IDS("CFID2")) ; //idx[5]=find_idx(IDS("CFPERIOD")) ;
	for (j=0 ; j<steps[i].nb_atomes ; j++)
	 {
		//if (datas[idx[5]][i]==1) return 0 ;
	    if      (steps[i].datas[idx[0]][j]==actions["grain-rayon-around"]["id"]) idxOK=steps[i].datas[idx[1]][j];
		else if (steps[i].datas[idx[1]][j]==actions["grain-rayon-around"]["id"]) idxOK=steps[i].datas[idx[0]][j];  //Tout va bien
		else {DISP_Warn("Il n'y a pas le bon ID dans la chaines !") ; continue ; }

	    radius+=rayons[idxOK] ; nbgrains++ ;
	 }

	nbstepok++ ;
	}

char chaine[500] ;
sprintf(chaine, "Nombre de ts utilises : %d sur %ld. Rayon : %.15lf.\n", nbstepok, (loop[2]-loop[0])/loop[1], radius/(double)nbgrains) ;
DISP_Info(chaine) ;

return 0 ;

}

//===========================================
int LcfDump::coupletot (string chemin, LucDump & ldump)
{
  int i ; 
  long int loop[3] ; 
  double omega, **result ; 
  Vector c, caxe, axis, couple ;
  Matrix3x3 Rot ; 
  Cylindre C ; 
  
  loopdat(loop) ;
  result=(double**)malloc(3*sizeof(double*)) ; 
  result[0]=(double*)malloc(sizeof(double)*(loop[2]-loop[0])/loop[1]) ; 
  result[1]=(double*)malloc(sizeof(double)*(loop[2]-loop[0])/loop[1]) ; 
  result[2]=(double*)malloc(sizeof(double)*(loop[2]-loop[0])/loop[1]) ; 
  
  cout << " Entrer la position du centre du cylindre (x y z) (constante) : \n" ;
  cin >> C.centre(1) >> C.centre(2) >> C.centre(3) ;
  //C.centre(1)=C.centre(2)=0 ; C.centre(3)=0.08 ;
  cout << "Coordonnées entrées : " << C.centre(1) << " " << C.centre(2) << " " << C.centre(3) << "\n";
  
  cout << " Entrer l'axe initial du cylindre (x y z) (premier tsdump) : \n" ;
  cin >> caxe(1) >> caxe(2) >> caxe(3) ;
  //caxe(1)=1 ; caxe(2)=caxe(3)=0 ;
  caxe.normalise() ; 
  cout << "Coordonnées entrées : " << caxe(1) << " " << caxe(2) << " " << caxe(3) << "\n";
  
  cout << " Entrer la vitesse de rotation (en degré par tsdump) : \n" ;
  cin >> omega ; 
  //omega=0.9 ; 
  omega=omega/180*M_PI ;
  cout << "Vitesse entrée : " << omega << "\n";
  
  cout << " Entrer l'axe de rotation du cylindre (x, y, z) : \n" ;
  cin >> axis(1) >> axis(2) >> axis(3) ;
  //axis(1)=axis(2)=0 ; axis(3)=-1 ;
  axis.normalise() ; 
  cout << "Coordonnées entrées : " << axis(1) << " " << axis(2) << " " << axis(3) << "\n";
  
  cout << "\nLcfDump::coupletot          " ;
  actions.total=(loop[2]-loop[0])/loop[1] ; actions.disp_progress() ;
  
  C.L=actions.Cst["L"] ; C.D=actions.Cst["D"] ; 
  
  for (i=loop[0] ; i<loop[2] ; i+=loop[1])
	{
	actions.valeur=(i-loop[0])/loop[1] ;
	check_timestep(i) ;
	C.axe=Geometrie::rotation(caxe, (i-loop[0])*omega, axis) ; 
	couple=0 ; 
	steps[i].LCFcouple (ldump.steps[i], couple, C) ;
	result[0][(i-loop[0])/loop[1]]=couple(1) ; 
	result[1][(i-loop[0])/loop[1]]=couple(2) ; 
	result[2][(i-loop[0])/loop[1]]=couple(3) ;
	//printf("%f %f %f\n", couple(1), couple(2), couple(3)) ; 
	}

ofstream out ;
string chemin2 ; 
chemin2=chemin+"-forcecyl.txt" ; 
out.open(chemin2.c_str(), ios::out) ; 
actions.ecrire.ASCIIwrite (out, result, 3, (loop[2]-loop[0])/loop[1]) ; 
return 1; 
}

//--------------------------------------------------------
int LcfDump::forcetank(string chemin, LucDump & ldump)
{
 double * result ;  double *resmean ; int slices;
 int i,j ; long int loop[3] ; int tsmean=0 ; 
 FILE *outx,*outy,*outz, *out ; 
 string chemin2 ; 
 if (actions["mean"].set)
 {chemin2=chemin+"-forcecyl.txt" ; out=fopen(chemin2.c_str(), "w") ;}
 else
 {
   chemin2=chemin+"-forcecyl-x.txt" ; outx=fopen(chemin2.c_str(), "w") ; 
   chemin2=chemin+"-forcecyl-y.txt" ; outy=fopen(chemin2.c_str(), "w") ; 
   chemin2=chemin+"-forcecyl-z.txt" ; outz=fopen(chemin2.c_str(), "w") ; 
 }
 
 result=(double*)malloc(sizeof(double)*actions["forcetank"]["slices"]*3) ; 
 if (actions["mean"].set) 
 {
   resmean=(double*)malloc(sizeof(double)*actions["forcetank"]["slices"]*3) ; 
   for (i=0 ; i<actions["forcetank"]["slices"]*3 ; i++) resmean[i]=0 ; 
 }
 
 loopdat(loop) ;
 actions.total=(loop[2]-loop[0])/loop[1] ; actions.disp_progress() ;  
 for (i=loop[0] ; i<loop[2] ; i+=loop[1])
 {
   actions.valeur=(i-loop[0])/loop[1] ;
   check_timestep(i) ;
   steps[i].LCFforce_by_tank(ldump.steps[i], result) ; 
   slices=(int)(actions["forcetank"]["slices"]) ; 
   if (actions["mean"].set)
   {
     for (j=0 ; j<actions["forcetank"]["slices"]*3 ; j++) resmean[j]+=result[j] ; 
     tsmean++ ; 
   }
   else
   {
     for (j=0 ; j<slices ; j++)
     {
      fprintf(outx, "%g ", result[0*slices+j]) ; 
      fprintf(outy, "%g ", result[1*slices+j]) ; 
      fprintf(outz, "%g ", result[2*slices+j]) ; 
     }
     fprintf(outx, "\n") ; 
     fprintf(outy, "\n") ; 
     fprintf(outz, "\n") ; 
   }
 }
 
 if (actions["mean"].set) 
 {
   for (i=0 ; i<3 ; i++) 
   {
     for (j=0 ; j<actions["forcetank"]["slices"] ; j++)
        fprintf(out, "%g ", resmean[slices*i+j]/tsmean) ; 
     fprintf(out, "\n") ; 
   } 
   fclose (out) ; 
 }
 else
 { fclose(outx); fclose(outy) ; fclose(outz) ; }
 
 return 0 ; 
}






























