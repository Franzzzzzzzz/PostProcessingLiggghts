#include "Headers/Dump.h"
#include <sys/stat.h>

// =====================================================
// Fonctions de la superclasse Dump ====================
//======================================================
Dump::Dump (void)
{
 nbsteps=0 ;
 if (actions["dumpall"].set)
   inscrit=-2 ;
 else
   inscrit=-1 ;
}

//-------------------------------------------------------
int Dump::open (string chem)
{
struct stat ficdmp, fictmp ;
string chemtmp(chem) ;
fstream dumptmp ;
int i, totalsize ;

chemtmp.append(".tmp") ;
dumptmp.open(chemtmp.c_str(), ios::in | ios::binary) ;
if (stat(chemtmp.c_str(), &fictmp)>=0)
 {
 stat(chem.c_str(), &ficdmp);
 if (ficdmp.st_mtime>fictmp.st_mtime)
  {
  DISP_Info(_("Tmp plus ancien que le dump. Recalculé")) ;
  dumptmp.close() ;
  remove(chemtmp.c_str()) ;
  }
 }

dumpin.open(chem.c_str()) ;
if (!dumpin)
  {DISP_Err(_("Erreur d'ouverture du dump\n")) ; return -1 ;}

if (dumptmp.is_open() && actions["chainforce"].set)
 {
 //DISP_Info("L'utilisation des fichiers tmp avec les chainforce est supprimé par précaution (difficile de distinguer le TL du TCF). Le tmp sera supprimé") ;
 //dumptmp.close() ;
 //remove(chemtmp.c_str()) ;
 }

//cout << "Reading TMP a priori" ; fflush(stdout) ;
if (dumptmp.is_open() && !actions["dumpall"].set )//&& !actions["chainforce"].set)
  {
  Step null_step ;
  dumptmp.read((char*)&nbsteps, sizeof(int)) ;
  Type = nbsteps >> (sizeof(int)*8-3) ;
  if (Type==0) DISP_Info(_("Le .tmp est ancien; il ne contient pas de donné de type ...\n")) ;
  else if (Type >3) DISP_Warn(_("La donnée de type est étrange ...")) ;
  nbsteps=nbsteps&(~(3<<(sizeof(int)*8-3))) ;
  steps.resize(nbsteps, null_step) ;
  for (i=0 ; i<nbsteps ; i++)
      {
	  // DOES NOT WORK !
	  //std::size_t found = chem.find("force");
	  //if (found!=std::string::npos) steps[i].Type=TCF ;
	  //else steps[i].Type=TL ;
      //if (actions["vtk"].set) steps[i].Type=TF ;
      //else if (actions["chainforce"].set) {DISP_Err("L'utilisation des fichiers tmp avec les chainforce est supprimé par précaution (difficile de distinguer le TL du TCF).") ;
      //	  	  	  	  	  	  	  	  	  steps[i].Type=TCF ; }
      //else steps[i].Type=TL ;

      dumptmp.read((char*)&steps[i].posinfile, sizeof(streampos)) ;
      dumptmp.read((char*)&steps[i].timestep, sizeof(int)) ;
      steps[i].Type=Type ;
      steps[i].multisphere = & multisphere ;
      }
  dumptmp.close() ;
  return 0 ;
  }
else
  {
  cout << "\nDump::open          " ;
  actions.disp_progress() ;

  actions.total=dumpin.size() ;

  if (actions["dumpall"].set)
    read(5, -1) ;
  else
    read(4, -1) ;

  totalsize=nbsteps*(sizeof(streampos)+sizeof(int)) ;
  if (totalsize/(double)(dumpin.size())>0.1 && !actions["forcetmp"].set)
  {
    DISP_Info(_("Le tmp ne sera pas ecrit car trop gros comparé au dump. Utiliser --forcetmp pour contrecarrer cela.\n")) ;
  }
  else
  {
  dumptmp.open(chemtmp.c_str(), ios::out | ios::binary) ;

  int nbstepsmod ;
  if (Type>0 && Type < 4) {nbstepsmod=(nbsteps | (Type << (sizeof(int)*8-3))) ; }
  if (!dumptmp.is_open()) cout << _("WARNING : impossible de créer le fichier tmp de dump\n") ;
  dumptmp.write((char*)&nbstepsmod, sizeof(int)) ;
  //dumptmp.write((char*)&nbsteps, sizeof(int)) ;
  for (i=0 ; i<nbsteps ; i++)
     {
     dumptmp.write((char*)&steps[i].posinfile, sizeof(streampos)) ;
     dumptmp.write((char*)&steps[i].timestep, sizeof(int)) ;
     }
  dumptmp.close() ;
  }
  return 1 ;
  }
}
//-----------------------------------------------------------------
int Dump::write_asVTKBase(string chemin, char type)
{
long int loop[3] ;
int ts, i;
Writing * W ;
bool stop ;

if (type == 1) {DISP_Err (_("L'écriture d'un seul dump n'a pas été reprogrammée (pb pour la décompression de dump.stress")) ; return -1 ; } //TODO

W=&(actions.ecrire) ;
W->chemin=&(chemin) ;

pthread_mutex_init(&W->mutex, NULL);
pthread_cond_init (&W->sigin, NULL);
pthread_cond_init (&W->sigout, NULL);
pthread_mutex_lock(&W->mutex);
W->startall() ;

stop=false ;
loopdat(loop) ;
cout << "LDump::write       " ;
actions.set_progress(loop) ; actions.disp_progress() ;

while (stop==false)
{
  pthread_cond_wait(&W->sigout, &W->mutex);
  switch (W->Signal)
  {
    case WHATAREYOU : W->sendin(LDUMP) ; break ;
    case FIRSTTS : ts=loop[0] ; W->i=&ts ; check_timestep(ts) ;
                   W->sendin(OK) ; break ;
    case ASKINGTS : W->i=&steps[ts].timestep ; W->sendin(OK) ; break ;
    case CHOOSELDUMP : W->sendin(OK) ; break ; //Always the case ...
    case CHOOSECFDUMP : W->sendin(UNABLE) ;  break ;
    case NEXTTS :
      ts+=loop[1] ; actions.valeur=ts ;
      if (ts>=loop[2])
	{ts=-1 ; W->i=&ts ; W->sendin(FINISH) ; }
      else
	{check_timestep(ts) ; W->i=&ts ; W->sendin(OK) ; }
      break ;
    case ASKING1D :
      W->d[0]=&(steps[ts].datas[steps[ts].find_idx(*(W->i))][0]) ; W->i=&steps[ts].nb_atomes ; W->sendin(OK) ; break ;
    case ASKING2D :
      for (i=0 ; i<3 ; i++)
        W->d[i]=&(steps[ts].datas[steps[ts].find_idx(*(W->i+i))][0]) ;
      W->i=&steps[ts].nb_atomes ; W->sendin(OK) ; break ;
    case ASKING3D :
      for (i=0 ; i<9 ; i++)
	W->d[i]=&(steps[ts].datas[steps[ts].find_idx(*(W->i+i))][0]) ;
      W->i=&steps[ts].nb_atomes ; W->sendin(OK) ; break ;
    case ASKINGND :
      for (i=0 ; *(W->i+i)!=-1 ; i++)
	W->d[i]=&(steps[ts].datas[steps[ts].find_idx(*(W->i+i))][0]) ;
      W->i=&steps[ts].nb_atomes ; W->sendin(OK) ; break ;
    case ASKINGMULTISPHERE :
      multisphere.prepare_Writing(steps[ts]) ;
      W->d=(double **)(&(multisphere)) ; // Trick: je passe un pointeur class multisphere comme si c'était un pointeur de pointeur de double c'est pas beau mais ça marche :)
      W->i=&(multisphere.ngp) ; W->sendin(OK) ;
      break ;
    case FINI : stop=true ; W->sendin(OK) ; break ;
    case ASKINGGRID : W->sendin(UNABLE) ; break ;
    default : DISP_Warn(_("Signal d'écriture inconnu")); W->sendin(UNABLE) ; W->disp_Signal() ; break ;
  }
}
return 1 ;
}

//------------------------------------------------------------------
void Dump::disp(void)
{
 cout << _("\nRésumé du dump :\n") ;
 cout << nbsteps << _("pas de dump.\n--------------------\n") ;
}
//------------------------------------------------------------------
int Dump::check_timestep(long int timestep)
{
int i ;
//if (actions["dumpall"].set) return 0 ;
if (inscrit==timestep) return 1 ;
if (inscrit>=0 && !actions["dumpall"].set) // il faut libérer la mémoire
  {
  inscrit=free_timestep(inscrit) ;
  }
if (!actions["dumpall"].set)
{
  dumpin.seekg(steps[timestep].posinfile, ios::beg) ;
  read(3, (int)timestep) ;
}

if (actions["multisphere"].set) multisphere.set_current_step(timestep) ;

if (steps[timestep].filtered==false)
{
  Filter filtre_tmp ;
  // ------ Application des pré-filtres ------

  if (prefiltre.size()>0)
  {
   //printf("PRE ") ; filtre_tmp.disp(prefiltre) ;
   filtre_tmp.do_filter(steps[timestep], prefiltre) ;
   //printf("PRE ") ; filtre_tmp.disp(prefiltre) ;
  }
  //------ Application des filtres globaux (ie. issus de la ligne de commande) ----
  if (actions.filtre_global.size()>0)
  {
   filtre_tmp.do_filter(steps[timestep], actions.filtre_global) ;
  //printf("GLOB%d %d ", steps[timestep].Type , steps[timestep].nb_atomes) ; filtre_tmp.disp(actions.filtre_global) ;
  }
  if (actions["anglebystep"].set)
  {
   double angle ; char filtrestr[50] ;
   vector <struct Op> tmp_op ;
   angle=-actions["anglebystep"]["anglebystep"]/180.*M_PI*timestep;
   // Création du filtre pour le dump
   sprintf(filtrestr, "null::zrotate::%.10f", angle) ;
   DISP_Err(_("Please check the type value here (anglebystep in Dump.cpp, check_timestep)\n")) ;
   tmp_op=filtre_tmp.parse_arg(filtrestr,Type) ;
   //tmp_op.valeur=angle ; //Fait automatiquement normalement
   filtre_tmp.do_filter(steps[timestep], tmp_op) ;
  }
  //------- Application des post-filtres --------
  if (postfiltre.size()>0)
  {
   filtre_tmp.do_filter(steps[timestep], postfiltre) ;
   //printf("POST ") ; filtre_tmp.disp(postfiltre) ;
  }
steps[timestep].filtered=true ;
}

if (actions["teststress"].set)
   {
	DISP_Err(_("teststress doit être reprogrammé correctement pour être utilisé\n")) ;
	/*
	int  i ;
	for (i=0 ; i<steps[timestep].nb_triangles ; i++)
	{ Vector centre ;
	centre=steps[timestep].get_tri_center(i) ;

    steps[timestep].datas[steps[timestep].find_idx(FORCEX)][i]=abs(centre(3)) ;
    steps[timestep].datas[steps[timestep].find_idx(FORCEY)][i]=0 ;
    steps[timestep].datas[steps[timestep].find_idx(FORCEZ)][i]=centre(3) ;

	}*/
   }

//--------- Fin des filtres --------
inscrit=timestep ;
return 2 ;
}

//-----------------------------------------------------------------
int Dump::free_timestep (long int inscrit)
{
int i ;

for (i=0 ; i<steps[inscrit].nb_idx ; i++)
     {steps[inscrit].datas[i].clear() ; }
  steps[inscrit].datas.clear() ;
  steps[inscrit].idx_col.clear() ;
  steps[inscrit].nb_idx=0 ;
  steps[inscrit].filtered=false ;
  inscrit=-1 ;
return inscrit ;
}
//---------------------------------------------------------------
void Dump::loopdat (long int *res)
{
if (actions["extract"].set)
   { res[0]=(long int)actions["extract"]["extract_deb"] ; res[2]=(long int)actions["extract"]["extract_fin"] ;
     if (isinf(actions["extract"]["extract_fin"])) {res[2]=nbsteps ;}
     if (res[0]<0) {res[0]=0 ; DISP_Warn(_("WARNING : le debut de l'extraction est <0. RAZ.\n")) ; }
     if (res[2]>nbsteps) {res[2]=nbsteps ; DISP_Warn(_("WARNING : la fin de l'extraction est >nbsteps. Remise à nbsteps.\n")) ; }
     if (res[0]>nbsteps) {DISP_Err(_("ERR: step initial demandée supérieure au nombre de steps\n\n")) ; exit(1) ; }
     if (res[2]<0) {DISP_Err(_("ERR: step final < 0\n\n")) ; exit(1) ;}
  }
else
   {res[0]=0 ; res[2]=nbsteps;}

if (actions["downsampling"].set)
   { res[1]=(long int)actions["downsampling"]["downsampling"] ;
   if (res[1]<=0 || res[1]>nbsteps) { res[1]=1 ; DISP_Warn(_("WARNING : downsampling incorrect. Remis à 1.\n")) ;}
   }
else
	res[1]=1 ;
}
