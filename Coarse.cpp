#include "Headers/Coarse.h"

// Semaphore et mutex
//sem_t nb_proc_disp ;
//pthread_mutex_t mutex_reading = PTHREAD_MUTEX_INITIALIZER;
//pthread_mutex_t mutex_wait_i = PTHREAD_MUTEX_INITIALIZER;

int CoarseDump::do_coarse (Dump &dump, Dump &dumpatm, int type)
{
int  i ; long int loop[3] ; 
Step stepnull ; 
Coarse coarse_null; 

dumpatm.loopdat(loop) ; 

nbsteps=(loop[2]-1 -loop[0])/loop[1] +1 ;

actions.set_progress(loop) ; 
cout << "\nCoarseDump::do_coarse          " ; actions.disp_progress() ; 
coarse.resize(nbsteps, coarse_null) ; 
steps.resize(nbsteps, stepnull) ; 

for (i=loop[0] ; i<loop[2] ; i+=loop[1])
  {
  actions.valeur=i ;

  if (type==0 || type==2) dump.check_timestep(i) ; 
  if (type==1) dumpatm.check_timestep(i) ; // Normalement, si le dumpatm n'est pas seul il est écrit lors de l'appel au dumpCF
  
  steps[(i-loop[0])/loop[1]].Type=TCOARSE ;    
  coarse[(i-loop[0])/loop[1]].cstep=&(steps[(i-loop[0])/loop[1]]) ;

  if(type==1) coarse[(i-loop[0])/loop[1]].coarse(dumpatm.steps[i]) ; 
  else        coarse[(i-loop[0])/loop[1]].coarse(dump.steps[i]) ; 
 
  switch(type) {
       case 0 : coarse[(i-loop[0])/loop[1]].do_coarse(dump.steps[i], stepnull, type) ; break ; 
       case 1 : coarse[(i-loop[0])/loop[1]].do_coarse(stepnull, dumpatm.steps[i], type) ; break ; 
       case 2 : coarse[(i-loop[0])/loop[1]].do_coarse(dump.steps[i], dumpatm.steps[i], type) ; break ; 
       }
  }

return 1 ;
}

//-----------------------------------------------
int CoarseDump::mean(int meantype)
{
Step nullstep ; 
Coarse & tmp = coarse[0] ; int i, j, k ; 
int idx[6] ;
double ** stdp, **stdt, **stdphi ; 
vector <int> weight ; 
vector <vector <double> > total ; 

// Pour calculer la standard deviation en prenant en compte les éventuelles correlations
stdp=(double**)malloc(sizeof(double*)*tmp.nb_boites_tot) ; 
for (i=0 ; i<tmp.nb_boites_tot ; i++) stdp[i]=(double*)malloc(sizeof(double)*nbsteps) ; 
stdt=(double**)malloc(sizeof(double*)*tmp.nb_boites_tot) ; 
for (i=0 ; i<tmp.nb_boites_tot ; i++) stdt[i]=(double*)malloc(sizeof(double)*nbsteps) ; 
stdphi=(double**)malloc(sizeof(double*)*tmp.nb_boites_tot) ; 
for (i=0 ; i<tmp.nb_boites_tot ; i++) stdphi[i]=(double*)malloc(sizeof(double)*nbsteps) ;
if ( !actions["is2D"].set) DISP_Warn("COARSTDTAU n'est calculé que pour sigmaxy sigmayx") ; 
idx[0]=tmp.cstep->find_idx(IDS("COARSIGTOTXX")) ; idx[1]=tmp.cstep->find_idx(IDS("COARSIGTOTYY")) ; idx[2]=tmp.cstep->find_idx(IDS("COARSIGTOTZZ")) ; 
idx[3]=tmp.cstep->find_idx(IDS("COARSIGTOTXY")) ; idx[4]=tmp.cstep->find_idx(IDS("COARSIGTOTYX")) ; idx[5]=tmp.cstep->find_idx(IDS("COARPHI")) ; 

if (meantype==1)
{
 weight.resize(tmp.cstep->nb_idx, 0) ; 
 total.resize(tmp.cstep->nb_idx) ; 
 for (i=0 ; i<tmp.cstep->nb_idx ; i++)
   total[i].resize(tmp.nb_boites_tot, 0) ; 
 
 for (i=0 ; i<tmp.cstep->nb_idx ; i++)
 {
   if (tmp.cstep->idx_col[i]==IDS("ID") || 
       tmp.cstep->idx_col[i]==IDS("POSX") || tmp.cstep->idx_col[i]==IDS("POSY") || tmp.cstep->idx_col[i]==IDS("POSZ"))
   { weight[i]=-2 ; for (k=0 ; k<tmp.nb_boites_tot ; k++) total[i][k]=1 ;}
   else if (tmp.cstep->idx_col[i]==IDS("VX") || tmp.cstep->idx_col[i]==IDS("VY") || tmp.cstep->idx_col[i]==IDS("VZ") ||
            tmp.cstep->idx_col[i]==IDS("COARPHI") || tmp.cstep->idx_col[i]==IDS("COARRAD") )
   {
     weight[i]=tmp.cstep->find_idx(IDS("COARATM")) ; 
     for (k=0 ; k<tmp.nb_boites_tot ; k++) 
        {
          tmp.cstep->datas[i][k]*=tmp.cstep->datas[weight[i]][k] ; 
          total[i][k] += tmp.cstep->datas[weight[i]][k] ;
        }
   }
   else if (tmp.cstep->idx_col[i]==IDS("COARATM") || tmp.cstep->idx_col[i]==IDS("COARCONTACTS") ||
            tmp.cstep->idx_col[i]==IDS("COARNBGP")||
            tmp.cstep->idx_col[i]==IDS("CFBREAKCOUNT") || tmp.cstep->idx_col[i]==IDS("CFFMAX") || tmp.cstep->idx_col[i]==IDS("CFKSPR") || tmp.cstep->idx_col[i]==IDS("CFLSPR")
           )
   {
     weight[i]=-1 ;
     for (k=0 ; k<tmp.nb_boites_tot ; k++) 
        total[i][k]+=1 ; 
   }
   else if ( tmp.cstep->idx_col[i]==IDS("SIGMAXX") || tmp.cstep->idx_col[i]==IDS("SIGMAXY") || tmp.cstep->idx_col[i]==IDS("SIGMAXZ") ||
             tmp.cstep->idx_col[i]==IDS("SIGMAYX") || tmp.cstep->idx_col[i]==IDS("SIGMAYY") || tmp.cstep->idx_col[i]==IDS("SIGMAYZ") ||
             tmp.cstep->idx_col[i]==IDS("SIGMAZX") || tmp.cstep->idx_col[i]==IDS("SIGMAZY") || tmp.cstep->idx_col[i]==IDS("SIGMAZZ") ||
             tmp.cstep->idx_col[i]==IDS("COARSIGKXX") || tmp.cstep->idx_col[i]==IDS("COARSIGKXY") || tmp.cstep->idx_col[i]==IDS("COARSIGKXZ") ||
             tmp.cstep->idx_col[i]==IDS("COARSIGKYX") || tmp.cstep->idx_col[i]==IDS("COARSIGKYY") || tmp.cstep->idx_col[i]==IDS("COARSIGKYZ") ||
             tmp.cstep->idx_col[i]==IDS("COARSIGKZX") || tmp.cstep->idx_col[i]==IDS("COARSIGKZY") || tmp.cstep->idx_col[i]==IDS("COARSIGKZZ") ||
             tmp.cstep->idx_col[i]==IDS("COARSIGTOTXX") || tmp.cstep->idx_col[i]==IDS("COARSIGTOTXY") || tmp.cstep->idx_col[i]==IDS("COARSIGTOTXZ") ||
             tmp.cstep->idx_col[i]==IDS("COARSIGTOTYX") || tmp.cstep->idx_col[i]==IDS("COARSIGTOTYY") || tmp.cstep->idx_col[i]==IDS("COARSIGTOTYZ") ||
             tmp.cstep->idx_col[i]==IDS("COARSIGTOTZX") || tmp.cstep->idx_col[i]==IDS("COARSIGTOTZY") || tmp.cstep->idx_col[i]==IDS("COARSIGTOTZZ") )
   { weight[i]=tmp.cstep->find_idx(IDS("COARCONTACTS")) ;
     for (k=0 ; k<tmp.nb_boites_tot ; k++) 
        {
          tmp.cstep->datas[i][k] *= tmp.cstep->datas[weight[i]][k] ; 
          total[i][k] += tmp.cstep->datas[weight[i]][k] ;
        }
   }
   else if (tmp.cstep->idx_col[i]==IDS("COARORIENTXX") || tmp.cstep->idx_col[i]==IDS("COARORIENTXY") || tmp.cstep->idx_col[i]==IDS("COARORIENTXZ") ||
            tmp.cstep->idx_col[i]==IDS("COARORIENTYX") || tmp.cstep->idx_col[i]==IDS("COARORIENTYY") || tmp.cstep->idx_col[i]==IDS("COARORIENTYZ") ||
            tmp.cstep->idx_col[i]==IDS("COARORIENTZX") || tmp.cstep->idx_col[i]==IDS("COARORIENTZY") || tmp.cstep->idx_col[i]==IDS("COARORIENTZZ") )
   { weight[i]=tmp.cstep->find_idx(IDS("COARNBGP")) ;  
     for (k=0 ; k<tmp.nb_boites_tot ; k++) 
        {
          tmp.cstep->datas[i][k] *= tmp.cstep->datas[weight[i]][k] ; 
          total[i][k] += tmp.cstep->datas[weight[i]][k] ;
        }
   }
 }
}

for (j=0 ; j<tmp.nb_boites_tot ; j++)
   {
     stdp[j][0]=(coarse[0].cstep->datas[idx[0]][j]+coarse[0].cstep->datas[idx[1]][j])/2. ; 
     stdt[j][0]=(coarse[0].cstep->datas[idx[3]][j]+coarse[0].cstep->datas[idx[4]][j])/2. ;
     stdphi[j][0]=coarse[0].cstep->datas[idx[5]][j] ;
   }

for (i=1 ; i<nbsteps ; i++)
  {
    if (meantype==0)
      tmp=tmp+coarse[i] ; // non-weighted mean
    else
    {
     for (j=0 ; j<tmp.cstep->nb_idx ; j++)
     {
      if (weight[j]==-2) // These one should not be averaged
      {
        for (k=0 ; k<tmp.nb_boites_tot ; k++)
          total[j][k]=1 ; 
      }
      else if (weight[j]==-1) // These one should be averaged just with the # of timestep
      {
        for (k=0 ; k<tmp.nb_boites_tot ; k++)
        {
          tmp.cstep->datas[j][k] += coarse[i].cstep->datas[j][k] ; 
          total[j][k]+=1 ; 
        }
      }
      else // Thes one are weighted average
      {
        for (k=0 ; k<tmp.nb_boites_tot ; k++) 
        {
          tmp.cstep->datas[j][k] += coarse[i].cstep->datas[weight[j]][k]*coarse[i].cstep->datas[j][k] ; 
          total[j][k] += coarse[i].cstep->datas[weight[j]][k] ;
        }
      }
     }
    }
   
   for (j=0 ; j<tmp.nb_boites_tot ; j++)
   {
     stdp[j][i]=(coarse[i].cstep->datas[idx[0]][j]+coarse[i].cstep->datas[idx[1]][j])/2. ; 
     stdt[j][i]=(coarse[i].cstep->datas[idx[3]][j]+coarse[i].cstep->datas[idx[4]][j])/2. ;
     stdphi[j][i]=coarse[i].cstep->datas[idx[5]][j] ;
   }
  }
for (j=0 ; j<tmp.cstep->nb_idx ; j++)
  {
   if (meantype==0)
   {
    if (tmp.cstep->idx_col[j]==IDS("ID") || 
        tmp.cstep->idx_col[j]==IDS("POSX") || tmp.cstep->idx_col[j]==IDS("POSY") || tmp.cstep->idx_col[j]==IDS("POSZ")) continue;
    for (i=0 ; i<tmp.nb_boites_tot ; i++) 
    {
      tmp.cstep->datas[j][i]/=nbsteps;
    }
   }
   else
   {
    for (i=0 ; i<tmp.nb_boites_tot ; i++) 
    {
      tmp.cstep->datas[j][i]/=total[j][i];
    }
   }
  }

printf("Averaged %d steps\n", nbsteps) ; 

tmp.cstep->idx_col.push_back(IDS("COARSTDPRES")) ; tmp.cstep->idx_col.push_back(IDS("COARSTDTAU")) ; tmp.cstep->idx_col.push_back(IDS("COARSTDPHI")) ; 
tmp.cstep->nb_idx+=3 ; tmp.cstep->datas.resize(tmp.cstep->idx_col.size()) ; 
tmp.cstep->datas[tmp.cstep->nb_idx-3].resize(tmp.cstep->nb_atomes,0) ; tmp.cstep->datas[tmp.cstep->nb_idx-2].resize(tmp.cstep->nb_atomes,0) ; tmp.cstep->datas[tmp.cstep->nb_idx-1].resize(tmp.cstep->nb_atomes,0) ; 
for (i=0 ; i<tmp.nb_boites_tot ; i++) 
  {
    tmp.cstep->datas[tmp.cstep->nb_idx-3][i]=Calcul::blockstd(stdp[i],nbsteps ) ;
    tmp.cstep->datas[tmp.cstep->nb_idx-2][i]=Calcul::blockstd(stdt[i],nbsteps ) ;    
    tmp.cstep->datas[tmp.cstep->nb_idx-1][i]=Calcul::blockstd(stdphi[i],nbsteps ) ;
  }
	  
nbsteps=1 ;

for (i=0 ; i<tmp.nb_boites_tot ; i++) 
{free(stdp[i]) ; free(stdt[i]) ; free(stdphi[i]) ; }
free(stdp) ; free(stdt) ; free(stdphi) ;  

return 1 ; 
}
//------------------------------------------------
int CoarseDump::write(string chemin)
{
long int loop[3] ;
int ts, i; 
Writing * W ;
bool stop ; 
int idx[9] ; 

W=&(actions.ecrire) ; 
W->chemin=&(chemin) ; 

pthread_mutex_init(&W->mutex, NULL);
pthread_cond_init (&W->sigin, NULL);
pthread_cond_init (&W->sigout, NULL);
pthread_mutex_lock(&W->mutex);
W->startall() ; 

stop=false ; 
loop[2]=nbsteps ; loop[0]=0 ; loop[1]=1 ; 
cout << "Coarse::write       " ;
open("") ;
actions.set_progress(loop) ; actions.disp_progress() ;

while (stop==false) 
{
  pthread_cond_wait(&W->sigout, &W->mutex); 
  switch (W->Signal)
  {
    case WHATAREYOU : W->sendin(COARSEDUMP) ; break ;
    case ASKINGGRID : W->i=coarse[ts].nb_boites ; W->sendin(OK) ; break ; 
    case ASKINGTS : W->sendin(UNABLE) ; break ;
    case FIRSTTS : ts=loop[0] ; W->i=&ts ; check_timestep(ts) ; fflush(stdout) ;  W->sendin(OK) ; break ;
    case CHOOSELDUMP : W->sendin(OK) ; break ; 
    case CHOOSECFDUMP : W->sendin(UNABLE) ;  break ; 
    case NEXTTS : 
      ts+=loop[1] ; actions.valeur=ts ; 
      if (ts>=loop[2]) 
	{ts=-1 ; W->i=&ts ; W->sendin(FINISH) ; } 
      else 
	{check_timestep(ts) ; W->i=&ts ; W->sendin(OK) ; }
      break ; 
    case ASKING1D : 
      idx[0]=steps[ts].find_idx(*(W->i)) ; 
      if (idx[0]==-1) {DISP_Warn("Donnée pas trouvée 1D\n") ; W->sendin(OK) ;}
      else {W->d[0]=&(steps[ts].datas[idx[0]][0]) ; W->i=&steps[ts].nb_atomes ; W->sendin(OK) ; }
      break ;  
    case ASKING2D : 
      for (i=0 ; i<3 ; i++)
      {
	idx[i]=steps[ts].find_idx(*(W->i+i)) ; 
	if(idx[i]==-1) {DISP_Warn("Donnée pas trouvée 2D: ") ; printf("%d \n",*(W->i+i)) ;  W->sendin(OK) ;break ; }
      }
      if (i<3) break ; 
      for (i=0 ; i<3 ; i++)
        W->d[i]=&(steps[ts].datas[idx[i]][0]) ; 
      W->i=&steps[ts].nb_atomes ; W->sendin(OK) ; break ;  
    case ASKING3D : 
      for (i=0 ; i<9 ; i++) 
      {
	idx[i]=steps[ts].find_idx(*(W->i+i)) ;
	if(idx[i]==-1) {DISP_Warn("Donnée pas trouvée 3D\n") ; W->sendin(OK) ;break ; }
      }
      if (i<9) break ; 
      for (i=0 ; i<9 ; i++) 
	W->d[i]=&(steps[ts].datas[idx[i]][0]) ; 
      W->i=&steps[ts].nb_atomes ; W->sendin(OK) ; break ;  
    case ASKINGND : 
      for (i=0 ; *(W->i+i)!=-1 ; i++) 
	W->d[i]=&(steps[ts].datas[steps[ts].find_idx(*(W->i+i))][0]) ;
      W->i=&steps[ts].nb_atomes ; W->sendin(OK) ; break ;    
    case FINI : stop=true ; W->sendin(OK) ; break ; 
    default : DISP_Warn("Signal d'écriture inconnu"); W->sendin(UNABLE) ; W->disp_Signal() ; break ; 
  }
}
return 1 ;   
}
//============================================================
//============================================================
Coarse::Coarse()
{
 nb_boites_tot=0 ; 
}
//------------------------------------------------------------
int Coarse::coarse (Step &step)
{
int j, k, l, m, v; Vector tailles ; 
bool specatm=false, speccf=false ; 
static bool infofenetre=true ;
static bool showwarn=true, showwarn2=true ;
static unsigned char infobox=0 ;

if (actions["use-box"].set)
    { 
	if (actions["use-box"]["box_xmin"]<step.box[0][0]) {infobox|=32 ; borders[0][0]=step.box[0][0]; } else {borders[0][0]=actions["use-box"]["box_xmin"] ;}
	if (actions["use-box"]["box_xmax"]>step.box[0][1]) {infobox|=16 ; borders[0][1]=step.box[0][1]; } else {borders[0][1]=actions["use-box"]["box_xmax"] ;}
	if (actions["use-box"]["box_ymin"]<step.box[1][0]) {infobox|=8  ; borders[1][0]=step.box[1][0]; } else {borders[1][0]=actions["use-box"]["box_ymin"] ;}
	if (actions["use-box"]["box_ymax"]>step.box[1][1]) {infobox|=4  ; borders[1][1]=step.box[1][1]; } else {borders[1][1]=actions["use-box"]["box_ymax"] ;}
	if (actions["use-box"]["box_zmin"]<step.box[2][0]) {infobox|=2  ; borders[2][0]=step.box[2][0]; } else {borders[2][0]=actions["use-box"]["box_zmin"] ;}
	if (actions["use-box"]["box_zmax"]>step.box[2][1]) {infobox|=1  ; borders[2][1]=step.box[2][1]; } else {borders[2][1]=actions["use-box"]["box_zmax"] ;}

	if ((infobox > 0)  && (!(infobox & 128))) {DISP_Info("\nUse-box plus grand que la simu : ") ; printf("0x%x\n", infobox) ;  infobox|=128 ; }

    }
else
    { 
    borders[0][0]=step.box[0][0] ; borders[0][1]=step.box[0][1] ; borders[1][0]=step.box[1][0] ; 
    borders[1][1]=step.box[1][1] ; borders[2][0]=step.box[2][0] ; borders[2][1]=step.box[2][1] ;
    }

nb_boites[0]=actions["coarse-graining"]["nbbox_x"] ; 
nb_boites[1]=actions["coarse-graining"]["nbbox_y"] ;
nb_boites[2]=actions["coarse-graining"]["nbbox_z"] ;
nb_boites_tot=nb_boites[0]*nb_boites[1]*nb_boites[2] ; 
cstep->nb_atomes=nb_boites_tot ;
VolumeO.resize(nb_boites_tot, 0) ; 
Volume.resize(nb_boites_tot, 0) ; 
// Ecriture des données de boîtes
// Initialisation du fenêtrage ...

// Par default : fenêtre0 est winsphere (moyennage des atomes)
//               fenetre est creneau3D                 

if (actions["wingauss"].set)
   {
   DISP_Warn("Ne devrait pas être utilisé. Par defaut, winsphere pour les coarse ATM, creneau3D pour les coarse CF\n") ;   
   specatm=speccf=true ; 
   fenetre=&Fonction::int_gaussienne3D ; 
   fenetreO=&Fonction::gaussienne3D ; 
   DISP_Warn("Volume0 et Volume ne sont pas calculé pour Wingauss\n") ; 
   if (actions["wingauss"]["sigma"]==-1)
      {
      DISP_Warn("Attention : le calcul automatique du sigma pour une wingauss simple ne se fait qu'à partir du nombre de boîte dans la dimension x\n") ;
      if (actions["sig-boites"].set) {tailles(1)=actions["sig-boites"]["sig-x"] ; tailles(2)=actions["sig-boites"]["sig-y"] ; tailles(3)=actions["sig-boites"]["sig-z"] ;}
      else tailles=taille_boite(step) ; 
      sigmafenetre(1)=sigmafenetre(2)=sigmafenetre(3)=1/sqrt(2)*tailles(1) ;
      }
   else
      {sigmafenetre(1)=sigmafenetre(2)=sigmafenetre(3)=actions["wingauss"]["sigma"] ; }
   sigmafenetreO=sigmafenetre ;
   }
else if (actions["wingauss3D"].set)
   { DISP_Err("WINGAUSS3D not implemented") ; } 
else if (actions["wincreneau"].set)
   {    
    specatm=true ; 
    if (actions["sig-boites"].set) {tailles(1)=actions["sig-boites"]["sig-x"] ; tailles(2)=actions["sig-boites"]["sig-y"] ; tailles(3)=actions["sig-boites"]["sig-z"] ;}
    else sigmafenetreO=taille_boite(step) ;
    if (actions["is2D"].set)
    {
      fenetreO=&Fonction::intersect_creneau2D ; 
      for (v=0 ; v<VolumeO.size() ; v++) VolumeO[v]=sigmafenetreO(1)*sigmafenetreO(2) ; 
    }
    else
    {
      if (infofenetre)
         DISP_Warn("Win creneau pas clair en 3D...\n") ;   
      if (nb_boites[0]==1 && nb_boites[1]==1)
      {
	if (infofenetre)
	    DISP_Warn ("... mais comme il n'y a qu'une boite en x et y, on dit que c'est une tranche ! EN TEST\n") ; 
	fenetreO=&Fonction::intersect_tranche_z ; 
      }
      else if (nb_boites[0]==1 && nb_boites[2]==1)
      {
	if (infofenetre)
	   DISP_Warn ("... mais comme il n'y a qu'une boite en x et z, on dit que c'est une tranche ! EN TEST\n") ; 
	fenetreO=&Fonction::intersect_tranche_y ; 
      }
      else if (nb_boites[1]==1 && nb_boites[2]==1)
      {
	if (infofenetre)
	   DISP_Warn ("... mais comme il n'y a qu'une boite en y et z, on dit que c'est une tranche ! EN TEST\n") ; 
	fenetreO=&Fonction::intersect_tranche_x ; 
      }
      else
	fenetreO=&Fonction::creneau3D ; 
      for (v=0 ; v<VolumeO.size() ; v++) VolumeO[v]=sigmafenetreO(1)*sigmafenetreO(2)*sigmafenetreO(3) ; 
    }
   }   
   
// Default case   
if (specatm==false)
   {
   if (infofenetre) {
	DISP_Info("Attention : pour le coarsening sur une fenêtre sphérique, le diamètre de la sphère est la plus petite dimension du pas spatial.\n") ;
	}
   if (actions["sig-boites"].set) {tailles(1)=actions["sig-boites"]["sig-x"] ; tailles(2)=actions["sig-boites"]["sig-y"] ; tailles(3)=actions["sig-boites"]["sig-z"] ;}
   else tailles=taille_boite(step) ;
   tailles(1)=tailles(1)>tailles(2)?tailles(2):tailles(1) ;
   if (!actions["is2D"].set) tailles(1)=tailles(1)>tailles(3)?tailles(3):tailles(1) ;
   sigmafenetreO(1)=sigmafenetreO(2)=sigmafenetreO(3)=tailles(1)/2.0 ;

   if (infofenetre) {
        DISP_Info("Attention : pour le coarsening sur une fenêtre sphérique, le diamètre de la sphère est la plus petite dimension du pas spatial : ") ;
	cout << sigmafenetreO(1)*1000*2 ; DISP_Info("mm.\n") ; }
   if (showwarn2 && sigmafenetreO(1)<actions.Cst["Radius"]) {DISP_Warn("Attention : la taille de la sphère de moyennage est plus petite que le rayon des atomes\n") ; showwarn2=false ;}

   //fenetre=&Fonction::int_intersect_sphere ;
   if (actions["is2D"].set) {fenetreO=&Fonction::intersect_cercle ; 
     for (int v=0 ; v<VolumeO.size() ; v++) VolumeO[v]=M_PI*sigmafenetreO(1)*sigmafenetreO(1) ; }
   else {fenetreO=&Fonction::intersect_sphere ; 
     for (int v=0 ; v<VolumeO.size() ; v++) VolumeO[v]=4/3.0*M_PI*sigmafenetreO(1)*sigmafenetreO(1)*sigmafenetreO(1) ; }
   
   specatm=true ;  
   }  
if (speccf==false)
   {
   // Fonction créneau 3D
   if (actions["is2D"].set) { if (infofenetre){DISP_Info("Chainforce fenetre creneau 2D\n") ;} fenetre = &Fonction::int_creneau2D ; }
   else {fenetre=&Fonction::int_creneau3D ; }

   if (actions["sig-boites"].set) {sigmafenetre(1)=actions["sig-boites"]["sig-x"] ; sigmafenetre(2)=actions["sig-boites"]["sig-y"] ; sigmafenetre(3)=actions["sig-boites"]["sig-z"] ;}
   else sigmafenetre=taille_boite(step) ;

   if (actions["winboxyper"].set)
   	   { if (nb_boites[1]>1 || (actions["use-box"].set && ((borders[1][1]-borders[1][0])<(step.box[1][1]-step.box[1][0]))))
   	     { DISP_Warn("A VERIFIER !!!!! Attention, pour une boite y periodique, la taille en y devrait être égale à l'épaisseur de la boite de simu, et il ne devrait y avoir qu'une seulle boite en y\n") ;
   	     }
   	   }
   if (showwarn && (sigmafenetre(1)<actions.Cst["Radius"] || sigmafenetre(2)<actions.Cst["Radius"] || (sigmafenetre(3)<actions.Cst["Radius"] && (!actions["is2D"].set))))
	   {DISP_Warn("Attention, l'une des dimension de la boîte de moyennage est inéfrieur au rayon des atomes"); showwarn=false ; }
	   
   if (!actions["is2D"].set) 
     for (int v=0 ; v<Volume.size() ; v++) Volume[v]=sigmafenetre(1)*sigmafenetre(2)*sigmafenetre(3) ; 
   else
     for (int v=0 ; v<Volume.size() ; v++) Volume[v]=sigmafenetre(1)*sigmafenetre(2) ;
   speccf=true ; 	   
   }   

// Inscrivons les boîtes de moyennage dans datas

cstep->idx_col.push_back(IDS("ID")) ;
cstep->idx_col.push_back(IDS("POSX")) ; cstep->idx_col.push_back(IDS("POSY")) ; cstep->idx_col.push_back(IDS("POSZ")) ;
cstep->idx_col.push_back(IDS("VX")) ; cstep->idx_col.push_back(IDS("VY")) ; cstep->idx_col.push_back(IDS("VZ")) ;
cstep->idx_col.push_back(IDS("COARPHI")); cstep->idx_col.push_back(IDS("COARATM")); cstep->idx_col.push_back(IDS("COARCONTACTS")); 
cstep->idx_col.push_back(IDS("COARRAD") );
cstep->idx_col.push_back(IDS("SIGMAXX")); cstep->idx_col.push_back(IDS("SIGMAXY")); cstep->idx_col.push_back(IDS("SIGMAXZ"));
cstep->idx_col.push_back(IDS("SIGMAYX")); cstep->idx_col.push_back(IDS("SIGMAYY")); cstep->idx_col.push_back(IDS("SIGMAYZ"));
cstep->idx_col.push_back(IDS("SIGMAZX")); cstep->idx_col.push_back(IDS("SIGMAZY")); cstep->idx_col.push_back(IDS("SIGMAZZ"));

cstep->idx_col.push_back(IDS("COARSIGKXX")); cstep->idx_col.push_back(IDS("COARSIGKXY")); cstep->idx_col.push_back(IDS("COARSIGKXZ"));
cstep->idx_col.push_back(IDS("COARSIGKYX")); cstep->idx_col.push_back(IDS("COARSIGKYY")); cstep->idx_col.push_back(IDS("COARSIGKYZ"));
cstep->idx_col.push_back(IDS("COARSIGKZX")); cstep->idx_col.push_back(IDS("COARSIGKZY")); cstep->idx_col.push_back(IDS("COARSIGKZZ"));

cstep->idx_col.push_back(IDS("COARSIGTOTXX") ); cstep->idx_col.push_back(IDS("COARSIGTOTXY") ); cstep->idx_col.push_back(IDS("COARSIGTOTXZ") );
cstep->idx_col.push_back(IDS("COARSIGTOTYX") ); cstep->idx_col.push_back(IDS("COARSIGTOTYY") ); cstep->idx_col.push_back(IDS("COARSIGTOTYZ") );
cstep->idx_col.push_back(IDS("COARSIGTOTZX") ); cstep->idx_col.push_back(IDS("COARSIGTOTZY") ); cstep->idx_col.push_back(IDS("COARSIGTOTZZ") );

if (actions["multisphere"].set || actions["multisphereflux"].set)
{
cstep->idx_col.push_back(IDS("COARNBGP"));
cstep->idx_col.push_back(IDS("COARORIENTXX") ); cstep->idx_col.push_back(IDS("COARORIENTXY") ); cstep->idx_col.push_back(IDS("COARORIENTXZ") );
cstep->idx_col.push_back(IDS("COARORIENTYX") ); cstep->idx_col.push_back(IDS("COARORIENTYY") ); cstep->idx_col.push_back(IDS("COARORIENTYZ") );
cstep->idx_col.push_back(IDS("COARORIENTZX") ); cstep->idx_col.push_back(IDS("COARORIENTZY") ); cstep->idx_col.push_back(IDS("COARORIENTZZ") );
}
if (IDS.exist("CFBREAKCOUNT")) cstep->idx_col.push_back(IDS("CFBREAKCOUNT") ); //steps[timestep].idx_col.push_back(IDS("CFFMAX") );
if (IDS.exist("CFFMAX")) cstep->idx_col.push_back(IDS("CFFMAX") );
if (IDS.exist("CFKSPR")) cstep->idx_col.push_back(IDS("CFKSPR") );
if (IDS.exist("CFLSPR")) cstep->idx_col.push_back(IDS("CFLSPR") );

cstep->nb_idx=cstep->idx_col.size() ;
cstep->datas.resize(cstep->idx_col.size()) ; 
for (j=0 ; j<cstep->idx_col.size() ; j++) cstep->datas[j].resize(cstep->nb_atomes,0) ; 
 
for (j=0 ; j<nb_boites[0] ; j++)
  {
    for (k=0 ; k<nb_boites[1] ; k++)
    {
     for (l=0 ; l<nb_boites[2] ; l++)
      {
      cstep->datas[0][linidx(j,k,l)]=j+k*nb_boites[0]+l*nb_boites[0]*nb_boites[1] ; 
	
      cstep->datas[1][linidx(j,k,l)]=(j+0.5)*(borders[0][1]-borders[0][0])/(double)nb_boites[0]+borders[0][0] ; 
      cstep->datas[2][linidx(j,k,l)]=(k+0.5)*(borders[1][1]-borders[1][0])/(double)nb_boites[1]+borders[1][0] ; 
      cstep->datas[3][linidx(j,k,l)]=(l+0.5)*(borders[2][1]-borders[2][0])/(double)nb_boites[2]+borders[2][0] ; 

      for (m=4 ; m<cstep->datas.size() ; m++)
         cstep->datas[m][linidx(j,k,l)]=0 ;
      }
    }
  }  

// Compute the volume in the case of delcyl et deltank
// Very basic programmation, I don't want to deal with that right now. 
if (actions["delcyl"].set)
{
  Vector xcyl (actions["delcyl"]["xcyl"], actions["delcyl"]["ycyl"], actions["delcyl"]["zcyl"]) ; 
  Vector dircyl (actions["delcyl"]["dirxcyl"], actions["delcyl"]["dirycyl"], actions["delcyl"]["dirzcyl"]) ; 
  Vector xsphere ;
  if (infofenetre) DISP_Warn ("delcyl est en phase de test. Fonctionne seulement en winsphere !!") ; 
  for (v=0 ; v<nb_boites_tot ; v++)
  {
   xsphere(1)=cstep->datas[1][v] ; xsphere(2)=cstep->datas[2][v] ; xsphere(3)=cstep->datas[3][v] ; 
   VolumeO[v]=Fonction::volume_intersect_sphere_cylinder (xsphere, sigmafenetreO(1), xcyl, dircyl, actions["delcyl"]["rayon"], actions["delcyl"]["longueur"]); 
  }
}
if (actions["deltank"].set) // VERY VERY CRUDE
{  
  DISP_Warn("Attention, l'utilisation de deltank est très simplifiée (pas de longueur du tank)\n") ; 
  Vector pos ; 
  for (v=0 ; v<nb_boites_tot ; v++)
  {
   pos(1)=cstep->datas[1][v] ; pos(2)=cstep->datas[2][v] ; pos(3)=cstep->datas[3][v] ;
   Volume[v]=Fonction::volume_intersect_pave_cylindre(pos, sigmafenetre, actions["deltank"]["rayon"]) ;  
  }
}
  
infofenetre=false ;   
return 1 ;
}

//-------------------------------------------------
Vector Coarse::taille_boite (Step &step)
{
double borders[3][2] ;
Vector retour ; 

if (actions["use-box"].set)
    {
	if (actions["use-box"]["box_xmin"]<step.box[0][0]) {borders[0][0]=step.box[0][0]; } else {borders[0][0]=actions["use-box"]["box_xmin"] ;}
	if (actions["use-box"]["box_xmax"]>step.box[0][1]) {borders[0][1]=step.box[0][1]; } else {borders[0][1]=actions["use-box"]["box_xmax"] ;}
	if (actions["use-box"]["box_ymin"]<step.box[1][0]) {borders[1][0]=step.box[1][0]; } else {borders[1][0]=actions["use-box"]["box_ymin"] ;}
	if (actions["use-box"]["box_ymax"]>step.box[1][1]) {borders[1][1]=step.box[1][1]; } else {borders[1][1]=actions["use-box"]["box_ymax"] ;}
	if (actions["use-box"]["box_zmin"]<step.box[2][0]) {borders[2][0]=step.box[2][0]; } else {borders[2][0]=actions["use-box"]["box_zmin"] ;}
	if (actions["use-box"]["box_zmax"]>step.box[2][1]) {borders[2][1]=step.box[2][1]; } else {borders[2][1]=actions["use-box"]["box_zmax"] ;}
    }
else
    {
    borders[0][0]=step.box[0][0] ; borders[0][1]=step.box[0][1] ; borders[1][0]=step.box[1][0] ; 
    borders[1][1]=step.box[1][1] ; borders[2][0]=step.box[2][0] ; borders[2][1]=step.box[2][1] ;
    }

(retour)(1)=(borders[0][1]-borders[0][0])/(double)actions["coarse-graining"]["nbbox_x"] ;
(retour)(2)=(borders[1][1]-borders[1][0])/(double)actions["coarse-graining"]["nbbox_y"] ;
(retour)(3)=(borders[2][1]-borders[2][0])/(double)actions["coarse-graining"]["nbbox_z"] ;
return (retour) ; 
}

//==========================================================================
Coarse operator+(Coarse c1, Coarse c2)
{
  int i, j ; Coarse temp ;  
 
 if (c1==c2)
   temp=c1 ;  
 else
   {
   cout << "ERR : tentative d'addition de deux coarse n'ayant pas des propriétés identiques" ;
   return temp ; 
   }
 
 for (j=0 ; j<temp.cstep->nb_idx ; j++)
 {
   if (temp.cstep->idx_col[j]==IDS("ID") || temp.cstep->idx_col[j]==IDS("POSX") || temp.cstep->idx_col[j]==IDS("POSY") || temp.cstep->idx_col[j]==IDS("POSZ")) continue;
   for (i=0 ; i<temp.nb_boites_tot ; i++)
   {
   temp.cstep->datas[j][i]=c1.cstep->datas[j][i]+c2.cstep->datas[j][i] ; 
   }
 }
 
 return temp ; 
}

//---------------------------------------------------
void Coarse::operator=(Coarse c)
{
 nb_boites[0]=c.nb_boites[0] ; 
 nb_boites[1]=c.nb_boites[1] ;
 nb_boites[2]=c.nb_boites[2] ;
 nb_boites_tot=c.nb_boites_tot ;
 cstep=c.cstep ; 
 useCF=c.useCF ; 
}

//----------------------------------------------------
bool operator==(Coarse c1, Coarse c2)
{
  if (c1.nb_boites[0]==c2.nb_boites[0] && c1.nb_boites[1]==c2.nb_boites[1] && c1.nb_boites[2]==c2.nb_boites[2]
    && c1.borders[0][0]==c2.borders[0][0] && c1.borders[0][1]==c2.borders[0][1] && c1.borders[1][0]==c2.borders[1][0]
    && c1.borders[1][1]==c2.borders[1][1] && c1.borders[2][0]==c2.borders[2][0] && c1.borders[2][1]==c2.borders[2][1])
    return true ; 
  else
    return false ;   
}

//==========================================================================
int Coarse::do_coarse(Step &step, Step &stepatm, int type) // On fait les 2 coarse en même temps.
{
int i, j, k, l ; 
double width, height, depth, rayon, masse, densite ;
int idx[16], idxatm[8], idxlin[2] ; 
static bool infoonce=true ; 
int bool_coarse=0 ; bool domultisphere=false ; 
bool est2D ; 
Vector r, f, xa, xb ; 

// 1 : checks du type de coarse 0 cfcoarse seulement (cannot do w/kinetic), 1 atmcoarse seul, 2 both
if (type==1 && actions["w/kinetic"].set )
   {DISP_Err("ERR: impossible de calculer les contraintes cinétiques sans dump atomique");
    return 1 ; }
switch(type) {
  case 0 : if (step.Type!=TCF) {DISP_Err("ERR: impossible de calculer les contraintes. Types de step incorrects"); return 2 ;} break;
  case 1 : if (stepatm.Type!=TL) {DISP_Err("ERR: impossible de calculer les vitesses. Types de step incorrects"); return 2 ;} break ; 
  case 2 : if (step.Type!=TCF && stepatm.Type!=TL) {DISP_Err("ERR: impossible de calculer les contraintes. Types de step incorrects"); return 2 ;} break ;
  }
  
// 2 : récupération des indexes  
if (type==0 || type==2)
  {
  idx[0]=step.find_idx(IDS("CFPOSX"))     ; idx[1]=step.find_idx(IDS("CFPOSY"))     ; idx[2]=step.find_idx(IDS("CFPOSZ")) ;
  idx[3]=step.find_idx(IDS("CFX"))      ; idx[4]=step.find_idx(IDS("CFY"))      ; idx[5]=step.find_idx(IDS("CFZ")) ;
  idx[6]=step.find_idx(IDS("CFFORCEX")) ; idx[7]=step.find_idx(IDS("CFFORCEY")) ; idx[8]=step.find_idx(IDS("CFFORCEZ")) ;
  if (IDS.exist("CFBREAKCOUNT")) bool_coarse|=1 ; 
  if (IDS.exist("CFFMAX"))       bool_coarse|=2 ; 
  if (IDS.exist("CFKSPR"))       bool_coarse|=4 ; 
  if (IDS.exist("CFLSPR"))       bool_coarse|=8 ; 

  if (bool_coarse&1) idx[9]=step.find_idx(IDS("CFBREAKCOUNT")) ; /*idx[10]=step.find_idx(IDS("CFFMAX")) ;*/
  if (bool_coarse&2) idx[10]=step.find_idx(IDS("CFFMAX")) ;
  if (bool_coarse&4) idx[11]=step.find_idx(IDS("CFKSPR")) ;
  if (bool_coarse&8) idx[12]=step.find_idx(IDS("CFLSPR")) ;
    
  for (j=0 ; j<9 ; j++) { if (idx[j]==-1) { printf("%d ", j) ; DISP_Err ("ERR: donnees manquantes pour le coarse graining\n") ; return 3 ;}}
  }
if (type==1 || type==2)
  {
  idxatm[0]=stepatm.find_idx(IDS("POSX")) ; idxatm[1]=stepatm.find_idx(IDS("POSY"))  ; idxatm[2]=stepatm.find_idx(IDS("POSZ")) ;
  idxatm[3]=stepatm.find_idx(IDS("VX"))   ; idxatm[4]=stepatm.find_idx(IDS("VY"))    ; idxatm[5]=stepatm.find_idx(IDS("VZ")) ; 
  for (j=0 ; j<6 ; j++) { if (idxatm[j]==-1) { printf("/%d/ ", j) ;  DISP_Err ("ERR: donnees manquantes pour le coarse graining\n") ; printf("[%d %d]\n",step.timestep, stepatm.timestep) ;  return 3 ;}}
  
  idxatm[6]=stepatm.find_idx(IDS("RAYON")) ;
  if (idxatm[6]==-1)
    {
    if (infoonce)
    {
      DISP_Info("Pas de données de RAYON. Le rayon atm utilisé pour le coarse graining est celui de Actions ") ; 
      infoonce=false ; 
    }
    rayon=actions.Cst["Radius"] ; 
    }
  densite=actions.Cst["Rhograin"] ; 
  est2D=actions["is2D"].set ;   
  }
  
  
// 3 : coarse graining pour vitesses si nécessaire pour commencer
double coeff ; Vector position ; 
int natomescompt=stepatm.nb_atomes  ;
int bidx[3], idxclosest, ii, jj, kk, bornes[6] ;
Vector vitesse, vitesset ;
//int DeltaBoites=calc_deltaboites(type) ; // TODO Compute deltaboite

// On récupère les vectors de datas qu'on met dans des références pour simplifier
vector <double> & vx=cstep->datas[cstep->find_idx(IDS("VX"))] ; 
vector <double> & vy=cstep->datas[cstep->find_idx(IDS("VY"))] ;
vector <double> & vz=cstep->datas[cstep->find_idx(IDS("VZ"))] ; 
vector <double> & nb_atm=cstep->datas[cstep->find_idx(IDS("COARATM"))] ;
vector <double> & mean_radius=cstep->datas[cstep->find_idx(IDS("COARRAD"))] ;
vector <double> & phi=cstep->datas[cstep->find_idx(IDS("COARPHI"))] ;
if (type==1 || type==2)
 {
 for (j=0 ; j<natomescompt ; j++)
  {
  if (isnan(stepatm.datas[idxatm[0]][j])) continue ; // L'atome a été rajouté mais n'existe pas.
  position(1)=stepatm.datas[idxatm[0]][j] ; position(2)=stepatm.datas[idxatm[1]][j] ; position(3)=stepatm.datas[idxatm[2]][j] ; 
  if (idxatm[6]!=-1) {rayon=stepatm.datas[idxatm[6]][j] ;}
  calc_bornes_v2(position(1), position(2), position(3), rayon, 1, bornes) ; //TODO
  for (ii=bornes[0] ; ii<=bornes[1] ; ii++)
   {
   for (jj=bornes[2] ; jj<=bornes[3] ; jj++)
    {
    for (kk=bornes[4] ; kk<=bornes[5] ; kk++)
     {
     idxclosest=linidx(ii, jj, kk) ;
     xa=position-position_vec(idxclosest) ; //boites[idxclosest].x0 ;
     coeff=VolumeO[idxclosest] ; 
     fenetreO(coeff, xa, sigmafenetreO, rayon) ;
     if (isnan(coeff)) 
     {DISP_Warn("Coeff a pris une valeur NAN") ; }
    
     if (est2D) masse=M_PI*densite*rayon*rayon ;
     else masse=4/3.*M_PI*densite*rayon*rayon*rayon ;
     
     nb_atm[idxclosest]+=coeff ;
     mean_radius[idxclosest]+=rayon*coeff ;
     phi[idxclosest]+=coeff*masse ; 
     vx[idxclosest]+=stepatm.datas[idxatm[3]][j]*coeff*masse ;
     vy[idxclosest]+=stepatm.datas[idxatm[4]][j]*coeff*masse ;
     vz[idxclosest]+=stepatm.datas[idxatm[5]][j]*coeff*masse ;

     }
    }
   }
  }
 for (j=0 ; j<nb_boites_tot ; j++)
  {
  phi[j]=phi[j]/densite ;
  if (nb_atm[j]>0)
     mean_radius[j]=mean_radius[j]/nb_atm[j] ;
  if (VolumeO[j]!=0) nb_atm[j]*=VolumeO[j] ; else DISP_Warn("Etrange d'avoir un volume nul ...") ; 
  if (phi[j]!=0)
    {
     vx[j]=vx[j]/(phi[j]*densite) ;
     vy[j]=vy[j]/(phi[j]*densite) ;
     vz[j]=vz[j]/(phi[j]*densite) ;
    }
  }
 } 

// 4 : calcul de la contrainte cinétique si besoin 
Matrix3x3 sigmak ; 
if ((type==1 || type==2) && actions["w/kinetic"].set)
 {
 for (j=0 ; j<natomescompt ; j++)
  {
  position(1)=stepatm.datas[idxatm[0]][j] ; position(2)=stepatm.datas[idxatm[1]][j] ; position(3)=stepatm.datas[idxatm[2]][j] ; 
  if (idxatm[6]!=-1) {rayon=stepatm.datas[idxatm[6]][j] ;}
  calc_bornes_v2(position(1), position(2), position(3), rayon, 1, bornes) ; 
  //printf("%d %d %d %d %d %d \n", bornes[0], bornes[1], bornes[2],bornes[3] ,bornes[4] ,bornes[5]) ; 
  for (ii=bornes[0] ; ii<=bornes[1] ; ii++)
   {
   for (jj=bornes[2] ; jj<=bornes[3] ; jj++)
    {
    for (kk=bornes[4] ; kk<=bornes[5] ; kk++)
     {
     idxclosest=linidx(ii, jj, kk) ;
     xa=position-position_vec(idxclosest) ;

     coeff=VolumeO[idxclosest] ; 
     fenetreO(coeff, xa, sigmafenetreO, rayon) ;
     if (isnan(coeff)) DISP_Warn("Coeff a pris une valeur NAN") ; 
     
     if (est2D) masse=M_PI*densite*rayon*rayon ;
     else masse=4/3.*M_PI*densite*rayon*rayon*rayon ;

     vitesse(1)=stepatm.datas[idxatm[3]][j]-vx[idxclosest] ; 
     vitesse(2)=stepatm.datas[idxatm[4]][j]-vy[idxclosest] ; 
     vitesse(3)=stepatm.datas[idxatm[5]][j]-vz[idxclosest] ;
     vitesset=vitesse.t() ;

     sigmak=(vitesse.multiply(vitesset))*coeff*masse*(+1) ; // NB : le signe est inversé pour être cohérent avec la contrainte dûe au contact, qui a aussi un signe inversé en raison d'une différence de convention dans les vecteurs
     addto_sigmak (idxclosest, sigmak) ; 
     }
    }
   }
  }
 }  
 
// 5 : coarse de forces si besoin
vector <double> & nb_contacts=cstep->datas[cstep->find_idx(IDS("COARCONTACTS"))] ;
vector <vector <double> * > cfother (10) ; 
if (bool_coarse&1) cfother[0]=&(cstep->datas[cstep->find_idx(IDS("CFBREAKCOUNT"))]) ;
if (bool_coarse&2) cfother[1]=&(cstep->datas[cstep->find_idx(IDS("CFFMAX"))]) ;
if (bool_coarse&4) cfother[2]=&(cstep->datas[cstep->find_idx(IDS("CFKSPR"))]) ;
if (bool_coarse&8) cfother[3]=&(cstep->datas[cstep->find_idx(IDS("CFLSPR"))]) ;
Matrix3x3 sigma ; 
natomescompt=step.nb_atomes  ;
f=f.t() ;  
if (type==0 || type==2)
{ 
 for (j=0 ; j<natomescompt ; j++)
  {
  position(1)=step.datas[idx[0]][j] ; position(2)=step.datas[idx[1]][j] ; position(3)=step.datas[idx[2]][j] ;  
  f(1)=step.datas[idx[6]][j] ; f(2)=step.datas[idx[7]][j] ; f(3)=step.datas[idx[8]][j] ;
  r(1)=step.datas[idx[3]][j] ; r(2)=step.datas[idx[4]][j] ; r(3)=step.datas[idx[5]][j] ;
  calc_bornes_v2(position(1), position(2), position(3), r.norm(), 0, bornes) ;
  for (ii=bornes[0] ; ii<=bornes[1] ; ii++)
   {
   for (jj=bornes[2] ; jj<=bornes[3] ; jj++)
    {
    for (kk=bornes[4] ; kk<=bornes[5] ; kk++)
     {
     idxclosest=linidx(ii, jj, kk) ;  
     xa=r/(-2.0)+position-position_vec(idxclosest) ; xb=r/2.+position-position_vec(idxclosest) ;
     coeff=Volume[idxclosest] ; 
     fenetre(coeff,xa, xb, sigmafenetre) ; 
     if (isnan(coeff)) DISP_Warn("Coeff a pris une valeur NAN\n") ; 
     sigma=(r.multiply(f))*coeff*(-1) ;
     addto_sigma (idxclosest, sigma) ;
     if (bool_coarse&1) (*(cfother[0]))[idxclosest]+=step.datas[idx[9]][j]*coeff ; //printf("%g ", (*cfbreak)[idxclosest]) ; 
     if (bool_coarse&2) (*(cfother[1]))[idxclosest]+=step.datas[idx[10]][j]*coeff ;
     if (bool_coarse&4) (*(cfother[2]))[idxclosest]+=step.datas[idx[11]][j]*coeff ;
     if (bool_coarse&8) (*(cfother[3]))[idxclosest]+=step.datas[idx[12]][j]*coeff ; 

     nb_contacts[idxclosest]+=coeff;
     }
    }
   }       
  } 
 if (bool_coarse>0)
 {
  for (j=0 ; j<nb_boites_tot ; j++)
   {
    if (nb_contacts[j]>0)
    {
     for (k=0,l=1 ; k<4 ; k++,l*=2)
       if (bool_coarse & l)
	 (*(cfother[k]))[j]/=nb_contacts[j] ; //printf("%g ", (*cfbreak)[idxclosest]) ; 
     //boites[j].cffmax/=boites[j].nb_contacts ;
    }
   }
   
  }
}  

// 6: compute coase in case of multisphere
if (actions["multisphere"].set || actions["multisphereflux"].set)
{
 int kxx=cstep->find_idx(IDS("COARORIENTXX")) ;
 vector <double> & nb_gp=cstep->datas[cstep->find_idx(IDS("COARNBGP"))] ;
 Vector3d Ksegment ; Matrix3d Kmatseg ; 
 int deb,end ; 
 
 stepatm.multisphere->get_orientations(stepatm) ; 

 if (actions["multisphereflux"].set) {deb=0 ; end=stepatm.multisphere->ngp ;}
 else {deb=1 ; end=stepatm.multisphere->ngp+1 ;}
 

 for (i=deb ; i<end ; i++) // NB: ngp start counting from 1 Fuck not always ......
 {
  if (stepatm.multisphere->data[0][i]!=GP_OK) {fflush(stdout) ; continue ;} // remove periodique, lost, etc.  

  position(1)=stepatm.multisphere->data[1][i] ; position(2)=stepatm.multisphere->data[2][i] ; position(3)=stepatm.multisphere->data[3][i] ; 

  calc_bornes_v0(position(1), position(2), position(3), bornes) ; 

  idxclosest=linidx(bornes[0], bornes[1], bornes[2]) ;
  Ksegment << stepatm.multisphere->data[4][i], stepatm.multisphere->data[5][i], stepatm.multisphere->data[6][i] ; 
  Ksegment=Ksegment/(Ksegment.norm()) ; 
  Kmatseg=Ksegment*(Ksegment.transpose());
  for (j=0 ; j<9 ; j++)
    {cstep->datas[kxx+j][idxclosest]+=Kmatseg(j/3,j%3) ; }
  nb_gp[idxclosest]++ ;  
 }
 for (j=0 ; j<nb_boites_tot ; j++)
 {
   if (nb_gp[j]>0)
   {
     for (k=0 ; k<9 ; k++) 
        cstep->datas[kxx+k][j]/=((double)(nb_gp[j])) ;
   }
 }
}

if (type==2) compute_sigmatot () ; 
return 1 ; 
}

//---------------------
void Coarse::addto_sigma (int idxclosest, Matrix3x3 Sig)
{
  static int idx=cstep->find_idx(IDS("SIGMAXX")) ; 
  int i ; 
  for (i=0 ; i<9 ; i++) {cstep->datas[idx+i][idxclosest]+=Sig[i] ;} 
}
void Coarse::addto_sigmak (int idxclosest, Matrix3x3 Sig)
{
  static int idx=cstep->find_idx(IDS("COARSIGKXX")) ; 
  int i ; 
  for (i=0 ; i<9 ; i++) {cstep->datas[idx+i][idxclosest]+=Sig[i] ; }
}
void Coarse::compute_sigmatot (void)
{
  static int idx0=cstep->find_idx(IDS("SIGMAXX")) , idx1=cstep->find_idx(IDS("COARSIGKXX")), idx2=cstep->find_idx(IDS("COARSIGTOTXX"))  ; 
  int i,j ; double check=0 ; 
  for (j=0 ; j<nb_boites_tot ; j++)
  {
   for (i=0 ; i<9 ; i++) 
   {
     if (isinf(cstep->datas[idx0+i][j]) || isnan(cstep->datas[idx0+i][j])) { cstep->datas[idx0+i][j]=0 ; }
     if (isinf(cstep->datas[idx1+i][j]) || isnan(cstep->datas[idx1+i][j])) { cstep->datas[idx1+i][j]=0 ; }
   }
   for (i=0 ; i<9 ; i++) {cstep->datas[idx2+i][j]=cstep->datas[idx0+i][j]+cstep->datas[idx1+i][j] ; //printf("{A%g}", cstep->datas[idx0+i][j]) ;
   //printf("{B%g}", cstep->datas[idx1+i][j]) ;
   check+=fabs(cstep->datas[idx1+i][j]) ; }
  }
  if (check<1e-100) DISP_Warn("La partie contrainte cinétique est très faible, il n'y a pas un oublie de w/kinetic ?\n") ; 
}
Vector Coarse::position_vec (int id)
{
  Vector res ;
  static int idx=cstep->find_idx(IDS("POSX")) ;
  res(1)=cstep->datas[idx][id] ;
  res(2)=cstep->datas[idx+1][id] ;
  res(3)=cstep->datas[idx+2][id] ;
  return res ;
}

//--------------------
void Coarse::calc_bornes_v2 (double x, double y, double z, double r, int curtype, int *bornes)
{
  int dx, dy, dz ; int bidx[3] ;
  static bool is2D=actions["is2D"].set ; 
  bidx[0]=int(floor((x-borders[0][0])/(double)((borders[0][1]-borders[0][0])/(double)nb_boites[0]))) ; 
  bidx[1]=int(floor((y-borders[1][0])/(double)((borders[1][1]-borders[1][0])/(double)nb_boites[1]))) ; 
  bidx[2]=int(floor((z-borders[2][0])/(double)((borders[2][1]-borders[2][0])/(double)nb_boites[2]))) ; 
  
  if (curtype==0)
  {
   dx=ceil(r/sigmafenetre(1)) ; 
   dy=ceil(r/sigmafenetre(2)) ; 
   if (!is2D) dz=ceil(r/sigmafenetre(3)) ;
   else dz=0 ; 
  }
  else if (curtype==1)
  {
   dx=ceil(r/sigmafenetreO(1));
   dy=ceil(r/sigmafenetreO(2));
   if (!is2D) dz=ceil(r/sigmafenetreO(3));
   else dz=0 ; 
  }
  else if (curtype==2) // not looking around
  {
    dx=dy=dz=0 ; 
  }
  else DISP_Err("Should never happen. Coarse::calc_bornes_v2") ; 
  
  if (dx<1 || dy<1 || (dz<1 && !is2D)) DISP_Warn("Coarse::calc_bornes_v2 very strange dx, dy or dz") ; 
  
  bornes[0]=((bidx[0]-dx)<0?0:(bidx[0]-dx)) ;
  bornes[2]=((bidx[1]-dy)<0?0:(bidx[1]-dy)) ;
  bornes[4]=((bidx[2]-dz)<0?0:(bidx[2]-dz)) ;
  bornes[1]=((bidx[0]+dx)>=nb_boites[0]?(nb_boites[0]-1):(bidx[0]+dx)) ;
  bornes[3]=((bidx[1]+dy)>=nb_boites[1]?(nb_boites[1]-1):(bidx[1]+dy)) ;
  bornes[5]=((bidx[2]+dz)>=nb_boites[2]?(nb_boites[2]-1):(bidx[2]+dz)) ;
  
}

void Coarse::calc_bornes_v0 (double x, double y, double z, int *bornes)
{
  int bidx[3] ;
  //static bool is2D=actions["is2D"].set ; 
  bidx[0]=floor((x-borders[0][0])/(double)((borders[0][1]-borders[0][0])/(double)nb_boites[0])) ; 
  bidx[1]=floor((y-borders[1][0])/(double)((borders[1][1]-borders[1][0])/(double)nb_boites[1])) ; 
  bidx[2]=floor((z-borders[2][0])/(double)((borders[2][1]-borders[2][0])/(double)nb_boites[2])) ; 
  
  
  bornes[0]=(bidx[0]<0?0:bidx[0]) ;
  bornes[1]=(bidx[1]<0?0:bidx[1]) ;
  bornes[2]=(bidx[2]<0?0:bidx[2]) ;
  bornes[0]=(bidx[0]>=nb_boites[0]?(nb_boites[0]-1):bidx[0]) ;
  bornes[1]=(bidx[1]>=nb_boites[1]?(nb_boites[1]-1):bidx[1]) ;
  bornes[2]=(bidx[2]>=nb_boites[2]?(nb_boites[2]-1):bidx[2]) ;
  
}

/*
void Coarse::calc_bornes (double x, double y, double z, int DeltaBoites, int * bornes)
{
int bidx[3] ; 

bidx[0]=floor((x-borders[0][0])/(double)((borders[0][1]-borders[0][0])/(double)nb_boites[0])) ; 
bidx[1]=floor((y-borders[1][0])/(double)((borders[1][1]-borders[1][0])/(double)nb_boites[1])) ; 
bidx[2]=floor((z-borders[2][0])/(double)((borders[2][1]-borders[2][0])/(double)nb_boites[2])) ; 
bornes[0]=((bidx[0]-DeltaBoites)<0?0:(bidx[0]-DeltaBoites)) ;
bornes[2]=((bidx[1]-DeltaBoites)<0?0:(bidx[1]-DeltaBoites)) ;
bornes[4]=((bidx[2]-DeltaBoites)<0?0:(bidx[2]-DeltaBoites)) ;
bornes[1]=((bidx[0]+DeltaBoites)>=nb_boites[0]?(nb_boites[0]-1):(bidx[0]+DeltaBoites)) ;
bornes[3]=((bidx[1]+DeltaBoites)>=nb_boites[1]?(nb_boites[1]-1):(bidx[1]+DeltaBoites)) ;
bornes[5]=((bidx[2]+DeltaBoites)>=nb_boites[2]?(nb_boites[2]-1):(bidx[2]+DeltaBoites)) ;
}

int Coarse::calc_deltaboites (int type)
{
 double deltamax=0 ; double  one=1 ; 
 double delta[6]={0,0,0,0,0,0} ; 
 static bool first=true ; 
 if (actions["setdeltaboite"].set) 
 { if (first) printf("Delta boites : %d\n", (int)ceil(std::max(one, deltamax))) ; 
   first=false ; 
   return (actions["setdeltaboite"]["valeur"]) ; 
 }
 //if (first) {DISP_Info("calc_deltaboites a été desactivé, le code ne marchait pas trop ......\n") ; first=false ; }
 if (type==0 || type==2)
 {
   delta[0]=(sigmafenetre(1)+actions.Cst["Radius"])/((borders[0][1]-borders[0][0])/(double)nb_boites[0]) ; delta[0]=std::min(delta[0],(double)nb_boites[0]) ;
   delta[1]=(sigmafenetre(2)+actions.Cst["Radius"])/((borders[1][1]-borders[1][0])/(double)nb_boites[1]) ; delta[1]=std::min(delta[1],(double)nb_boites[1]) ;
   if (!actions["is2D"].set)
      delta[2]=(sigmafenetre(3)+actions.Cst["Radius"])/((borders[2][1]-borders[2][0])/(double)nb_boites[2]) ; delta[2]=std::min(delta[2],(double)nb_boites[2]) ;
 }
 if (type==1 || type==2)
 {
   delta[3]=(sigmafenetreO(1)+actions.Cst["Radius"])/((borders[0][1]-borders[0][0])/(double)nb_boites[0]) ; delta[3]=std::min(delta[3],(double)nb_boites[0]) ;
   delta[4]=(sigmafenetreO(2)+actions.Cst["Radius"])/((borders[1][1]-borders[1][0])/(double)nb_boites[1]) ; delta[4]=std::min(delta[4],(double)nb_boites[1]) ;
   if (!actions["is2D"].set)
      delta[5]=(sigmafenetreO(3)+actions.Cst["Radius"])/((borders[2][1]-borders[2][0])/(double)nb_boites[2]) ; delta[5]=std::min(delta[5],(double)nb_boites[2]) ;
 }
 deltamax=std::max(delta[0],delta[1]) ; deltamax=std::max(delta[2],deltamax) ; deltamax=std::max(delta[3],deltamax) ; 
 deltamax=std::max(delta[4],deltamax) ; deltamax=std::max(delta[5],deltamax) ; 
 if (first) printf("Delta boites : %d\n", (int)ceil(std::max(one, deltamax))) ; 
 first=false ; 
 return ((int)ceil(std::max(one, deltamax))-1) ;  
}
*/
//==========================================================================
//==========================================================================
/*Boites operator+ (Boites b1, Boites b2)
{
  Boites temp ; 
  
 if (b1==b2)
   temp=b1 ; 
 else
 {cout << "ERR : tentative d'ajout de boîtes non compatibles\n" ; return temp ; }

 temp.nb_atomes=b1.nb_atomes+b2.nb_atomes ; 
 temp.nb_contacts=b1.nb_contacts+b2.nb_contacts ; 
 temp.vx=b1.vx+b2.vx ; 
 temp.vy=b1.vy+b2.vy ; 
 temp.vz=b1.vz+b2.vz ; 
 temp.sigma=b1.sigma+b2.sigma ; 
 temp.sigmak=b1.sigmak+b2.sigmak ;
 temp.sigmatot=b1.sigmatot+b2.sigmatot ;
 temp.phi=b1.phi+b2.phi ;
 temp.mean_radius=b1.mean_radius+b2.mean_radius ;
 return temp ; 
}

//--------------------------------------------------------------------------
void Boites::operator= (Boites b)
{
  x0=b.x0 ; 
  nb_atomes=b.nb_atomes ; 
  vx=b.vx ; vy=b.vy ; vz=b.vz ; 
  phi=b.phi ; 
  mean_radius=b.mean_radius ;
  //volume=b.volume ; 
  sigma=b.sigma ; nb_contacts=b.nb_contacts ; 
  sigmak=b.sigmak ;
  sigmatot=b.sigmatot ;
}

//----------------------------------------------------------------------
bool operator== (Boites b1, Boites b2)
{
 if (b1.x0(1)==b2.x0(1) && b1.x0(2)==b2.x0(2) && b1.x0(3)==b2.x0(3))
   return true ; 
 else
   return false ; 
}
*/

/*=========================================================*/
int CoarseDump::do_coarse_basic (Dump &dump, string chemin)   // Coarse basic : seulement du LDump, avec boîte créneau de base, toujours moyennée en temps
{
// TODO
  DISP_Err("Il faudrait réimplémenter le do_coarse_basic") ;
  return 1 ; 
/*long int debut, fin, step ;
int  i,j, idx[6] ; Coarse coarse_null ;
int bidx[3], idxclosest ;
long int loop[3] ;

actions.copyarg("coarse-graining-basic", "coarse-graining") ;
dump.loopdat(loop) ;

cout << "\nCoarseDump::do_coarse_basic          " ;
actions.total=loop[2]-loop[0] ; actions.valeur=i ;actions.disp_progress() ;
coarse.resize(1, coarse_null) ;

i=loop[0] ;
dump.check_timestep(loop[0]) ;
coarse[0].coarse(dump.steps[i]) ;  // First corse, on initialise les boîtes.
idx[0]=dump.steps[i].find_idx(IDS("POSX")) ; idx[1]=dump.steps[i].find_idx(IDS("POSY")) ; idx[2]=dump.steps[i].find_idx(IDS("POSZ")) ;
idx[3]=dump.steps[i].find_idx(IDS("VX")) ; idx[4]=dump.steps[i].find_idx(IDS("VY")) ; idx[5]=dump.steps[i].find_idx(IDS("VZ")) ;

for (i=loop[0] ; i<loop[2] ; i+=loop[1])
    {
    actions.valeur=i ;
    dump.check_timestep(i) ;
    for (j=0 ; j<dump.steps[i].nb_atomes ; j++)
     {
     bidx[0]=floor((dump.steps[i].datas[idx[0]][j]-coarse[0].borders[0][0])/(double)((coarse[0].borders[0][1]-coarse[0].borders[0][0])/(double)coarse[0].nb_boites[0])) ;
     bidx[1]=floor((dump.steps[i].datas[idx[1]][j]-coarse[0].borders[1][0])/(double)((coarse[0].borders[1][1]-coarse[0].borders[1][0])/(double)coarse[0].nb_boites[1])) ;
     bidx[2]=floor((dump.steps[i].datas[idx[2]][j]-coarse[0].borders[2][0])/(double)((coarse[0].borders[2][1]-coarse[0].borders[2][0])/(double)coarse[0].nb_boites[2])) ;
     bidx[0]=bidx[0]<0?0:bidx[0] ; bidx[1]=bidx[1]<0?0:bidx[1] ; bidx[2]=bidx[2]<0?0:bidx[2] ;
     bidx[0]=bidx[0]>=coarse[0].nb_boites[0]?coarse[0].nb_boites[0]-1:bidx[0] ;
     bidx[1]=bidx[1]>=coarse[0].nb_boites[1]?coarse[0].nb_boites[1]-1:bidx[1] ;
     bidx[2]=bidx[2]>=coarse[0].nb_boites[2]?coarse[0].nb_boites[2]-1:bidx[2] ;

     idxclosest=coarse[0].linidx(bidx[0], bidx[1],bidx[2]) ;

     coarse[0].boites[idxclosest].vx+=dump.steps[i].datas[idx[3]][j] ;
     coarse[0].boites[idxclosest].vy+=dump.steps[i].datas[idx[4]][j] ;
     coarse[0].boites[idxclosest].vz+=dump.steps[i].datas[idx[5]][j] ;
     coarse[0].boites[idxclosest].nb_atomes+=1 ;
     }
    }
for (i=0 ; i<coarse[0].nb_boites[0]*coarse[0].nb_boites[1]*coarse[0].nb_boites[2] ; i++)
{
	coarse[0].boites[i].vx/=coarse[0].boites[i].nb_atomes ;
	coarse[0].boites[i].vy/=coarse[0].boites[i].nb_atomes ;
	coarse[0].boites[i].vz/=coarse[0].boites[i].nb_atomes ;
	coarse[0].boites[i].phi=coarse[0].boites[i].nb_atomes ;
}

ofstream out ;
stringstream chemin2 ;

// Write as VTK
chemin2.str("") ;
chemin2 << chemin << "-BasicCoarse.vtk" ;
chemin2 << "Basic-coarse.vtk" ;
out.open(chemin2.str().c_str(), ios::out) ;
coarse[0].write_asVTK(out) ;
out.close() ;

// Write as matlab
chemin2.str("") ;
chemin2 << chemin << "-BasicCoarse.txt" ;
out.open(chemin2.str().c_str(), ios::out) ;
for (i=0 ; i<coarse[0].nb_boites[0]*coarse[0].nb_boites[1]*coarse[0].nb_boites[2] ; i++)
{
	out << coarse[0].boites[i].x0(1) << " " << coarse[0].boites[i].x0(2) << " " << coarse[0].boites[i].x0(3) << " " ;
	out << coarse[0].boites[i].vx << " " << coarse[0].boites[i].vy << " " << coarse[0].boites[i].vz << " " << coarse[0].boites[i].phi << "\n"  ;
}
out.close() ;


return 1 ;*/
}




//===============================================
// Réécriuture des fonctions de Dump pour faire l'encapsulation pour le Writing

int CoarseDump::open (string chemin) {
int i ; Step null_step ;
if (chemin.length() >0) DISP_Warn("Chemin devrait être vide ...") ;
// Ne rien faire, normalement tout a été fait lors du calcul du coarse 
return 1; 
}

int CoarseDump::check_timestep (long int timestep) 
{
int i ;
if (inscrit==timestep) return 1 ;
if (inscrit>=0) // il faut libérer la mémoire
  {
  inscrit=free_timestep(inscrit) ;
  } 

//steps[timestep].nb_atomes=coarse[timestep].nb_boites_tot ;  
//steps[timestep].nb_idx=steps[timestep].idx_col.size() ; 

// Réorganisation des cellules
Filter::filtrer(steps[timestep], "id::sort::null", TCOARSE) ; 
inscrit=timestep ;
return 2 ; 
}


int CoarseDump::free_timestep (long int inscrit) 
{
int i ; 
steps[inscrit].filtered=false ; 
inscrit=-1 ; 
return inscrit ; 
}












