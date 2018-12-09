#include "Headers/Filter.h"

Filter::Filter()
{
// Liste des opérations à ajouter pour initialiser les filtres
add_operation ("null",&Filter::do_nothing) ;
add_operation ("<", &Filter::op_compare) ;
add_operation (">", &Filter::op_compare) ;
add_operation ("<=", &Filter::op_compare) ;
add_operation (">=", &Filter::op_compare) ;
add_operation ("=", &Filter::op_compare) ;
add_operation ("!=", &Filter::op_compare) ;
add_operation ("::+::", &Filter::op_operation) ;
add_operation ("::-::", &Filter::op_operation) ;
add_operation ("::/::", &Filter::op_operation) ;
add_operation ("::*::", &Filter::op_operation) ;
add_operation ("::sum::", &Filter::op_collapse) ;
add_operation ("::mean::", &Filter::op_collapse) ;
add_operation ("::max::", &Filter::op_collapse) ;
add_operation ("::min::", &Filter::op_collapse) ;
add_operation ("::sumby::", &Filter::op_collapseby) ;
add_operation ("::meanby::", &Filter::op_collapseby) ;
add_operation ("::sumbyt0::", &Filter::op_collapseby) ;
add_operation ("::meanbyt0::", &Filter::op_collapseby) ;
add_operation ("::2Dto3D::", &Filter::make2D3D) ;
add_operation ("::pbc::", &Filter::swappbc) ;
add_operation ("::lieedwardpbc::", &Filter::lieedwardswap) ;
add_operation ("::moveref::", &Filter::moveref) ; // Change referentiel function of time ex: vx::moveref::0.1
add_operation ("::changeref::", &Filter::translate_ref) ;
add_operation ("::changerefx::", &Filter::translate_refx) ;
add_operation ("::changerefy::", &Filter::translate_refy) ;
add_operation ("::changerefz::", &Filter::translate_refz) ;
add_operation ("::sort::", &Filter::quicksort_base) ;
add_operation ("::fill::", &Filter::fill) ;
add_operation ("::zrotate::", &Filter::zrotate) ;
add_operation ("::quartx+y+::", &Filter::moveto_xplusyplus) ;
add_operation ("::%sup::", &Filter::pourcentage) ;
add_operation ("::%inf::", &Filter::pourcentage) ;
add_operation ("::.cfpp.::", &Filter::cf_postproc) ;
add_operation ("::.cpplink.::", &Filter::copypartiallink) ;
add_operation ("::.wallforceatm.::", &Filter::wallforceatm) ;
add_operation ("::multisphere::", &Filter::multisphere) ;
add_operation ("::multisphere_remove_atm::", &Filter::multisphere_remove_atm) ;
add_operation ("::cluster::", &Filter::cluster) ;
add_operation ("::quat2orient::", &Filter::quat2orient);
add_operation ("::.nopbc.::", &Filter::no_periodic_chain) ;
}

//============================================
int Filter::add_operation (const char nom[], pf fonction)
{string tmp ; tmp=nom ; Operations.insert ( pair <string,pf> (tmp,fonction) ); return 1 ; }

int Filter::filtrer (Step &step, const char * commande, int TypeDump)
{
 Filter tmpfiltre ;
 vector <struct Op> operations ;

 operations=tmpfiltre.parse_arg(commande, TypeDump) ;
 tmpfiltre.do_filter(step, operations) ;
 return 1 ;
}

//=============================================
vector <struct Op> Filter::parse_arg (const char arg[], int TypeDump)
{
int nb_item, i, j ;
char tmp_commande[500] ;
vector <struct Op> retour ;
j=0 ; nb_item=0 ;
for (i=0 ; i<2000 ; i++)
  {
  if (arg[i]==';' || arg[i]=='"' || arg[i]==0)
    {
    nb_item++ ;
    tmp_commande[j]=0 ;
    retour.push_back(parse_operation(tmp_commande, j, TypeDump)) ;
    j=0 ;
    }
  else
    {
    tmp_commande[j]=arg[i] ;
    j++ ;
    }
  if (arg[i]==0) break ;
  }
return retour ;
}
//--------------------------------
struct Op Filter::parse_operation (char com[500], int length, int TypeDump)
{
int i, j=0, change=0 ; unsigned char type ;
struct Op op_tmp ;
char id[100], op[100], valeur[100] ;

for (i=0 ; i<length ; i++)
  {
  if (com[i]=='<' || com[i]=='>' || com[i]=='=' || com[i]=='!')
     {
     if (change!=4) {change=4 ; id[j]=0 ; j=0 ;}
     }
  else
     {
     if (change==4)
        {change++ ; op[j]=0 ; j=0 ; }
     if (com[i]==':')
       change++ ;
     if (change==1) {id[j]=0 ; j=0 ;}
     }

  switch (change)
    {
    case 0     : id[j]=com[i] ; break ;
    case 1 ... 4 : op[j]=com[i] ; break ;
    case 5     : valeur[j]=com[i] ; break ;
    }
  j++ ;
  }
valeur[j]=0 ;

int res ;
res=IDS(id) ;
if (res==-1)
{DISP_Warn ("Un type de donné n'est pas référencé, IDS va l'ajouter (avec un type TL) \n") ;
 res=IDS.new_id(id, TL) ;
 if (res!=-1) op_tmp.idx1=res ;
 else
 {
   DISP_Err("Le type inconnu n'a pas pu être ajouté !!\n");
   op_tmp.idx1=IDS("UNKNOWN") ;
 }
}
else
 op_tmp.idx1=res ;


if (((valeur[0]>65 && valeur[0]<90) || (valeur[0]>97 && valeur[0]<122)) && (strncmp(valeur, "nan",3)!=0) && (strncmp(valeur,"NAN",3)!=0) && (strncmp(valeur,"NaN",3)!=0))
  {
  type=1 ;
  res=IDS(valeur) ;
  if (res==-1)
    { DISP_Warn ("Un type de donné n'est pas référencé, IDS va l'ajouter\n") ;
      res=IDS.new_id(id, TL) ;
      if (res!=-1) op_tmp.idx2=res ;
      else
      {
	DISP_Err("Le type inconnu n'a pas pu être ajouté !!\n");
	op_tmp.idx2=IDS("UNKNOWN") ;
      }
    }
  else
    op_tmp.idx2=res ;
  }
else
  {
  type=0 ;
  sscanf(valeur,"%lf",&op_tmp.valeur) ; if(isnan(op_tmp.valeur)) printf("Nan in filter value, proceeding\n") ;
  }

if (type==0)
  {
  if (Operations.count(op)>0) op_tmp.operation=op ;
  else {cerr << "WARN2 : le type d'operation" << op <<"est inconnu pour le filtre\n" ; }
  }
else if (type==1)
  {
  if (Operations.count(op)>0) op_tmp.operation=op ;
  else {cerr << "WARN2 : le type d'operation" << op <<"est inconnu pour le filtre\n" ; }
  }
op_tmp.type=type ;
op_tmp.TypeDump=TypeDump ;
return op_tmp ;
}
//----------------------------------------------------
void Filter::disp (vector <struct Op> lst_op)
{
unsigned int i ;
for (i=0 ; i<lst_op.size() ; i++)
  {
  cout << "Operation ["<< lst_op[i].operation << "] sur [" << lst_op[i].idx1 << "] par [" << lst_op[i].idx2 << "] ou ["<< lst_op[i].valeur<<"] type de dump [" ;
  if      (lst_op[i].TypeDump==1) cout << "TF].\n" ;
  else if (lst_op[i].TypeDump==2) cout << "TL].\n" ;
  else if (lst_op[i].TypeDump==3) cout << "TCF].\n" ;
  else if (lst_op[i].TypeDump==16) cout << "TCOARSE].\n" ;
  else    cout << "??? (bizarre...)].\n" ;
  }
}

//----------------------------------------------------
int Filter::do_filter (Step &step, vector <struct Op> lst_op)
{
  unsigned int i ;
  int retour ;

  for (i=0 ; i<lst_op.size() ; i++)
  {
    if (step.Type!=lst_op[i].TypeDump) {/*DISP_Info ("Un filtre n'a pas été appliqué, les types de dump ne correspondaient pas\n") ;*/ continue ; }
  // Execution de la commande nécessaire
  retour = (this->*(Operations[lst_op[i].operation]))(step, lst_op[i]) ;
  }
return 1 ;
}
//=======================================================
int Filter::swappbc (Step &step, struct Op op)
{
int i, idx ;

idx=step.find_idx(op.idx1) ;
for (i=0 ; i<step.nb_atomes; i++)
    {
	if (step.datas[idx][i]>op.valeur)
		step.datas[idx][i]-=2*op.valeur ;
	else if (step.datas[idx][i]<-op.valeur)
		step.datas[idx][i]+=2*op.valeur ;
    }
return 1 ;
}
//-------------------------------------------------------
int Filter::lieedwardswap (Step &step, struct Op op)
{
 double box[6] ; int idx[3],i ; double a,b,c ;
 // Conversion de la bounding box from liggghts documentation
 box[0]=step.box[0][0]-fmin(0, fmin(step.triclinic[0], fmin(step.triclinic[1], step.triclinic[0]+step.triclinic[1]))) ;
 box[1]=step.box[0][1]-fmax(0, fmax(step.triclinic[0], fmax(step.triclinic[1], step.triclinic[0]+step.triclinic[1]))) ;
 box[2]=step.box[1][0]-fmin(0, step.triclinic[2]) ;
 box[3]=step.box[1][1]-fmax(0, step.triclinic[2]) ;
 box[4]=step.box[2][0] ; box[5]=step.box[2][1] ;
 //printf("%g %g %g %g %g %g %g\n", box[0], box[1], box[2], box[3], box[4], box[5]) ; fflush(stdout) ;
 if (step.Type!=TL) {DISP_Err("LieEdward PBC non implémenté en type non atomique") ; return -1 ; }

 idx[0]=step.find_idx(IDS("POSX")) ;
 idx[1]=step.find_idx(IDS("POSY")) ;
 idx[2]=step.find_idx(IDS("POSZ")) ;

 step.box[0][0]=box[0] ; step.box[0][1]=box[1] ; step.box[1][0]=box[2] ; step.box[1][1]=box[3] ; step.box[2][0]=box[4] ; step.box[2][1]=box[5] ;
 step.triclinic[0]=step.triclinic[1]=step.triclinic[2]=0 ;

 for (i=0 ; i<step.nb_atomes ; i++)
 {
  // printf("%d{%g %g %g %g}", i, step.datas[idx[0]][i],step.datas[idx[1]][i], step.datas[idx[2]][i] ) ; fflush(stdout) ;
  a= (step.datas[idx[0]][i]-box[0])/(box[1]-box[0]) ;
  b= (step.datas[idx[1]][i]-box[2])/(box[3]-box[2]) ;
  c= (step.datas[idx[2]][i]-box[4])/(box[5]-box[4]) ;

  if (a>=1 || a<0) step.datas[idx[0]][i]-=floor(a)*(box[1]-box[0]) ;

  if (b>=1 || b<0) step.datas[idx[1]][i]-=floor(b)*(box[3]-box[2]) ;

  if (c>=1 || c<0) step.datas[idx[2]][i]-=floor(c)*(box[5]-box[4]) ;

  if (step.datas[idx[0]][i]<box[0] || step.datas[idx[0]][i]>box[1] || step.datas[idx[1]][i]<box[2] || step.datas[idx[1]][i]>box[3] || step.datas[idx[2]][i]<box[4] || step.datas[idx[2]][i]>box[5])
    DISP_Err("Erreur dans le filtre LieEdward") ;
 }
 return 1 ;
}
//--------------------------------------------------------
int Filter::moveref (Step &step, struct Op op)
{
 int idx[9], id ; int i, ts ;
 static int initial_ts=-1, warn=-1 ;

 if (initial_ts==-1) initial_ts=step.timestep ;

 id=step.find_idx(op.idx1) ;
 ts=step.timestep-initial_ts ;

 // constante translation in positions
 idx[0]=step.find_idx(IDS("POSX")) ;
 idx[1]=step.find_idx(IDS("POSY")) ;
 idx[2]=step.find_idx(IDS("POSZ")) ;
 if (id==idx[0] || id==idx[1] || id==idx[2])
 {
  for (i=0 ; i<step.nb_atomes ; i++) step.datas[id][i]+=op.valeur ;
  return 0 ;
 }

 // referentiel moving at constante velocity
 idx[3]=step.find_idx(IDS("VX")) ;
 idx[4]=step.find_idx(IDS("VY")) ;
 idx[5]=step.find_idx(IDS("VZ")) ;
 if (id==idx[3] || id==idx[4] || id==idx[5])
 {
  if ( id==idx[3] )
  {for (i=0 ; i<step.nb_atomes ; i++)
      {
      step.datas[id][i]-=op.valeur ;
      step.datas[idx[0]][i]-=op.valeur*ts ;
      }}
  else if ( id==idx[4] )
  {for (i=0 ; i<step.nb_atomes ; i++)
      {
      step.datas[id][i]-=op.valeur ;
      step.datas[idx[1]][i]-=op.valeur*ts ;
      }}
  else
  {for (i=0 ; i<step.nb_atomes ; i++)
      {
      step.datas[op.idx1][i]-=op.valeur ;
      step.datas[idx[2]][i]-=op.valeur*ts ;
      }}

 return 0 ;
 }

 //referential moving at constant ACCELERATION (yes ...)
 if (op.idx1==IDS("FX") || op.idx1==IDS("FY") || op.idx1==IDS("FZ")) // referentiel moving at constante velocity
 {
  if (warn==-1) {DISP_Warn ("C'est un filtre SPECIAL PURPOSE. Il faut le vérifier avant de l'utiliser.\n") ; warn=1 ; }
  if ( op.idx1==IDS("FX") )
  {for (i=0 ; i<step.nb_atomes ; i++)
      {
      step.datas[idx[3]][i]-=op.valeur*ts/1000 ;
      step.datas[idx[0]][i]-=op.valeur*ts/1000*ts/1000/2 ;
      }}
  else if ( op.idx1==IDS("FY") )
  {for (i=0 ; i<step.nb_atomes ; i++)
      {
      step.datas[idx[4]][i]-=op.valeur*ts/1000 ;
      step.datas[idx[1]][i]-=op.valeur*ts/1000*ts/1000/2 ;
      }}
  else
  {for (i=0 ; i<step.nb_atomes ; i++)
      {
      step.datas[idx[5]][i]-=op.valeur*ts/1000 ;
      step.datas[idx[2]][i]-=op.valeur*ts/1000*ts/1000/2 ;
      }}

 return 0 ;
 }

 DISP_Warn ("Impossible d'appliquer le filtre de referentiel mobile\n") ;
 return -1 ;

}


//-------------------------------------------------------
int Filter::translate_refgen (Step &step, struct Op op, int dir)
{
double pos[3] ;
int idx[4] ;
int i ; bool ok=false ;

if (step.Type==TL)
{
 idx[0]=step.find_idx(IDS("POSX")) ;
 idx[1]=step.find_idx(IDS("POSY")) ;
 idx[2]=step.find_idx(IDS("POSZ")) ;
}
else if (step.Type==TCF)
{
 idx[0]=step.find_idx(IDS("CFPOSX")) ;
 idx[1]=step.find_idx(IDS("CFPOSY")) ;
 idx[2]=step.find_idx(IDS("CFPOSZ")) ;
}
else DISP_Warn("Filter::translate_refgen : type de dump inconnu pour ce filtre\n") ;

idx[3]=step.find_idx(op.idx1) ;
for (i=0 ; i<step.nb_atomes ; i++)
 {
 if (step.datas[idx[3]][i]==op.valeur)
    {
	ok=true ;
	pos[0]=step.datas[idx[0]][i] ;
	pos[1]=step.datas[idx[1]][i] ;
	pos[2]=step.datas[idx[2]][i] ;
	break ;
    }
 }

if (ok==false)
 {
 DISP_Warn("Erreur : impossible d'effectuer la translation de referentiel")  ;
 return 0 ;
 }
//printf("%f %f %f|", pos[0], pos[1], pos[2]) ;

for (i=0 ; i<step.nb_atomes ; i++)
 {
 if (dir & 1) step.datas[idx[0]][i]-=pos[0] ;
 if (dir & 2) step.datas[idx[1]][i]-=pos[1] ;
 if (dir & 4) step.datas[idx[2]][i]-=pos[2] ;
 }
return 1 ;
}
//-------------------------------------------------------
int Filter::op_operation (Step &step, struct Op op)
{
int idx, i ;
idx=step.find_idx(op.idx1) ;
for (i=0 ; i<step.nb_atomes ; i++)
    {
    if (op.operation=="::+::") {step.datas[idx][i]+=op.valeur ; }
    if (op.operation=="::-::") {step.datas[idx][i]-=op.valeur ; }
    if (op.operation=="::/::") {step.datas[idx][i]/=op.valeur ; }
    if (op.operation=="::*::") {step.datas[idx][i]*=op.valeur ; }
    }
return 1;
}
//--------------------------------------------------------
int Filter::op_collapseby (Step &step, struct Op op)
{
int i, j, k ;
static bool t0=true ;
static vector <double> keys ;
static vector < vector <int> > ids1 ;
static vector < vector <int> > ids2 ;
int idx[3] ; int group ; int nogroup=0 ;

if (step.Type==TL) idx[0]=step.find_idx(IDS("ID")) ;
else if (step.Type==TCF) { idx[0]=step.find_idx(IDS("CFID1")) ; idx[1]=step.find_idx(IDS("CFID2")) ; }
else {DISP_Err ("Dump inconnu pour l'opération collapseby\n") ; return 3 ; }

// Création initiale des groupes d'atomes
if (op.operation=="::sumby::" || op.operation=="::meanby::")
{
  DISP_Err ("sumby and meanby not implemented (yet?)") ; return 2 ;
}
else if (t0 && (op.operation=="::sumbyt0::" || op.operation=="::meanbyt0::"))
{
 DISP_Info ("Filter sumbyt0 ou meanbyt0 : Les valeurs individuelles sont calculées au premier pas de temps. Les ID des particules / chaînes sont supposés constants au cours du temps.\n") ;

 IDS.new_id("group",step.Type) ;
 IDS.new_id("groupnumber",step.Type) ;

 idx[2]=step.find_idx(op.idx1) ;
 t0=false ;
 // Création des clefs et des IDs d'atomes correspondant à ces clefs.
 for (i=0 ; i<step.nb_atomes ; i++)
 {
   for (j=0 ; j<keys.size() && keys[j]!=step.datas[idx[2]][i] ; j++) ;
   if (j==keys.size()) // Il faut ajouter la nouvelle clef
   {
    keys.push_back(step.datas[idx[2]][i]) ;
    ids1.resize(ids1.size()+1) ;
    ids1[ids1.size()-1].push_back(step.datas[idx[0]][i]) ;
    if (step.Type==TCF)
    {
      ids2.resize(ids2.size()+1) ;
      ids2[ids2.size()-1].push_back(step.datas[idx[1]][i]) ;
    }
   }
   else
   {
     ids1[j].push_back(step.datas[idx[0]][i]) ;
     if (step.Type==TCF) ids2[j].push_back(step.datas[idx[1]][i]) ; ;
   }
 }
 // Pour accélérer la recherche ultérieure, on va ranger les IDs. Heu ... non avec 2 clefs c'est relou
 /*for (i=0 ; i<keys.size() ; i++)
 {
  int change ; double tmp ;
  do {
  change=0 ;
  for (j=0 ; j<ids[i].size()-1 ; j++)
    if (ids[i][j]>ids[i][j+1])
    { tmp=ids[i][j] ; ids[i][j]=ids[i][j+1] ; ids[i][j+1]=tmp ; change++ ; }
  } while (change!=0) ;
 }*/
 //for (i=0 ; i<keys.size() ; i++) {printf("%f :", keys[i]) ; for (j=0 ; j<ids1[i].size() ; j++) printf("%d/%d ", ids1[i][j],ids2[i][j]) ; printf("\n") ; fflush(stdout) ; }
 printf("%lu groupes\n", keys.size() ) ;

}

// Ajout des N groupes au tableau data par envoie des N premières lignes en fin de tableau
step.idx_col.push_back(IDS("GROUP")) ;
step.idx_col.push_back(IDS("GROUPNUMBER")) ;
step.datas.resize(step.datas.size()+2) ;
step.datas[step.datas.size()-2].resize(keys.size(),0) ;
step.datas[step.datas.size()-1].resize(keys.size(),0) ;
// step.nb_idx+=2 <---- Fait seulement à la fin, en même temps que le changement de nombre d'atomes.
for (i=0 ; i<keys.size() ; i++)
{
  step.copy_atm_end(i) ;
  step.datas[step.datas.size()-2][i]=keys[i] ;
  for (j=0 ; j<step.nb_idx ; j++)
    step.datas[j][i]=0 ;
}


for (i=keys.size() ; i<step.nb_atomes+keys.size() ; i++)
{
  group=-1 ;
  for (j=0 ; j<keys.size() ; j++)
   for (k=0 ; k<ids1[j].size() && ids1[j][k]<=step.datas[idx[0]][i] ; k++)
     if ((ids1[j][k]==step.datas[idx[0]][i] && ids2[j][k]==step.datas[idx[1]][i])
       ||(ids1[j][k]==step.datas[idx[1]][i] && ids2[j][k]==step.datas[idx[0]][i]))
     {group=j ; j=keys.size() ; break ; }
  if (group==-1) {nogroup++ ; continue ; }

  for (j=0 ; j<step.nb_idx ; j++)
   step.datas[j][group]+=step.datas[j][i] ;
  step.datas[step.datas.size()-1][group]++ ;

}
if (nogroup>0)
{
 printf("%d ", nogroup) ; DISP_Info ("atomes n'ont pas eu de groupe attribués\n") ;
}

if (op.operation=="::meanbyt0::")
{
 for (i=0 ; i<keys.size() ; i++)
   for (j=0 ; j<step.nb_idx ; j++)
   {
     if (step.datas[step.datas.size()-1][i]>0)
       step.datas[j][i]/=step.datas[step.datas.size()-1][i] ;
     else
       step.datas[j][i]=0 ;
   }
}
step.nb_atomes=keys.size() ;
step.nb_idx+=2 ;
return 1;
}
//-------------------------------------------------------
int Filter::op_collapse(Step &step, struct Op op)
{
int i, j ; int optype = 0 ;
if (op.operation=="::sum::") optype=1 ;
else if (op.operation=="::mean::") optype=2 ;
else if (op.operation=="::max::") optype=3 ;
else if (op.operation=="::min::") optype=4 ;
else DISP_Warn ("Undefined op.operation (not implemented) \n") ;

for (j=0 ; j<step.nb_idx ; j++)
{
for (i=1 ; i<step.nb_atomes ; i++)
    {
      switch (optype) {
        case 1 :
        case 2 : step.datas[j][0]+=step.datas[j][i] ; break ;
        case 3   : if (step.datas[j][0] < step.datas[j][i]) step.datas[j][0]=step.datas[j][i] ; break ;
        case 4   : if (step.datas[j][0] > step.datas[j][i]) step.datas[j][0]=step.datas[j][i] ; break ;
        default : step.datas[j][0]=0 ; break ;
      }
    }
}
if (optype==1) // op.operation=="::mean::"
{
 for (j=0 ; j<step.nb_idx ; j++)
    step.datas[j][0]/=step.nb_atomes ;
}
step.nb_atomes=1 ;
return 1;
}

//-------------------------------------------------------
int Filter::zrotate (Step &step, struct Op op)
{
int i;
int idx[3] ;
Matrix3x3 rot ; Vector pt1, axe ;

//if (step.Type != TF) { cout << "ERR : le filtre ::zrotate:: ne peut être appliqué qu'au FDump\n" ; return -1 ; }
axe(1)=axe(2)= 0.0 ; axe(3)=1.0 ;
rot=Geometrie::get_rot_matrix(op.valeur, axe) ;
// Pour les type F
if (step.Type==TF)
 {
 for (i=0 ; i<step.nb_pts ; i++)
   step.atm_rotate(rot, i+step.nb_pts) ;
 for(i=0 ; i<step.nb_triangles ; i++)
   step.atm_rotate(rot, i) ;
 }
else
 {
 for(i=0 ; i<step.nb_atomes ; i++)
    step.atm_rotate(rot, i) ;
 }

return 1 ;
}

//----------------------------------------------------
int Filter::op_compare (Step &step, struct Op op)
{
int i,j=0 ;
int idx, idx2 ;
if (op.idx1>=128 && op.idx1<192) step.nb_atomes=step.nb_triangles ;
if (op.type==0)
{
  idx=step.find_idx(op.idx1) ;
  if (isnan(op.valeur))
  {
   for (i=0 ; i<step.nb_atomes-j ; i++)
    {
    if (op.operation=="<") {DISP_Warn("Filter incompatible avec nan") ;}
    if (op.operation==">") {DISP_Warn("Filter incompatible avec nan") ;}
    if (op.operation=="<="){DISP_Warn("Filter incompatible avec nan") ;}
    if (op.operation==">="){DISP_Warn("Filter incompatible avec nan") ;}
    if (op.operation=="=") {if (!isnan(step.datas[idx][i])) {step.swap_atm(i, step.nb_atomes-1-j) ; j++ ; i-- ; }}
    if (op.operation =="!="){if (isnan(step.datas[idx][i])) {step.swap_atm(i, step.nb_atomes-1-j) ; j++ ; i-- ; }}
    }
  }
  else
  {
   for (i=0 ; i<step.nb_atomes-j ; i++)
    {
    if (op.operation=="<") {if (step.datas[idx][i]>=op.valeur) {step.swap_atm(i, step.nb_atomes-1-j) ; j++ ; i-- ; }}
    if (op.operation==">") {if (step.datas[idx][i]<=op.valeur) {step.swap_atm(i, step.nb_atomes-1-j) ; j++ ; i-- ; }}
    if (op.operation=="<="){if (step.datas[idx][i]>op.valeur) {step.swap_atm(i, step.nb_atomes-1-j) ; j++ ; i-- ; }}
    if (op.operation==">="){if (step.datas[idx][i]<op.valeur) {step.swap_atm(i, step.nb_atomes-1-j) ; j++ ; i-- ; }}
    if (op.operation=="=") {if (step.datas[idx][i]!=op.valeur) {step.swap_atm(i, step.nb_atomes-1-j) ; j++ ; i-- ; }}
    if (op.operation =="!="){if (step.datas[idx][i]==op.valeur) {step.swap_atm(i, step.nb_atomes-1-j) ; j++ ; i-- ; }}
    }
  }
}

else if (op.type==1)
  {
  idx=step.find_idx(op.idx1) ; idx2=step.find_idx(op.idx2) ;
  for (i=0 ; i<step.nb_atomes-j ; i++)
    {
    if (op.operation=="<") {if (step.datas[idx][i]>=step.datas[idx2][i]) {step.swap_atm(i, step.nb_atomes-1-j) ; j++ ; i-- ; }}
    if (op.operation==">") {if (step.datas[idx][i]<=step.datas[idx2][i]) {step.swap_atm(i, step.nb_atomes-1-j) ; j++ ; i-- ; }}
    if (op.operation=="<=") {if (step.datas[idx][i]> step.datas[idx2][i]) {step.swap_atm(i, step.nb_atomes-1-j) ; j++ ; i-- ; }}
    if (op.operation==">=") {if (step.datas[idx][i]< step.datas[idx2][i]) {step.swap_atm(i, step.nb_atomes-1-j) ; j++ ; i-- ; }}
    if (op.operation=="=") {if (step.datas[idx][i]!=step.datas[idx2][i]) {step.swap_atm(i, step.nb_atomes-1-j) ; j++ ; i-- ; }}
    if (op.operation=="!=") {if (step.datas[idx][i]==step.datas[idx2][i]) {step.swap_atm(i, step.nb_atomes-1-j) ; j++ ; i-- ; }}
    }
  }
else
  {
  cout << "ERR : le nom de l'opération à effectuer n'est pas correct\n" ;
  }

if (j>0)
   {
   step.del_end_atms(j) ;
   step.nb_atomes-=j ;
   if (op.idx1>=128 && op.idx1<192) step.nb_triangles=step.nb_atomes ;
   }
return 1 ;
}
//=====================================================

void Filter::echanger (vector< vector<double> > &tableau, int idx1, int idx2, int nbidx)
{
int i ; double tmp ;

for (i=0 ; i<nbidx ; i++)
 {
 tmp=tableau[i][idx1] ;
 tableau[i][idx1]=tableau[i][idx2] ;
 tableau[i][idx2]=tmp ;
 }
}

//-----------------------------------
int Filter::partitionner (vector< vector<double> > &tableau, int premier, int dernier, int pivot, int idx, int nbidx)
{
 int i, j;

 j=premier ;
 echanger (tableau, pivot, dernier, nbidx) ;
 for (i=premier ; i<dernier ; i++)
 {
  if (tableau[idx][i]<=tableau[idx][dernier])
  {
    echanger(tableau, i, j, nbidx) ;
    j++ ;
  }
 }
 echanger (tableau, j, dernier, nbidx) ;
 return j ;
}

//------------------------------------
int Filter::quicksort_base (Step & step, struct Op op)
{
int tmp_idx1 ; tmp_idx1=step.find_idx(op.idx1) ;
int maxrecurse ;
maxrecurse=ceil(log2(step.nb_atomes))*10 ;
quicksort (step.datas, 0, step.nb_atomes-1, tmp_idx1, step.nb_idx, maxrecurse) ;

// Check temporaire
int i ;
for (i=1 ; i<step.nb_atomes ; i++)
 {
 if(step.datas[tmp_idx1][i-1]>step.datas[tmp_idx1][i]) {DISP_Err("ERR SORT\n") ;}
 }

return 0 ;
}

int Filter::quicksort (vector<vector<double> > &tableau, int premier, int dernier, int idx, int nbidx, int recurseleft)
{
 int pivot ;
 //pivot=round((premier+(rand()/(double)RAND_MAX)*(dernier-premier))) ;
 if (recurseleft<1)
   {
   //DISP_Info("Too many nested quicksort. Let's Heapsort !\n") ;
   heapsort(tableau, premier, dernier, idx, nbidx) ;
   return 1 ;
   }

 if (premier<dernier)
 {
 //pivot=premier+(dernier-premier)/2 ;
 pivot=premier ;
 pivot=partitionner(tableau, premier, dernier, pivot, idx, nbidx) ;
 //if (pivot-1-premier>15)
 quicksort(tableau, premier, pivot-1, idx, nbidx, recurseleft-1) ;
 //else selectsort(tableau, premier, pivot-1, idx, nbidx) ;
 //if (dernier-pivot-1>15)
 quicksort(tableau, pivot+1, dernier, idx, nbidx, recurseleft-1) ;
 //else selectsort(tableau, pivot+1, dernier, idx, nbidx) ;
 }
return 0 ;
}

int Filter::tamiser(vector <vector<double> > &tableau, int noeud, int max, int deb, int idx, int nbidx)
{
int k=noeud ;
int j=2*k ;
while (j<=max)
  {
  if (j<max && tableau[idx][j-1+deb]<tableau[idx][j+deb])
    j++ ;

  if(tableau[idx][k-1+deb]<tableau[idx][j-1+deb])
    {
    echanger(tableau, k-1+deb, j-1+deb, nbidx) ;
    k=j ;
    j=2*k ;
    }
  else
    return 0;
  }
return 1 ;
}

int Filter::heapsort(vector <vector<double> > &tableau, int deb, int fin, int idx, int nbidx)
{
int i, n ;
n=fin-deb+1 ;
for (i=n/2 ; i>=1 ; i--)
  tamiser(tableau,i,n,deb,idx,nbidx) ;

for (i=n ; i>=2 ; i--)
 {
 echanger(tableau, i-1+deb, 0+deb,nbidx) ;
 tamiser(tableau, 1, i-1, deb,idx,nbidx) ;
 }
return 1 ;
}


int Filter::selectsort (vector<vector<double> > &tableau, int premier, int dernier, int idx, int nbidx)
{// Tri basique par sélection dès que le tableau est suffisament petit (pour éviter de dxfaire quicksort jusqu'au bout, d'après wikipedia c'est plus rapide)
int i, j, min ;

for (i=premier ; i<dernier ; i++)
 {
	min=i ;
	for (j=i+1 ; j<dernier ; j++)
	{ if (tableau[idx][j]<tableau[idx][min]) min=j ;}
	if (min>i) {echanger (tableau, i, min, nbidx) ; }
 }
return 0 ;
}

int Filter::fill(Step & step, struct Op op)
{
 int i,j ;
 int idx1 ; idx1=step.find_idx(op.idx1) ;
 std::vector<double>::iterator it;
 for (i=0 ; i<step.nb_atomes ; i++)
    { if (step.datas[idx1][i]!=i+1)
    	{
	//printf("F") ;
    	for (j=0; j<step.nb_idx ; j++)
    	  {
    	  it = step.datas[j].begin()+i;
    	  step.datas[j].insert(it,NAN) ;
    	  }
    	it = step.datas[idx1].begin()+i;
    	step.datas[idx1][i]=i+1 ;
    	step.nb_atomes++ ;
    	}
    }


}

//=======================================
int Filter::do_nothing(Step & step, struct Op op)
{
// This function do nothing (just to initialize the operation NULL for forward compatibility).
return 0 ;
}
//=========================================
int Filter::pourcentage (Step & step, struct Op op)
{
// Dans tous les cas on commence par ranger la colonne demandée
struct Op range ;
range.operation="::sort::" ;
range.idx1=op.idx1 ; range.idx2=UNKNOWN ; range.valeur=0 ;
quicksort_base(step, range) ;

// Ensuite on ne garde que ce qui intéresse
unsigned int i ; int stop ;
stop=round(step.datas[0].size()*op.valeur/100) ;
if (op.operation=="::%inf::")
 {
 for (i=0 ; i<step.datas.size() ; i++)
   {
   step.datas[i].resize(stop) ;
   }
 }
else if (op.operation=="::%sup::")
 {
 for (i=0 ; i<step.datas.size() ; i++)
   {
   step.datas[i].erase(step.datas[i].begin(),step.datas[i].begin()+stop) ;
   }
 }
return 1 ;
}
//=========================================
int Filter::cf_postproc (Step &step, struct Op op)
{
unsigned int i,j, ts ; int idx[17] ; //Matrix cart(3,1), sph(3,1) ;
Vector r, f, fnv, x1, x2, xch ;
double fn, ft ;
static bool first=true ; static double samplingcf, intersampling  ;
Step * atmstep = NULL ;

// Ajoute des données à celles extraites pour les chaînes de force.
//step.idx_col.push_back(CFMAG) ;
step.idx_col.push_back(IDS("CFPOSX")) ;
step.idx_col.push_back(IDS("CFPOSY")) ;
step.idx_col.push_back(IDS("CFPOSZ")) ;
step.idx_col.push_back(IDS("CFX")) ;
step.idx_col.push_back(IDS("CFY")) ;
step.idx_col.push_back(IDS("CFZ")) ;
step.idx_col.push_back(IDS("CFMAG")) ;
//step.idx_col.push_back(CFR) ;
//step.idx_col.push_back(CFTHETA) ;
//step.idx_col.push_back(CFPHI) ;
step.nb_idx+=7 ;
step.datas.resize(step.idx_col.size()) ;
//step.datas[step.find_idx(CFR)].resize(step.nb_atomes) ;
//step.datas[step.find_idx(CFTHETA)].resize(step.nb_atomes) ;
//step.datas[step.find_idx(CFPHI)].resize(step.nb_atomes) ;
step.datas[step.find_idx(IDS("CFMAG"))].resize(step.nb_atomes) ;
step.datas[step.find_idx(IDS("CFPOSX"))].resize(step.nb_atomes) ;
step.datas[step.find_idx(IDS("CFPOSY"))].resize(step.nb_atomes) ;
step.datas[step.find_idx(IDS("CFPOSZ"))].resize(step.nb_atomes) ;
step.datas[step.find_idx(IDS("CFX"))].resize(step.nb_atomes) ;
step.datas[step.find_idx(IDS("CFY"))].resize(step.nb_atomes) ;
step.datas[step.find_idx(IDS("CFZ"))].resize(step.nb_atomes) ;
// Récupère le bon pas de dump atomique
//i=0 ; plutôt que de partir de i=0, on va essayer de se rapprocher en supposant que le coefficient de downsampling est constant
// Attention, ne fonctionne que si le,premier appel est effectué avec le step[0] du cfdump.

// Si la position des atomes n'est pas dans le dump des contacts, onutilise le dump atomique pour les trouver
if (!actions["donotusecfpos"].set && (step.find_idx(IDS("CFID1X"))!=-1 && step.find_idx(IDS("CFID1Y"))!=-1 && step.find_idx(IDS("CFID1Z"))!=-1 &&
	step.find_idx(IDS("CFID2X"))!=-1 && step.find_idx(IDS("CFID2Y"))!=-1 && step.find_idx(IDS("CFID2Z"))!=-1))
 {
 idx[0]=step.find_idx(IDS("CFID1X")) ; idx[1]=step.find_idx(IDS("CFID1Y")) ; idx[2]=step.find_idx(IDS("CFID1Z")) ;
 idx[3]=step.find_idx(IDS("CFFORCEX")) ; idx[4]=step.find_idx(IDS("CFFORCEY")) ; idx[5]=step.find_idx(IDS("CFFORCEZ")) ;
 idx[6]=step.find_idx(IDS("CFID2X")) ; idx[7]=step.find_idx(IDS("CFID2Y")) ; idx[16]=step.find_idx(IDS("CFID2Z")) ;

 idx[8]=step.find_idx(IDS("CFPOSX")) ; idx[9]=step.find_idx(IDS("CFPOSY")) ; idx[10]=step.find_idx(IDS("CFPOSZ")) ;
 idx[11]=step.find_idx(IDS("CFX")) ; idx[12]=step.find_idx(IDS("CFY")) ; idx[13]=step.find_idx(IDS("CFZ")) ;
 idx[14]=step.find_idx(IDS("CFMAG")) ;// idx[9]=step.find_idx(IDS("CFR")) ; idx[10]=step.find_idx(IDS("CFTHETA")) ; idx[11]=step.find_idx(IDS("CFPHI")) ;
 idx[15]=step.find_idx(IDS("CFPERIOD")) ;
 fflush(stdout) ;
 for (i=0 ; i<step.nb_atomes ; i++)
   {
   step.datas[idx[14]][i]=sqrt(pow(step.datas[idx[3]][i],2)+pow(step.datas[idx[4]][i],2)+pow(step.datas[idx[5]][i],2)) ;
   step.datas[idx[11]][i]=step.datas[idx[6]][i] -step.datas[idx[0]][i] ;
   step.datas[idx[12]][i]=step.datas[idx[7]][i] -step.datas[idx[1]][i] ;
   step.datas[idx[13]][i]=step.datas[idx[16]][i]-step.datas[idx[2]][i] ;
   step.datas[idx[8]][i]=(step.datas[idx[6]][i]+step.datas[idx[0]][i])/2. ;
   step.datas[idx[9]][i]=(step.datas[idx[7]][i]+step.datas[idx[1]][i])/2. ;
   step.datas[idx[10]][i]=(step.datas[idx[16]][i]+step.datas[idx[2]][i])/2. ;
   }
 }
// Si la position des atomes n'est pas dans le dump des contacts, on utilise le dump atomique pour les trouver
else
 {
 if (first==true)
  {
  op.alter_dump->check_timestep(0) ;
  intersampling=step.timestep/(double)op.alter_dump->steps[0].timestep ;
  samplingcf=step.timestep ;
  first=false ; i=0 ;
  }
 else
  {
  i=floor(step.timestep/samplingcf)*intersampling-1 ; // Le -1 donne une marge de maneuvre, c'est mieux ...
  }
//op.alter_dump->check_timestep(i) ;
 while (op.alter_dump->steps[i].timestep != step.timestep)
   {
   if (op.alter_dump->steps[i].timestep > step.timestep)
     {cout << "ERR : impossible de trouver des timestep correspondant entre le LucDump et le LcfDump" ; return -1 ; }
    //else
    //   {i++ ; op.alter_dump->check_timestep(i) ;}
   i++ ;
   }
   op.alter_dump->check_timestep(i) ;
   ts=i ;
   // Copie de la boîte
   step.box[0][0]=op.alter_dump->steps[i].box[0][0] ;
   step.box[0][1]=op.alter_dump->steps[i].box[0][1] ;
   step.box[1][0]=op.alter_dump->steps[i].box[1][0] ;
   step.box[1][1]=op.alter_dump->steps[i].box[1][1] ;
   step.box[2][0]=op.alter_dump->steps[i].box[2][0] ;
   step.box[2][1]=op.alter_dump->steps[i].box[2][1] ;

 // Création des liens de step (important !!)

 (op.alter_dump)->steps[i].otherstep=&step ;
 step.otherstep=&((op.alter_dump->steps)[i]) ;

 // Check que les ID dans le dump atomique avancent par pas de 1
 atmstep=&(op.alter_dump->steps[i]) ;
 idx[0]=atmstep->find_idx(IDS("ID")) ;
 int warnings=0 ; static bool warn=true ;
 std::vector<double>::iterator it;
 for (i=0 ; i<atmstep->datas[idx[0]].size() ; i++)
    { if (atmstep->datas[idx[0]][i]!=i+1)
    	{
    	warnings++ ;
    	if (warn==true)
    	  { warn=false ;
    	  DISP_Warn("WARN : les ID n'augmentent pas par pas de 1 dans le LucDump. Tentative de correction par ajout d'un atome nul, qui ne devrait normalement pas être appelé") ;
	  //printf("{%d}", i+1) ;
    	  }
    	for (j=0; j<atmstep->nb_idx ; j++)
    	  {
    	  it = atmstep->datas[j].begin()+i;
    	  atmstep->datas[j].insert(it,NAN) ;
    	  }
    	it = atmstep->datas[idx[0]].begin()+i;
    	atmstep->datas[idx[0]][i]=i+1 ;
    	atmstep->nb_atomes++ ;
    	}
    }
 // Second check, ne devrait plus y avoir de problème !
 if (warnings>0)
  {
  for (i=0 ; i<atmstep->datas[idx[0]].size() ; i++)
    { if (atmstep->datas[idx[0]][i]!=i+1)
    	{DISP_Err("ERR : les ID n'augmentent pas par pas de 1 dans le LucDump. Cela aurait dû être corrigé à la boucle précédente !") ;}
    }
  }

 idx[0]=atmstep->find_idx(IDS("POSX")) ; idx[1]=atmstep->find_idx(IDS("POSY")) ; idx[2]=atmstep->find_idx(IDS("POSZ")) ;
 idx[3]=step.find_idx(IDS("CFFORCEX")) ; idx[4]=step.find_idx(IDS("CFFORCEY")) ; idx[5]=step.find_idx(IDS("CFFORCEZ")) ;
 idx[6]=step.find_idx(IDS("CFID1")) ; idx[7]=step.find_idx(IDS("CFID2")) ;
 idx[8]=step.find_idx(IDS("CFPOSX")) ; idx[9]=step.find_idx(IDS("CFPOSY")) ; idx[10]=step.find_idx(IDS("CFPOSZ")) ;
 idx[11]=step.find_idx(IDS("CFX")) ; idx[12]=step.find_idx(IDS("CFY")) ; idx[13]=step.find_idx(IDS("CFZ")) ;
 idx[14]=step.find_idx(IDS("CFMAG")) ;// idx[9]=step.find_idx(IDS("CFR")) ; idx[10]=step.find_idx(IDS("CFTHETA")) ; idx[11]=step.find_idx(IDS("CFPHI")) ;
 idx[15]=step.find_idx(IDS("CFPERIOD")) ;

 for (i=0 ; i<step.nb_atomes ; i++)
   {
   step.datas[idx[14]][i]=sqrt(pow(step.datas[idx[3]][i],2)+pow(step.datas[idx[4]][i],2)+pow(step.datas[idx[5]][i],2)) ;
   step.datas[idx[11]][i]=(atmstep->datas[idx[0]][(int)step.datas[idx[7]][i]-1])-(atmstep->datas[idx[0]][(int)step.datas[idx[6]][i]-1]) ;
   step.datas[idx[12]][i]=(atmstep->datas[idx[1]][(int)step.datas[idx[7]][i]-1])-(atmstep->datas[idx[1]][(int)step.datas[idx[6]][i]-1]) ;
   step.datas[idx[13]][i]=(atmstep->datas[idx[2]][(int)step.datas[idx[7]][i]-1])-(atmstep->datas[idx[2]][(int)step.datas[idx[6]][i]-1]) ;
   step.datas[idx[8]][i]=(atmstep->datas[idx[0]][(int)step.datas[idx[7]][i]-1]+atmstep->datas[idx[0]][(int)step.datas[idx[6]][i]-1])/2. ;
   step.datas[idx[9]][i]=(atmstep->datas[idx[1]][(int)step.datas[idx[7]][i]-1]+atmstep->datas[idx[1]][(int)step.datas[idx[6]][i]-1])/2. ;
   step.datas[idx[10]][i]=(atmstep->datas[idx[2]][(int)step.datas[idx[7]][i]-1]+atmstep->datas[idx[2]][(int)step.datas[idx[6]][i]-1])/2.;
   // TEST
   //printf("%.15f ", sqrt(step.datas[idx[11]][i]*step.datas[idx[11]][i]+step.datas[idx[12]][i]*step.datas[idx[12]][i])
   //		        -atmstep->datas[atmstep->find_idx(IDS("RAYON"))][(int)step.datas[idx[7]][i]-1]-atmstep->datas[atmstep->find_idx(IDS("RAYON"))][(int)step.datas[idx[6]][i]-1]) ;
   }

 // Si on veut prendre en compte les wallforce, il est temps de le faire !!!!
 if (actions["wallchainforce"].set)
  {
  Step * wallstep ; double chaineforce[3], position[3], norme ;
  int idxwall[7], added=0 ;
  op.wall_dump->check_timestep(ts) ; // Nécessite le même sampling entre tous les dump. On va supposer que le Ldump et le Walldump ont les même ts
  wallstep=&(op.wall_dump->steps[ts]) ;

  idxwall[0]=wallstep->find_idx(IDS("POSX")) ; idxwall[1]=wallstep->find_idx(IDS("POSY")) ; idxwall[2]=wallstep->find_idx(IDS("POSZ")) ;
  idxwall[3]=wallstep->find_idx(IDS("FORCEWALLX")) ; idxwall[4]=wallstep->find_idx(IDS("FORCEWALLY")) ; idxwall[5]=wallstep->find_idx(IDS("FORCEWALLZ")) ;
  idxwall[6]=wallstep->find_idx(IDS("ID")) ;

  for (i=0 ; i<wallstep->nb_atomes ; i++)
 	 {
	 if (wallstep->datas[idxwall[3]][i]==0 && wallstep->datas[idxwall[4]][i]==0 && wallstep->datas[idxwall[5]][i]==0) continue ;
	 norme=sqrt(pow(wallstep->datas[idxwall[0]][i]-actions["wallchainforce"]["xcyl"],2)+pow(wallstep->datas[idxwall[2]][i]-actions["wallchainforce"]["zcyl"],2)) ;
	 chaineforce[0]=2*actions.Cst["Radius"]*(wallstep->datas[idxwall[0]][i]-actions["wallchainforce"]["xcyl"])/norme ;
	 chaineforce[1]=0 ;
	 chaineforce[2]=2*actions.Cst["Radius"]*(wallstep->datas[idxwall[2]][i]-actions["wallchainforce"]["zcyl"])/norme ;
	 position[0]=wallstep->datas[idxwall[0]][i]-chaineforce[0]/2 ;
	 position[1]=wallstep->datas[idxwall[1]][i] ;
	 position[2]=wallstep->datas[idxwall[2]][i]-chaineforce[2]/2 ;

	 for (j=0 ; j<step.nb_idx ; j++)
	 	 {
		 if (step.idx_col[j]==IDS("CFID1")) { step.datas[j].push_back(wallstep->datas[idxwall[6]][i]) ; break ;}
		 else if (step.idx_col[j]==IDS("CFID2")) { step.datas[j].push_back(-1) ; break ; }// ID d'un mur mis à -1 par soucis de reconnaissance
		 else if (step.idx_col[j]==IDS("CFPERIOD")) { step.datas[j].push_back(0) ; break ;}
		 else if (step.idx_col[j]==IDS("CFFORCEX")) { step.datas[j].push_back(-wallstep->datas[idxwall[3]][i]) ; break ;}
		 else if (step.idx_col[j]==IDS("CFFORCEY")) { step.datas[j].push_back(-wallstep->datas[idxwall[4]][i]) ; break ;}
		 else if (step.idx_col[j]==IDS("CFFORCEZ")) { step.datas[j].push_back(-wallstep->datas[idxwall[5]][i]) ; break ;}
		 else if (step.idx_col[j]==IDS("CFMAG")) { step.datas[j].push_back(sqrt(pow(wallstep->datas[idxwall[3]][i],2)+pow(wallstep->datas[idxwall[4]][i],2)+pow(wallstep->datas[idxwall[5]][i],2))) ; break ;}
		 else if (step.idx_col[j]==IDS("CFX")) { step.datas[j].push_back(chaineforce[0]) ; break ;}
		 else if (step.idx_col[j]==IDS("CFY")) { step.datas[j].push_back(chaineforce[1]) ; break ;}
		 else if (step.idx_col[j]==IDS("CFZ")) { step.datas[j].push_back(chaineforce[2]) ; break ;}
		 else if (step.idx_col[j]==IDS("CFPOSX")) { step.datas[j].push_back(position[0]) ; break ;}
		 else if (step.idx_col[j]==IDS("CFPOSY")) { step.datas[j].push_back(position[1]) ; break ;}
		 else if (step.idx_col[j]==IDS("CFPOSZ")) { step.datas[j].push_back(position[2]) ; break ;}
		 else {DISP_Err("Erreur : lors de l'ajout de wallchainforce, un idxcol n'a pas pu être rempli") ; break ;}
	 	 }
	 added++ ;
 	 }
  // Check
  for (j=1 ; j<step.nb_idx ; j++)
 	 {
	 if (step.datas[j].size()!=step.datas[0].size()) DISP_Err("Erreur : l'ajout de wallchainforce ne s'est pas fait de la même manière sur tous les indexes.") ;
 	 }
  printf("Ajoutés : %d %d %d \n", step.nb_atomes, added, (int) step.datas[0].size()) ; fflush(stdout) ;
  if (step.nb_atomes+added != step.datas[0].size()) DISP_Err ("Erreur : l'ajout de wallforce a généré un problème ! ") ;
  step.nb_atomes+=added ;
  }
 }

 if (actions["wallchainforcenodump"].set) // ajout des forces tank particules d'après équilibre des forces
 {
  step.LCFforce_by_particle_v2(*atmstep,3) ;
 }

if (atmstep)
{
 idx[0]=atmstep->find_idx(IDS("IDMULTISPHERE")) ; idx[1]=step.find_idx(IDS("CFID1")) ; idx[2]=step.find_idx(IDS("CFID2")) ;
 if (idx[0]!=-1 && !actions["multisphere"].set) DISP_Warn("Attention, le champ id_mutisphere est défini dans le dump, mais aucune action n'est effectuée. utiliser --mutisphere peut aider\n") ;
 else if (idx[0]==-1 && actions["multisphere"].set) DISP_Warn("Etrange de definir --multisphere sans champ multisphere\n") ;
 else if (idx[0]!=-1 && actions["multisphere"].set)
 {//Suppression des chainforce internes aux multisphere
   int todelete=0 ;
   for (i=0 ; i<step.nb_atomes-todelete ; i++)
   {
     if ((atmstep->datas[idx[0]][(int)step.datas[idx[1]][i]-1])==-1 || (atmstep->datas[idx[0]][(int)step.datas[idx[2]][i]-1])==-1) continue ;
     if ((atmstep->datas[idx[0]][(int)step.datas[idx[1]][i]-1])==(atmstep->datas[idx[0]][(int)step.datas[idx[2]][i]-1]))
     {
      step.swap_atm(i,step.nb_atomes-todelete-1);
      i-- ;
      todelete++ ;
     }
   }
   step.nb_atomes-=todelete ;
   //DISP_Info("Nombre de chaînes internes aux multisphere supprimées :") ; printf("%d sur %d totales\n", todelete, step.nb_atomes+todelete) ;
 }
}

return 1 ;
}
//--------------------------------------------------------------------------
int Filter::wallforceatm (Step & step, struct Op op)
{
bool first ; int intersampling ; int samplingatm ;
int idx[8] ; Step * wstep ; int i ;

if (step.find_idx(IDS("FORCEWALLX")) != -1 || step.find_idx(IDS("FORCEWALLY")) != -1 || step.find_idx(IDS("FORCEWALLZ")) != -1)
{DISP_Info ("Le dump atomique contient déjà les forces avec le mur, le dump de wallforce n'a pas été utilisé") ; return 0 ; }

if (first==true)
 {
 op.wall_dump->check_timestep(0) ;
 intersampling=step.timestep/(double)op.wall_dump->steps[0].timestep ;
 samplingatm=step.timestep ;
 first=false ; i=0 ;
 }
else
  i=floor(step.timestep/samplingatm)*intersampling-1 ; // Le -1 donne une marge de maneuvre, c'est mieux ...

while (op.wall_dump->steps[i].timestep != step.timestep)
   {
   if (op.wall_dump->steps[i].timestep > step.timestep)
     {cout << "ERR : impossible de trouver des timestep correspondant entre le LucDump et le LcfDump" ; return -1 ; }
   i++ ;
   }
op.wall_dump->check_timestep(i) ;
wstep=&(op.wall_dump->steps[i]) ;

step.idx_col.push_back(IDS("FORCEWALLX")) ;
step.idx_col.push_back(IDS("FORCEWALLY")) ;
step.idx_col.push_back(IDS("FORCEWALLZ")) ;
step.nb_idx+=3 ;
step.datas.resize(step.idx_col.size()) ;
step.datas[step.find_idx(IDS("FORCEWALLX"))].resize(step.nb_atomes) ;
step.datas[step.find_idx(IDS("FORCEWALLY"))].resize(step.nb_atomes) ;
step.datas[step.find_idx(IDS("FORCEWALLZ"))].resize(step.nb_atomes) ;
idx[0]=step.nb_idx-3 ; idx[1]=step.nb_idx-2 ; idx[2]=step.nb_idx-1 ;
idx[3]=step.find_idx(IDS("ID") ) ;
idx[4]=wstep->find_idx(IDS("FORCEWALLX")) ;
idx[5]=wstep->find_idx(IDS("FORCEWALLY")) ;
idx[6]=wstep->find_idx(IDS("FORCEWALLZ")) ;
idx[7]=wstep->find_idx(IDS("ID")) ;

if (idx[0]==-1 || idx[1]==-1 || idx[2]==-1 || idx[3]==-1 || idx[4]==-1 || idx[5]==-1 || idx[6]==-1)
{DISP_Warn("Attention: wallchainforce non trouvées dans le walldump...") ; return 0 ; }

for (int i=0 ; i<step.nb_atomes ; i++)
{
 for (int j=0 ; j<wstep->nb_atomes ; j++)
 {
   if (wstep->datas[idx[7]][j]==step.datas[idx[3]][i])
   {
     step.datas[idx[0]][i]=wstep->datas[idx[4]][j] ;
     step.datas[idx[1]][i]=wstep->datas[idx[5]][j] ;
     step.datas[idx[2]][i]=wstep->datas[idx[6]][j] ;
     break ;
   }
 }
}

}

//=========================
int Filter::no_periodic_chain (Step & step, struct Op op)
{
int idx, i, j, tmp, old; static bool warn_once=true;
int idxs[6] ;
	// Suppression de toutes les chaînes periodiques
idx=step.find_idx(IDS("CFPERIOD")) ;

old=step.nb_atomes ;

i=0 ; j=step.nb_atomes-1 ;
while(i<j)
	{
	while (j>0 && step.datas[idx][j]==1) j-- ;
	while (i<step.nb_atomes && step.datas[idx][i]==0) i++ ;
	//if (i<j) { step.crush_atm(i,j) ; step.datas[idx][j]=1 ; }
	if (i<j) { step.swap_atm(i,j) ; }
	}

//for (j=0 ; j < step.nb_idx ; j++)
//	{step.datas[j].erase (step.datas[j].begin()+i, step.datas[j].end()) ; }
step.nb_atomes_supp=step.nb_atomes-i ;
step.nb_atomes=i ;

//printf("%d atomes supprimés sur %d initiaux.\n", old-step.nb_atomes, old) ;
if (step.nb_atomes_supp>0) step.has_periodic_chains=true ;
else step.has_periodic_chains=false ;

// Correction des chaines pour pas qu'elles soient trop fausses
idxs[0]=step.find_idx(IDS("CFPOSX")) ; idxs[1]=step.find_idx(IDS("CFPOSY")) ; idxs[2]=step.find_idx(IDS("CFPOSZ")) ;
idxs[3]=step.find_idx(IDS("CFX"))    ; idxs[4]=step.find_idx(IDS("CFY"))    ; idxs[5]=step.find_idx(IDS("CFZ")) ;
for (i=step.nb_atomes ; i<step.nb_atomes + step.nb_atomes_supp ; i++)
{
 if (warn_once==true)
 {warn_once=false ;
 DISP_Warn("Attention, le filtrage des chainforce periodiques n'a pas vraiment été testé. Ces chaines devraient idéalement ne pas être utilisées.\n") ;}

 if (step.datas[idxs[3]][i]>(step.box[0][1]-step.box[0][0])/2)
 {
   step.datas[idxs[3]][i]-=(step.box[0][1]-step.box[0][0]) ;
   step.datas[idxs[0]][i]-=(step.box[0][1]-step.box[0][0])/2 ;
 }
 if (step.datas[idxs[3]][i]<-(step.box[0][1]-step.box[0][0])/2)
 {
   step.datas[idxs[3]][i]+=(step.box[0][1]-step.box[0][0]) ;
   step.datas[idxs[0]][i]+=(step.box[0][1]-step.box[0][0])/2 ;
 }

 if (step.datas[idxs[4]][i]>(step.box[1][1]-step.box[1][0])/2)
 {
   step.datas[idxs[4]][i]-=(step.box[1][1]-step.box[1][0]) ;
   step.datas[idxs[1]][i]-=(step.box[1][1]-step.box[1][0])/2 ;
 }
 if (step.datas[idxs[4]][i]<-(step.box[1][1]-step.box[1][0])/2)
 {
   step.datas[idxs[4]][i]+=(step.box[1][1]-step.box[1][0]) ;
   step.datas[idxs[1]][i]+=(step.box[1][1]-step.box[1][0])/2 ;
 }

  if (step.datas[idxs[5]][i]>(step.box[2][1]-step.box[2][0])/2)
 {
   step.datas[idxs[5]][i]-=(step.box[2][1]-step.box[2][0]) ;
   step.datas[idxs[2]][i]-=(step.box[2][1]-step.box[2][0])/2 ;
 }
 if (step.datas[idxs[5]][i]<-(step.box[2][1]-step.box[2][0])/2)
 {
   step.datas[idxs[5]][i]+=(step.box[2][1]-step.box[2][0]) ;
   step.datas[idxs[2]][i]+=(step.box[2][1]-step.box[2][0])/2 ;
 }
}


//for (i=0 ; i<step.nb_atomes ; i++) {if (step.datas[idx][i]==1) { printf(".") ; fflush(stdout) ; }}
return 0 ;
}

//=========================
int Filter::moveto_xplusyplus (Step & step, struct Op op)
{
int i ;
int idx[3] ;
Matrix3x3 rot90, rotm90, rot180 ;
Vector axe (0.0,0.0,1.0) ;

idx[0]=step.find_idx(IDS("POSX")) ; idx[1]=step.find_idx(IDS("POSY")) ; idx[2]=step.find_idx(IDS("POSZ")) ;

rot90=Geometrie::get_rot_matrix(M_PI/2.0, axe) ;
rotm90=Geometrie::get_rot_matrix(-M_PI/2.0, axe) ;
rot180=Geometrie::get_rot_matrix(M_PI, axe) ;

for (i=0 ; i<step.nb_atomes ; i++)
	{
	if (step.datas[idx[0]][i]<0 && step.datas[idx[1]][i]>0) // rotation de -90deg
	   step.atm_rotate(rotm90, i) ;
	else if (step.datas[idx[0]][i]>0 && step.datas[idx[1]][i]<0) // rotation de +90deg
       step.atm_rotate(rot90,i) ;
	else if (step.datas[idx[0]][i]<0 && step.datas[idx[1]][i]<0) // rotation de 180deg
       step.atm_rotate(rot180,i) ;
	}
return 1 ;
}
//---------------------------------
int Filter::copypartiallink (Step & step, struct Op op)
{
int i ;
int idx[6] ;

idx[0]=step.find_idx(IDS("CFPOSX")) ; idx[1]=step.find_idx(IDS("CFPOSY")) ; idx[2]=step.find_idx(IDS("CFPOSZ")) ;
idx[3]=step.find_idx(IDS("CFX")) ; idx[4]=step.find_idx(IDS("CFY")) ; idx[5]=step.find_idx(IDS("CFZ")) ;
Matrix3x3 rot90, rotm90, rot180 ;
Vector axe (0.0,0.0,1.0) ;

idx[0]=step.find_idx(IDS("CFPOSX")) ; idx[1]=step.find_idx(IDS("CFPOSY")) ; idx[2]=step.find_idx(IDS("CFPOSZ")) ;
rot90=Geometrie::get_rot_matrix(M_PI/2.0, axe) ;
rotm90=Geometrie::get_rot_matrix(-M_PI/2.0, axe) ;
rot180=Geometrie::get_rot_matrix(M_PI, axe) ;

int nb_atomes_old=step.nb_atomes ;
for (i=0 ; i<nb_atomes_old ; i++)
	{
	if (step.datas[idx[1]][i]<actions.Cst["Radius"])
	   {
	   if (step.datas[idx[1]][i]-step.datas[idx[4]][i]/2<0 || step.datas[idx[1]][i]+step.datas[idx[4]][i]/2<0)
		   {step.copy_atm_end(i) ; step.nb_atomes++ ; step.atm_rotate(rot90, step.nb_atomes-1) ; }
	   }
	else if (step.datas[idx[0]][i]<actions.Cst["Radius"])
	   {
	   if (step.datas[idx[0]][i]-step.datas[idx[3]][i]/2<0 || step.datas[idx[0]][i]+step.datas[idx[3]][i]/2<0)
		   {step.copy_atm_end(i) ; step.nb_atomes++ ; step.atm_rotate(rotm90, step.nb_atomes-1) ; }
	   }
	}
return 1 ;
}
//--------------------------------
int Filter::make2D3D (Step & step, struct Op op)
{
 int idx[3] ;

 idx[0]=step.find_idx(IDS("POSX")) ; idx[1]=step.find_idx(IDS("POSY")) ; idx[2]=step.find_idx(IDS("POSZ")) ;
 if ((idx[0] !=-1 || idx[1] !=-1) && idx[2]==-1)
 {
    step.nb_idx ++ ;
    step.idx_col.push_back(IDS("POSZ")) ;
    step.datas.resize(step.datas.size()+1) ;
    step.datas[step.datas.size()-1].resize(step.nb_atomes,0) ;
 }

 idx[0]=step.find_idx(IDS("VX")) ; idx[1]=step.find_idx(IDS("VY")) ; idx[2]=step.find_idx(IDS("VZ")) ;
 if ((idx[0] !=-1 || idx[1] !=-1) && idx[2]==-1)
 {
    step.nb_idx ++ ;
    step.idx_col.push_back(IDS("VZ")) ;
    step.datas.resize(step.datas.size()+1) ;
    step.datas[step.datas.size()-1].resize(step.nb_atomes,0) ;
 }

 idx[0]=step.find_idx(IDS("FX")) ; idx[1]=step.find_idx(IDS("FY")) ; idx[2]=step.find_idx(IDS("FZ")) ;
 if ((idx[0] !=-1 || idx[1] !=-1) && idx[2]==-1)
 {
    step.nb_idx ++ ;
    step.idx_col.push_back(IDS("FZ")) ;
    step.datas.resize(step.datas.size()+1) ;
    step.datas[step.datas.size()-1].resize(step.nb_atomes,0) ;
 }

}

//--------------------------------------------------------------------
int Filter::multisphere_remove_atm (Step &step, struct Op op)
{
  step.multisphere->remove_atoms(step) ;
  return 0 ;
}
//---------------------------------------------------------------------
int Filter::quat2orient (Step &step, struct Op op)
{
 Quaternion<double> quat ;
 Matrix<double,3,3> rot, orient ;
 int idx[7] ;

 static Matrix<double,3,3> frame ; static bool first=true ;
 if (first)
 {
  if (!actions["superquadric"].set) DISP_Warn("Expecting --superquadric as arguments") ;
  frame(0,0)=actions["superquadric"]["semiaxisx"]*2 ; frame(1,1)=actions["superquadric"]["semiaxisy"]*2 ; frame(2,2)=actions["superquadric"]["semiaxisz"]*2 ; // *2 because we are giving the semi-axis
  frame(0,1)=frame(0,2)=frame(1,0)=frame(1,2)=frame(2,0)=frame(2,1)=0 ;
 }
 step.nb_idx +=9 ;
 step.idx_col.push_back(IDS("orientxx")) ; step.datas.resize(step.datas.size()+1) ; step.datas[step.datas.size()-1].resize(step.nb_atomes,0) ;
 step.idx_col.push_back(IDS("orientxy")) ; step.datas.resize(step.datas.size()+1) ; step.datas[step.datas.size()-1].resize(step.nb_atomes,0) ;
 step.idx_col.push_back(IDS("orientxz")) ; step.datas.resize(step.datas.size()+1) ; step.datas[step.datas.size()-1].resize(step.nb_atomes,0) ;
 step.idx_col.push_back(IDS("orientyx")) ; step.datas.resize(step.datas.size()+1) ; step.datas[step.datas.size()-1].resize(step.nb_atomes,0) ;
 step.idx_col.push_back(IDS("orientyy")) ; step.datas.resize(step.datas.size()+1) ; step.datas[step.datas.size()-1].resize(step.nb_atomes,0) ;
 step.idx_col.push_back(IDS("orientyz")) ; step.datas.resize(step.datas.size()+1) ; step.datas[step.datas.size()-1].resize(step.nb_atomes,0) ;
 step.idx_col.push_back(IDS("orientzx")) ; step.datas.resize(step.datas.size()+1) ; step.datas[step.datas.size()-1].resize(step.nb_atomes,0) ;
 step.idx_col.push_back(IDS("orientzy")) ; step.datas.resize(step.datas.size()+1) ; step.datas[step.datas.size()-1].resize(step.nb_atomes,0) ;
 step.idx_col.push_back(IDS("orientzz")) ; step.datas.resize(step.datas.size()+1) ; step.datas[step.datas.size()-1].resize(step.nb_atomes,0) ;

 idx[4]=step.find_idx(IDS("orientxx")) ;
 idx[0]=step.find_idx(IDS("QUATERNION1")) ; idx[1]=step.find_idx(IDS("QUATERNION2")) ; idx[2]=step.find_idx(IDS("QUATERNION3")) ; idx[3]=step.find_idx(IDS("QUATERNION4")) ;
 idx[5]=step.find_idx(IDS("TYPE")) ; idx[6]=step.find_idx(IDS("RAYON")) ;
 if (idx[5]!=-1 && idx[6]!=-1 && first) DISP_Info("If any atom has a type > 1, they are considered as spherical with the radius given by the dump\n") ;
 if (idx[5]!=-1 && idx[6]==-1 && first) DISP_Warn ("If some superquadric have type>1, things will break as the radius are not available\n") ; 
 if (idx[0]==-1 || idx[1]==-1 || idx[2]==-1 || idx[3]==-1 || idx[4]==-1) DISP_Err("IDX not found Filter::quat2orient\n") ;
 for (int i=0 ;i<step.nb_atomes ; i++)
 {
   if (idx[5]!=-1 && step.datas[idx[5]][i]>1)
   {
     for (int j=0 ; j<9 ; j++) step.datas[idx[4]+j][i]=(j/3==j%3)?2*step.datas[idx[6]][i]:0 ;
   }
   else
   {
    quat.w() = step.datas[idx[0]][i] ; quat.x() = step.datas[idx[1]][i] ; quat.y() = step.datas[idx[2]][i] ; quat.z() = step.datas[idx[3]][i] ;
    rot=quat.toRotationMatrix() ;
    orient=rot*frame*rot.transpose() ;
    for (int j=0 ; j<9 ; j++) step.datas[idx[4]+j][i]=orient(j/3,j%3) ;
   }
 }
first=false ;
}
//============================================================
int Filter::cluster (Step &step, struct Op op)
{
step.nb_idx++ ;
step.idx_col.push_back (IDS("group")) ; step.datas.resize(step.datas.size()+1) ; step.datas[step.datas.size()-1].resize(step.nb_atomes,0) ;

double maxdst=actions.Cst["Radius"]*2 ;
int gpcur=0 ;
int index=step.find_idx(IDS("GROUP")) ;
int idx[3] ;
idx[0]=step.find_idx(IDS("POSX")) ; idx[1]=step.find_idx(IDS("POSY")) ; idx[2]=step.find_idx(IDS("POSZ")) ;

auto distance = [&] (int i, int j) { return sqrt((step.datas[idx[0]][i]-step.datas[idx[0]][j])*(step.datas[idx[0]][i]-step.datas[idx[0]][j])+ (step.datas[idx[1]][i]-step.datas[idx[1]][j])*(step.datas[idx[1]][i]-step.datas[idx[1]][j]) + (step.datas[idx[2]][i]-step.datas[idx[2]][j])*(step.datas[idx[2]][i]-step.datas[idx[2]][j])) ; } ;

std::function<void(int,int)> paint ;
paint = [&](int i, int gp)
{
 for (int j=0 ; j<step.nb_atomes ; j++)
 {
  if (i==j) continue ;
  if (step.datas[index][j]==0 && distance(i,j) < maxdst)
  {
   step.datas[index][j]=gp ;
   paint(j, gp) ;
  }
 }
} ;

for (int i=0 ; i<step.nb_atomes ; i++)
{
 gpcur++ ;
 if (step.datas[index][i]==0)
 {
  step.datas[index][i]=gpcur ;
  paint(i, gpcur) ;
 }
}

return 0 ;
}




//=============================================================
int Filter::multisphere (Step &step, struct Op op)
{
  DISP_Err("REMOVED FROM IMPLEMENTATION:DO NO USE (filter::multisphere)") ;
/*  int i, j, k, l, n ;
  int idx[8] ; static vector < vector <int> > gps ; static bool first=true ; static int ngp=-1 ;
  static vector < bool > activegp ;
  static vector < Vector > pts, segments ;
  Vector t ; double box[6] ;

  Vector centroid ; int longest ; double maxlen ; static int type ;
  Vector vsph, null(0) ; int idmax ;
  idx[0]=step.find_idx(IDS("POSX")) ; idx[1]=step.find_idx(IDS("POSY")) ;  idx[2]=step.find_idx(IDS("POSZ")) ;
  idx[3]=step.find_idx(IDS("IDMULTISPHERE")) ; idx[4]=step.find_idx(IDS("ID")) ;

  step.idx_col.push_back(IDS("KX")) ; step.idx_col.push_back(IDS("KY")) ; step.idx_col.push_back(IDS("KZ")) ;
  step.nb_idx+=3 ;
  step.datas.resize(step.datas.size()+3) ;
  step.datas[step.datas.size()-1].resize(step.nb_atomes,0) ;
  step.datas[step.datas.size()-2].resize(step.nb_atomes,0) ;
  step.datas[step.datas.size()-3].resize(step.nb_atomes,0) ;
  idx[5]=step.datas.size()-3 ;
  idx[6]=idx[5]+1 ; idx[7]=idx[6]+1 ;

  if (actions["use-box"].set)
  {
   box[0] =  actions["use-box"]["box_xmin"] ;
   box[1] =  actions["use-box"]["box_xmax"] ;
   box[2] =  actions["use-box"]["box_ymin"] ;
   box[3] =  actions["use-box"]["box_ymax"] ;
   box[4] =  actions["use-box"]["box_zmin"] ;
   box[5] =  actions["use-box"]["box_zmax"] ;
  }
  else
  {box[0]=box[2]=box[4]=-std::numeric_limits <double>::infinity() ;
   box[1]=box[3]=box[5]= std::numeric_limits <double>::infinity()  ; }

  if (first)
  {
    first=false ;
    type=actions["multisphere"]["type"] ;
    for (i=0 ; i<step.nb_atomes ; i++)
    {
      if (step.datas[idx[3]][i]>=0)
      {
	if (step.datas[idx[3]][i]>ngp)
	{
	  gps.resize(step.datas[idx[3]][i]+1) ;
	  for (j=ngp+1 ; j<=step.datas[idx[3]][i] ; j++)
	    gps[j].push_back(0) ;
	  ngp=step.datas[idx[3]][i] ;
	}
	gps[(int)(step.datas[idx[3]][i])][0]++ ;
	gps[(int)(step.datas[idx[3]][i])].push_back((int)step.datas[idx[4]][i]) ;
      }
    }
  for (i=1, longest=0 ; i<=ngp ; i++) if (longest<gps[i][0]) longest=gps[i][0] ;
  pts.resize(longest, null) ;
  segments.resize(longest*(longest-1)/2, null) ;
  activegp.resize(ngp,true) ;
  }
 for (j=1 ; j<=ngp ; j++)
 {
  centroid=0 ;
  if (!activegp[j]) continue ;
  for (k=0 ; k<gps[j][0] ; k++)
  {
    if (step.datas[idx[4]][gps[j][k+1]-1]!=gps[j][k+1])
    {
      DISP_Err("Probleme in multisphere ici\n") ;
    }
    t.set(step.datas[idx[0]][gps[j][k+1]-1], step.datas[idx[1]][gps[j][k+1]-1], step.datas[idx[2]][gps[j][k+1]-1]);
    pts[k]=t ;
    if (t.isnan())
    {
      activegp[j]=false ; printf("Le groupe %d a été perdu. Atomes:", j) ; for (l=0 ; l<gps[j][0] ; l++) {printf("%d ", gps[j][l+1]-1) ; } printf("\n") ; break ;
    }
    centroid=centroid+pts[k] ;
  }
     if (!activegp[j]) continue ;
     centroid=centroid/4 ;
     if (centroid(1)<box[0] || centroid(1)>box[1] || centroid(2)<box[2] || centroid[2]>box[3] || centroid(3)<box[4] || centroid(3)>box[5]) { continue ; }
     for (k=0, n=0, maxlen=0, idmax=0 ; k<gps[j][0]-1 ; k++)
     {
       for (l=k+1 ; l<gps[j][0] ; l++, n++)
       {
	 segments[n]=pts[l]-pts[k] ;
	 if (maxlen<segments[n].norm())
	 {
	   maxlen=segments[n].norm() ;
	   idmax=n ;
	 }
       }
     }
     if (type==1) //Flat particles, have to do more
     {
       Vector crossp ;
       n=0 ;
       do
       {
	 crossp=segments[idmax].cross(segments[n]) ;
	 n++ ;
       } while (crossp.norm() < 0.000001 || crossp.isnan()) ;

       crossp=segments[idmax].norm()/crossp.norm()*crossp ;
       segments[idmax]=crossp ;

     }
     //else
     //  DISP_Err("Unknown particle shape for multisphere") ;
     vsph=Geometrie::cart2sph(segments[idmax]) ;
     if (type==0)
       {
	if (vsph(1)>actions.Cst["Radius"]*gps[j][0]*2)
	{segments[idmax](1)=NAN ; segments[idmax](2)=NAN ; segments[idmax](3)=NAN ;}
       }
     else if (type==1)
       {
	if (vsph(1)>actions.Cst["Radius"]*(floor(log2((gps[j][0]-1)/3.))*2)) {segments[idmax](1)=NAN ; segments[idmax](2)=NAN ; segments[idmax](3)=NAN ;}
       }
     for (k=0 ; k<gps[j][0] ; k++)
     {
	step.datas[idx[5]][gps[j][k+1]-1]=segments[idmax](1) ;
	step.datas[idx[6]][gps[j][k+1]-1]=segments[idmax](2) ;
	step.datas[idx[7]][gps[j][k+1]-1]=segments[idmax](3) ;
     }
  }*/
}



//-------------------------------------------------------------------
/*
void reorder_atm (int timestep)
{
  int idx ;

  // On va essayer d'utiliser un QuickSort pour le rangement ...

(".. En cours de rangement, ts %d -->", timestep) ; fflush(stdout) ;
  idx=find_idx(timestep, ID) ;
  quicksort(steps[timestep].datas, 0, steps[timestep].nb_atomes-1, idx, steps[timestep].nb_idx) ;
  printf("Done.    ") ; fflush(stdout) ;

  return ;
}*/
