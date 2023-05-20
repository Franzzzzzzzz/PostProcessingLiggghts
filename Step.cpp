#include "Headers/Step.h"

// =====================================================
// Fonctions générales de la classe Step================
// =====================================================
Step::Step() : nb_idx(0), nb_atomes(0), nb_atomes_supp(0), filtered(false) {triclinic[0]=triclinic[1]=triclinic[2]=0 ; }

//-----------------------------------------------
vector <int> Step::find_idx(vector<int>id, bool errcatch)
{
  vector<int> res ; bool bad ;
  for (auto i:id)
  {
    res.push_back(find_idx(i)) ;
    if (res[res.size()-1]==-1) bad ;
  }
  if (bad && errcatch) DISP_Err (_("An expected id could not be found, problems will happen.\n")) ;
  return (res) ;
}
//-----------------------------------------------
int Step::find_idx (int id)
{
int i ; static bool messages[5]={false,false,false,false,false} ;
// Recherche primaire de l'id
for (i=0 ; i<nb_idx ; i++)
  { /*printf("[%d]",idx_col[i]) ;  fflush(stdout) ; */if (idx_col[i]==id) {return i ;} }

// S'il n'existe pas, est-ce qu'on tente une reconstruction ?
if (!actions["norebuild"].set)
 {
  // On regarde s'il est possible de reconstruire la donnée demandée
  if (id==IDS("FORCEX") || id==IDS("FORCEY") || id==IDS("FORCEZ"))
  { if (find_idx(IDS("STRESSX"))!=-1 && find_idx(IDS("STRESSY"))!=-1 && find_idx(IDS("STRESSZ"))!=-1 )
      { if (messages[1]==false)
	{DISP_Info(_("Les id FORCE ont été cherchés mais non trouvés. Les id STRESS ont été trouvés. Les id de forces vont être reconstruits ... (indiqué une seule fois)\n")) ;
	 messages[1]=true ; }
	 buildtridata(1) ;
	 return (find_idx(id)) ;
      }
  }
  else if (id==IDS("NORMALEX") || id==IDS("NORMALEY") || id==IDS("NORMALEZ"))
  {
    if (messages[2]==false)
      {DISP_Info(_("Les id NORMALE ont été cherchés mais non trouvés. Ils vont être reconstruits ... (indiqué une seule fois)\n")) ;
       messages[2]=true ; }
      buildtridata(2) ;
      return (find_idx(id)) ;
  }
  else if (id==IDS("CENTREX") || id==IDS("CENTREY") || id==IDS("CENTREZ"))
  {
    if (messages[3]==false)
      {DISP_Info(_("Les id CENTRE ont été cherchés mais non trouvés. Ils vont être reconstruits ... (indiqué une seule fois)\n")) ;
       messages[3]=true ; }
      buildtridata(3) ;
      return (find_idx(id)) ;
  }
  else if ((id==IDS("FX") || id==IDS("FY") || id==IDS("FZ")) && (actions["chainforce"].set))
  {
    if (messages[2]==false)
    {DISP_Info(_("Les id F[XYZ] ont été cherchés mais non trouvés. Ils vont être reconstruits à partir du CFSTEP ... (indiqué une seule fois)\n")) ;
      messages[2]=true ; }
    otherstep->LCFforce_by_particle_v2(*this, 1) ;
    return (find_idx(id)) ;
  }
  else if ((id==IDS("FORCEWALLX") || id==IDS("FORCEWALLY") || id==IDS("FORCEWALLZ")) && (actions["chainforce"].set))
  {
    if (messages[2]==false)
    {DISP_Info(_("Les id FORCEWALL[XYZ] ont été cherchés mais non trouvés. Ils vont être reconstruits à partir du CFSTEP ... (indiqué une seule fois)\n")) ;
      messages[2]=true ; }
    otherstep->LCFforce_by_particle_v2(*this, 1) ;
    return (find_idx(id)) ;
  }
  else if (id==IDS("ATMPRESSURE") && (actions["chainforce"].set))
  {
    if (messages[2]==false)
    {DISP_Info(_("Les id ATMPRESSURE ont été cherchés mais non trouvés. Ils vont être reconstruits à partir du CFSTEP ... (indiqué une seule fois)\n")) ;
      messages[2]=true ; }
    otherstep->LCFpressure_by_particle(*this) ;
    return (find_idx(id)) ;
  }
 }
else // Pas de rebuild => on quitte
  {
  if (messages[0]==false)
	   {char erreur[100] ; sprintf(erreur,_("\nWARNING : l'id n° %d a été cherché mais pas trouvé. Enlever --norebuild peut aider. \n"), id) ; messages[0]=true ;
	   DISP_Warn(erreur) ; }
	return -1 ;
  }

return -1 ;
}
//--------------------------------------------------------
int Step::check_idx (int id) // Test l'existence des idx mais ne reconstruit pas par defaut (identique à Step.find_idx with actions["norebuild"]
{
  int i ;
  for (i=0 ; i<nb_idx ; i++)
  { /*printf("[%d]",idx_col[i]) ;  fflush(stdout) ; */if (idx_col[i]==id) {return i ;} }
  return -1 ;
}

//-------------------------------------------
bool Step::operator == (Step &step)
{
int i, idx;

if (idx_col.size()!=step.idx_col.size())
  return false ;

if (datas.size()!=step.datas.size())
  return false ;

for (i=0 ; i<(int)idx_col.size() ; i++)
  {
  idx=step.find_idx(idx_col[i]);
  if (datas[i].size()!=step.datas[idx].size()) return false ;
  if (datas[i]!=step.datas[idx]) return false ;
  }

return true ;
}
//------------------------------------------
Step & Step::operator = (const Step & step)
{
unsigned int i, j ;

Type=step.Type ;
posinfile=step.posinfile ;
nb_idx=step.nb_idx ;
idx_col=step.idx_col ;
datas=step.datas ;

for (i=0 ; i<datas.size() ; i++)
    {
    datas[i]=step.datas[i] ;
    }

timestep=step.timestep ;
for (i=0 ; i<3 ; i++) {for (j=0 ; j<2 ; j++) {box[i][j]=step.box[i][j];}}
nb_atomes=step.nb_atomes ;
nb_triangles=step.nb_triangles ;
nb_pts=step.nb_pts ;

cout << datas.size() ;
return *this ;
}


//-------------------------------------------
void Step::write_asDUMP (ofstream & out)
{
int i, j ;

// check_timestep(timestep) ; TODO : au dump de s'en occuper
out << "ITEM: TIMESTEP\n" << timestep << "\n" ;

if (Type==TL)
	{
	out << "ITEM: NUMBER OF ATOMS\n" << nb_atomes <<"\n" ;
	out << "ITEM: BOX BOUNDS\n" <<box[0][0]<<" "<< box[0][1]<<"\n"<<box[1][0]<<" "<<box[1][1]<<"\n" << box[2][0] << " " << box[2][1] << "\n" ;
	out << "ITEM: ATOMS " ;
	}
else if (Type==TCF)
	{
	out << "ITEM: NUMBER OF ENTRIES\n" << nb_atomes<< "\n" ;
	out << "ITEM: ENTRIES " ;
	}

bool flag=false ;
vector <string> aliases ;
for (i=0 ; i<nb_idx ; i++)
 {
 aliases=IDS.alias(idx_col[i]) ;
 if (aliases.size()==0) {DISP_Warn(_("WARN5 : index inconnu pour l'écriture des titres de colonnes\n")) ; continue ;}
 out << aliases[0] << " " ;
 }
 /*switch ()
   {
   case ID : out << "id " ; break ;
   case TYPE : out << "type " ; break ;
   case POSX : out << "x " ; break ;
   case POSY : out << "y " ; break ;
   case POSZ : out << "z " ; break ;
   case VX : out << "vx " ; break ;
   case VY : out << "vy " ; break ;
   case VZ : out << "vz " ; break ;
   case FX : out << "fx " ; break ;
   case FY : out << "fy " ; break ;
   case FZ : out << "fz " ; break ;
   case RAYON : out << "radius " ; break ;
   case MASSE : out << "mass " ; break ;
   case CFID1 : out << "c_cout[1] " ; break ;
   case CFID2 : out << "c_cout[2] " ; break ;
   case CFPERIOD : out << "c_cout[3] "; flag=true ; break ;
   case CFFORCEX : if (flag) out << "c_cout[4] "; else out << "c_cout[3] " ; break ;
   case CFFORCEY : if (flag) out << "c_cout[5] "; else out << "c_cout[4] " ; break ;
   case CFFORCEZ : if (flag) out << "c_cout[6] "; else out << "c_cout[5] " ; break ;
   case FORCEWALLX : out << "f_force_cyl[1] " ; break ;
   case FORCEWALLY : out << "f_force_cyl[2] " ; break ;
   case FORCEWALLZ : out << "f_force_cyl[3] " ; break ;
   default : DISP_Warn("WARN5 : index inconnu pour l'écriture des titres de colonnes\n") ; break ;
   }    */

out << "\n" ;

for (i=0 ; i<nb_atomes ; i++)
 {
 for (j=0 ; j<nb_idx ; j++)
   {
   out << datas[j][i] << " " ;
   }
 out << "\n" ;
 }
}
