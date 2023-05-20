#include "Headers/Step.h"

// ========================================================
// Fonctions de la classe Step, type CFStep================
// ========================================================
// ========================================================
// Fonctions orientées LCFStep ============================
//=========================================================
void Step::LCFdisp(void)
{
  cout << "Timestep : " << timestep << "\n" ;
  cout << "Nb entrées : " << nb_atomes << "\n" ;
  cout << "Nb idx : " << nb_idx << "\n\n" ;
}
//-------------------------------------------------------------
int Step::grain_force(double ** meanforces, double * meangrains, ofstream & out)
{
double dtheta, angle, anglecentre, Fn ; int boites ;
int i, j, idx[9], fen ;
double **forces, factor, sigma ; double *grains ;
static bool info=true ;

// Initialisation des boites
boites=actions["grainforce-by-angle"]["nbbox_theta"] ;
dtheta=2*M_PI/(double)boites ;
// Initialisation des tableaux
forces=(double **)malloc(boites*sizeof(double *)) ;
grains=(double *)malloc(boites*sizeof(double)) ;
for (i=0 ; i<boites ; i++)
   {forces[i]=(double *)malloc(3*sizeof(double)) ; forces[i][0]=forces[i][1]=forces[i][2]=0 ; grains[i]=0 ; }

// Initialisation des index
idx[0]=find_idx(IDS("CFFORCEX")) ; idx[1]=find_idx(IDS("CFFORCEY")) ; idx[2]=find_idx(IDS("CFFORCEZ")) ;
idx[3]=find_idx(IDS("CFX")) ; idx[4]=find_idx(IDS("CFY")) ; idx[5]=find_idx(IDS("CFZ")) ;
idx[6]=find_idx(IDS("CFID1")) ; idx[7]=find_idx(IDS("CFID2")) ; idx[8]=find_idx(IDS("CFPERIOD")) ;
if (idx[0]==-1||idx[1]==-1||idx[2]==-1||idx[3]==-1||idx[4]==-1||idx[5]==-1)
   {DISP_Err(_("Les valeurs CFFORCEX/Y/Z et CFX/Y/Z sont necessaires avec --grainforce_by_angle et un CFdump\n")) ; return 0 ; }

//if (has_periodic_chains==true) return 0 ;

if (actions["grainforce-by-angle"]["sigma"]==0) {fen=1; if (info) {DISP_Info(_("Utilisation d'une fenêtre créneau\n")) ; info=false ; } } // Fenêtre créneau
else {sigma=actions["grainforce-by-angle"]["sigma"]/180*M_PI ; fen=2; if(info) {DISP_Info (_("Utilisation d'une fenêtre gaussienne\n")) ; info=false ; } } // Fenêtre gaussienne

// Parcours des atomes, ajout des forces aux angles adaptés
int sgn, idatom ;
idatom=actions["grainforce-by-angle"]["id"] ; 
for (i=0 ; i<nb_atomes+nb_atomes_supp ; i++)
	{
	//if (datas[idx[8]][i]==1) return 0 ;
        if (datas[idx[7]][i]==idatom) sgn=-1;
	else if (datas[idx[6]][i]==idatom) sgn=1 ;  //Tout va bien
	else {DISP_Warn(_("Il n'y a pas le bon ID dans la chaines !\n")) ; continue ; }

	angle=Calcul::arctan(datas[idx[3]][i]*sgn, datas[idx[4]][i]*sgn) ;
	angle=Calcul::angle_0_2pi(angle) ;

	for (j=0 ; j<boites ; j++)
	   {
	   anglecentre=angle-(j*dtheta+dtheta/2) ;
	   if (anglecentre>M_PI) anglecentre=anglecentre-2*M_PI ;
	   if (anglecentre<-M_PI) anglecentre=anglecentre+2*M_PI ;
	   if (fen==1) Fonction::creneau1D (factor, anglecentre, dtheta) ;
	   else     Fonction::gaussienne1D (factor, anglecentre, sigma ) ;

	   factor*=dtheta ;
	   // Force totale
	   forces[j][0]+=(factor*datas[idx[0]][i]*sgn) ; forces[j][1]+=(factor*datas[idx[1]][i]*sgn) ; forces[j][2]+=(factor*datas[idx[2]][i]*sgn) ;

	   // Force normale
	   /*Fn=datas[idx[0]][i]*cos(angle)+datas[idx[2]][i]*sin(angle) ;
	   forces[j][3]+=factor*Fn ;

	   // Force tangentielle
	   forces[j][4]+=factor*(datas[idx[0]][i]-Fn*cos(angle)) ;
	   forces[j][5]+=factor*(datas[idx[1]][i]) ;
	   forces[j][6]+=factor*(datas[idx[2]][i]-Fn*sin(angle)) ;
       */
	   grains[j]+=factor ;
	   }

	 }

// Ecriture des données et éventuellement ajout pour moyennage
bool domean ;
domean=actions["mean"].set ; 
for (i=0 ; i<boites ; i++)
   {
    if (domean)
      { for (j=0 ; j<3 ; j++) meanforces[i][j]+=forces[i][j] ;
        meangrains[i]+=grains[i] ; }
    if (out.is_open())
    {
     out << (i*dtheta+dtheta/2) << " " << grains[i] << " "<< forces[i][0] << " " << forces[i][1] << " " << forces[i][2] << "\n" ; 
     //out << forces[i][4] << " " << forces[i][5] << " " << forces[i][6] << " " ;
    }
   }

// Cleanning
for (i=0 ; i<boites ; i++)
   {free(forces[i]) ; }
free(forces) ;
free(grains) ;

return 1 ;
}
//-----------------------------
int Step::LCFwrite_asVTK(ofstream & out, Step & step)
{
// C'est là qu'on rigole ...
int i , idxcf[12], idx[4];
int nb_liens ;

idxcf[0]=find_idx(IDS("CFID1")) ;
idxcf[1]=find_idx(IDS("CFID2")) ;
idxcf[2]=find_idx(IDS("CFFORCEX")) ;
idxcf[3]=find_idx(IDS("CFFORCEY")) ;
idxcf[4]=find_idx(IDS("CFFORCEZ")) ;
idxcf[5]=find_idx(IDS("CFMAG")) ;
idxcf[6]=find_idx(IDS("CFX")) ;
idxcf[7]=find_idx(IDS("CFY")) ;
idxcf[8]=find_idx(IDS("CFZ")) ;
idxcf[9]=find_idx(IDS("CFPOSX")) ;
idxcf[10]=find_idx(IDS("CFPOSY")) ;
idxcf[11]=find_idx(IDS("CFPOSZ")) ;
// On vérifie que les id augmentent linéairement par pas de 1 (pour simplifier ...) dans le LucDump
/*idx[0]=step.find_idx(IDS("ID")) ;
for (i=0 ; i<step.nb_atomes ; i++)
    {
    if (step.datas[idx[0]][i]!=i+1) {cout << "ERR : les ID n'augmentent pas par pas de 1 dans le LucDump. Programme à revoir ..." ; return -1 ; }
    }*/

// Si cutoff

if (actions["cutoff"].set)
   {
   nb_liens=floor((double)(actions["cutoff"]["cutoff"]/100.*nb_atomes)) ;
   }
else
   {
   nb_liens=nb_atomes ;
   }

// Header
out << "# vtk DataFile Version 2.0\n" ;
out << "Generated by AVFF. Chainforces !\nASCII\nDATASET POLYDATA\n" ;
out << "POINTS " << nb_liens*2 << " double\n" ;

idx[1]=step.find_idx(IDS("POSX")) ;
idx[2]=step.find_idx(IDS("POSY")) ;
idx[3]=step.find_idx(IDS("POSZ")) ;

/*for (i=nb_atomes-nb_liens ; i<nb_atomes ; i++)
{
 //if (actions.between)
 //   {
 //   if (step.datas[idx[3]][datas[idxcf[0]][i]-1] > actions.params["betweenlo"] && step.datas[idx[3]][datas[idxcf[1]][i]-1] > actions.params["betweenlo"]
 //      && step.datas[idx[3]][datas[idxcf[0]][i]-1] < actions.params["betweenhi"] && step.datas[idx[3]][datas[idxcf[1]][i]-1] < actions.params["betweenho"])
 //      {
 //      out << step.datas[idx[1]][datas[idxcf[0]][i]-1] << " " << step.datas[idx[2]][datas[idxcf[0]][i]-1] << " " << step.datas[idx[3]][datas[idxcf[0]][i]-1] << "\n" ;
 //      out << step.datas[idx[1]][datas[idxcf[1]][i]-1] << " " << step.datas[idx[2]][datas[idxcf[1]][i]-1] << " " << step.datas[idx[3]][datas[idxcf[1]][i]-1] << "\n" ;
 //      }
 //   }
 //else
 //   {
 //   }
out << step.datas[idx[1]][(int)datas[idxcf[0]][i]-1] << " " << step.datas[idx[2]][(int)datas[idxcf[0]][i]-1] << " " << step.datas[idx[3]][(int)datas[idxcf[0]][i]-1] << "\n" ;
out << step.datas[idx[1]][(int)datas[idxcf[1]][i]-1] << " " << step.datas[idx[2]][(int)datas[idxcf[1]][i]-1] << " " << step.datas[idx[3]][(int)datas[idxcf[1]][i]-1] << "\n" ;

}*/

for (i=nb_atomes-nb_liens ; i<nb_atomes ; i++)
	{
	out << datas[idxcf[9]][i]-datas[idxcf[6]][i]/2 << " " << datas[idxcf[10]][i]-datas[idxcf[7]][i]/2 << " " << datas[idxcf[11]][i]-datas[idxcf[8]][i]/2 << "\n" ;
	out << datas[idxcf[9]][i]+datas[idxcf[6]][i]/2 << " " << datas[idxcf[10]][i]+datas[idxcf[7]][i]/2 << " " << datas[idxcf[11]][i]+datas[idxcf[8]][i]/2 << "\n" ;
	}

out << "LINES " << nb_liens << " "<< nb_liens*3 << "\n" ;

for (i=nb_atomes-nb_liens ; i<nb_atomes ; i++)
  {out << "2 " << (i-nb_atomes+nb_liens)*2 << " " << (i-nb_atomes+nb_liens)*2+1 <<"\n" ; out.flush() ;}

out << "POINT_DATA " << nb_liens*2 << "\n" ;
out << "SCALARS magnitude double 1\n" ;
out << "LOOKUP_TABLE default\n" ;
for (i=nb_atomes-nb_liens ; i<nb_atomes ; i++)
 {out << datas[idxcf[5]][i] << " " << datas[idxcf[5]][i] << " " ;}

out << "SCALARS deltan double 1\n" ;
out << "LOOKUP_TABLE default\n" ;
double deltan ;
for (i=nb_atomes-nb_liens ; i<nb_atomes ; i++)
 {deltan=sqrt(datas[idxcf[6]][i]*datas[idxcf[6]][i]+datas[idxcf[7]][i]*datas[idxcf[7]][i]+datas[idxcf[8]][i]*datas[idxcf[8]][i]) ;
	out << 0.0015-deltan << " " << 0.0015-deltan << " " ;}


out << "\nVECTORS force double\n" ;
for (i=0 ; i<nb_atomes ; i++)
  {
  out <<datas[idxcf[2]][i] << " " << datas[idxcf[3]][i] << " " <<datas[idxcf[4]][i] << "\n" ;
  out <<datas[idxcf[2]][i] << " " << datas[idxcf[3]][i] << " " <<datas[idxcf[4]][i] << "\n" ;
  }
out << "\nVECTORS r double\n" ;
for (i=0 ; i<nb_atomes ; i++)
  {
  out <<datas[idxcf[6]][i] << " " << datas[idxcf[7]][i] << " " <<datas[idxcf[8]][i] << "\n" ;
  out <<datas[idxcf[6]][i] << " " << datas[idxcf[7]][i] << " " <<datas[idxcf[8]][i] << "\n" ;
  }

/*out << "SCALARS pz double 1\n" ;
out << "LOOKUP_TABLE default\n" ;
for (i=nb_atomes-nb_liens ; i<nb_atomes ; i++)
 {out << step.datas[step.find_idx(IDS("POSZ"))][(int)datas[idxcf[0]][i]-1] << " " << step.datas[step.find_idx(IDS("POSZ"))][(int)datas[idxcf[1]][i]-1]  << " " ;}


idxcf[5]=find_idx(IDS("POSZ")) ;
out << "CELL_DATA " << nb_liens << "\n" ;
out << "SCALARS POSZ double 1\n" ;
out << "LOOKUP_TABLE default\n" ;
for (i=nb_atomes-nb_liens ; i<nb_atomes ; i++)
 {out << datas[idxcf[5]][i] <<" " ;}

idxcf[5]=find_idx(IDS("CFR")) ;
out << "SCALARS CFR double 1\n" ;
out << "LOOKUP_TABLE default\n" ;
for (i=nb_atomes-nb_liens ; i<nb_atomes ; i++)
 {out << datas[idxcf[5]][i] <<" " ;}

idxcf[5]=find_idx(IDS("CFTHETA")) ;
out << "SCALARS CFTHETA double 1\n" ;
out << "LOOKUP_TABLE default\n" ;
for (i=nb_atomes-nb_liens ; i<nb_atomes ; i++)
 {out << datas[idxcf[5]][i] <<" " ;}

idxcf[5]=find_idx(IDS("CFPHI")) ;
out << "SCALARS CFPHI double 1\n" ;
out << "LOOKUP_TABLE default\n" ;
for (i=nb_atomes-nb_liens ; i<nb_atomes ; i++)
 {out << datas[idxcf[5]][i] <<" " ;}
*/
return 1 ;
}
//-----------------------------------------
int Step::LCFatm_rotate(Matrix3x3 & rot, int id)
{
int idx[3] ;
Vector pt1(3) ;

/*idx[0]=find_idx(IDS("POSX")) ; idx[1]=find_idx(IDS("POSY")) ; idx[2]=find_idx(IDS("POSZ")) ;
if (idx[0]!=-1 && idx[1]!=-1 && idx[2]!=-1)
   {
   pt1=Convert::toVect(datas[idx[0]][id], datas[idx[1]][id], datas[idx[2]][id]) ;
   pt1=Geometrie::rotation(pt1, rot) ;
   datas[idx[0]][id]=pt1(1,1) ; datas[idx[1]][id]=pt1(2,1) ; datas[idx[2]][id]=pt1(3,1) ;
   }

idx[0]=find_idx(IDS("CFX")) ; idx[1]=find_idx(IDS("CFY")) ; idx[2]=find_idx(IDS("CFZ")) ;
if (idx[0]!=-1 && idx[1]!=-1 && idx[2]!=-1)
   {
   pt1=Convert::toVect(datas[idx[0]][id], datas[idx[1]][id], datas[idx[2]][id]) ;
   pt1=Geometrie::rotation(pt1, rot) ;
   datas[idx[0]][id]=pt1(1,1) ; datas[idx[1]][id]=pt1(2,1) ; datas[idx[2]][id]=pt1(3,1) ;
   }*/

idx[0]=find_idx(IDS("CFFORCEX")) ; idx[1]=find_idx(IDS("CFFORCEY")) ; idx[2]=find_idx(IDS("CFFORCEZ")) ;
if (idx[0]!=-1 && idx[1]!=-1 && idx[2]!=-1)
   {
   pt1.set(datas[idx[0]][id], datas[idx[1]][id], datas[idx[2]][id]) ;
   pt1=Geometrie::rotation(pt1, rot) ;
   datas[idx[0]][id]=pt1(1) ; datas[idx[1]][id]=pt1(2) ; datas[idx[2]][id]=pt1(3) ;
   }
return 1;
}

//-----------------------------------------------------
int Step::LCFforce_by_tank (Step &lstep, double * resultat) // Calcule la force exercée par le mur du tank cylindrique sur les particules.
{
 int i, j, idx[11], zidx ; double dz, gravite, angle ; 
 int nbslice=actions["forcetank"]["slices"] ; double b[2] ; 
 
 gravite=4/3.*M_PI*actions.Cst["Radius"]*actions.Cst["Radius"]*actions.Cst["Radius"]*actions.Cst["Rhograin"]*actions.Cst["G"] ;
 if (actions["use-box"].set)
 {
   if (actions["is2D"].set) 
   {
     if (isinf(actions["use-box"]["box_ymin"])) b[0]=lstep.box[1][0] ;
     if (isinf(actions["use-box"]["box_ymax"])) b[1]=lstep.box[1][1] ;
   }
   else
   {
     if (isinf(actions["use-box"]["box_ymin"])) b[0]=lstep.box[2][0] ;
     if (isinf(actions["use-box"]["box_ymax"])) b[1]=lstep.box[2][1] ;
   }
 }
 else
 {
   if (actions["is2D"].set) {b[0]=lstep.box[1][0] ; b[1]=lstep.box[1][1] ; }
   else 		    {b[0]=lstep.box[2][0] ; b[1]=lstep.box[2][1] ; }
 }
 dz=(b[1]-b[0])/(double)(nbslice) ; 
 
  // Calcul des forces exercées par les murs sur les particules (tous les murs) ; 
 LCFforce_by_particle_v2 (lstep, 2) ; 
 idx[0]=lstep.find_idx(IDS("FORCEWALLX")) ; idx[1]=lstep.find_idx(IDS("FORCEWALLY")) ; idx[2]=lstep.find_idx(IDS("FORCEWALLZ")) ; 
 idx[5]=lstep.find_idx(IDS("POSX")) ; idx[6]=lstep.find_idx(IDS("POSY")) ; idx[7]=lstep.find_idx(IDS("POSZ")) ;
 //idx[8]=lstep.find_idx(IDS("FX")) ; idx[9]=lstep.find_idx(IDS("FY")) ; idx[10]=lstep.find_idx(IDS("FZ")) ;
 
 for (i=0 ; i<3 ; i++)
   for (j=0 ; j<nbslice ; j++)
     resultat[i*nbslice+j]=0 ; 
 
 for (i=0 ; i<lstep.nb_atomes ; i++ )
 {
    if (actions["is2D"].set) zidx=lstep.datas[idx[6]][i]/dz ;
    else zidx=lstep.datas[idx[7]][i]/dz ;
    
    if (zidx<0 || zidx>=nbslice) {continue ; }
    //xyz force
    //resultat[0*nbslice+zidx]+=(lstep.datas[idx[8]][i]-lstep.datas[idx[0]][i]) ;
    //resultat[1*nbslice+zidx]+=(lstep.datas[idx[9]][i]-lstep.datas[idx[1]][i]) ;
    //resultat[2*nbslice+zidx]+=(lstep.datas[idx[10]][i]-lstep.datas[idx[2]][i]+gravite) ;
    
    // r theta z force
    angle=atan2(lstep.datas[idx[6]][i],lstep.datas[idx[5]][i]) ; 
    resultat[0*nbslice+zidx]+=lstep.datas[idx[0]][i] ; 
    resultat[1*nbslice+zidx]+=lstep.datas[idx[1]][i] ; 
    resultat[2*nbslice+zidx]+=lstep.datas[idx[2]][i] ; 
    //resultat[0*nbslice+zidx]+=((lstep.datas[idx[8]][i]-lstep.datas[idx[0]][i])*cos(angle)+(lstep.datas[idx[9]][i]-lstep.datas[idx[1]][i])*sin(angle)) ;
    //resultat[1*nbslice+zidx]+=(-(lstep.datas[idx[8]][i]-lstep.datas[idx[0]][i])*sin(angle)+(lstep.datas[idx[9]][i]-lstep.datas[idx[1]][i])*cos(angle)) ;
    //resultat[2*nbslice+zidx]+=(lstep.datas[idx[10]][i]-lstep.datas[idx[2]][i]+gravite) ;
 } 
 return 0 ; 
}

//-----------------------------------------------
int Step::LCFcouple (Step & lstep, Vector &torque, Cylindre C) 
{
 int i, idx[12], natm=0 ; 
 Vector x, contact, force ; 
 double dst, gravite ; 
 
 torque(1)=torque(2)=torque(3)=0;
 
 LCFforce_by_particle_v2 (lstep, 2) ; 
 DISP_Err(_("Il faut verifier cette fonction avant de l'utiliser")) ; 
 // Ecarts de force (nul pour la plupart (ou epsilon), seuls les grains en contact avec une paroi doivent montrer quelque chose.
 idx[0]=lstep.find_idx(IDS("FORCEWALLX")) ; idx[1]=lstep.find_idx(IDS("FORCEWALLY")) ; idx[2]=lstep.find_idx(IDS("FORCEWALLZ")) ;
 idx[8]=lstep.find_idx(IDS("FX")) ; idx[9]=lstep.find_idx(IDS("FY")) ; idx[10]=lstep.find_idx(IDS("FZ")) ;
 gravite=4/3.*M_PI*actions.Cst["Radius"]*actions.Cst["Radius"]*actions.Cst["Radius"]*actions.Cst["Rhograin"]*actions.Cst["G"] ; 
 for (i=0 ; i<lstep.nb_atomes ; i++ )
 {
    lstep.datas[idx[0]][i]=lstep.datas[idx[8]][i]-lstep.datas[idx[0]][i] ;
    lstep.datas[idx[1]][i]=lstep.datas[idx[9]][i]-lstep.datas[idx[1]][i] ;
    lstep.datas[idx[2]][i]=lstep.datas[idx[10]][i]-lstep.datas[idx[2]][i]+gravite ;
 }
 
 // Grain en contact avec le cylindre ??
 idx[5]=lstep.find_idx(IDS("POSX")) ; idx[6]=lstep.find_idx(IDS("POSY")) ; idx[7]=lstep.find_idx(IDS("POSZ")) ;
 for (i=0 ; i<lstep.nb_atomes ; i++)
 {   
  x(1)=lstep.datas[idx[5]][i]; x(2)=lstep.datas[idx[6]][i]; x(3)=lstep.datas[idx[7]][i]; 
  dst=Geometrie::distance_cylindre(C, x) ;
  if (dst < 0 ) DISP_Warn ("Erreur : des grains à l'intérieur du cylindre ????\n") ; 
  if (Geometrie::distance_cylindre(C, x) < 4*actions.Cst["Radius"]) // *4 pour avoir une marge de sécurité
  {
    //lstep.datas[idx[0]][i]=lstep.datas[idx[8]][i]-lstep.datas[idx[0]][i] ;
    //lstep.datas[idx[1]][i]=lstep.datas[idx[9]][i]-lstep.datas[idx[1]][i] ;
    //lstep.datas[idx[2]][i]=lstep.datas[idx[10]][i]-lstep.datas[idx[2]][i] ;
    force(1)=lstep.datas[idx[0]][i] ; force(2)=lstep.datas[idx[1]][i] ; force(3)=lstep.datas[idx[2]][i] ;  
    contact=Geometrie::contact_cylindre(C, x) ;  
    torque=torque+force.cross(contact-C.centre) ; 
    natm++ ; 
  } 
 }
 return natm ;  
}

//------------------------------------------------------
/*int Step::LCFforce_by_particle (Step &lstep, int type)
{
  int idx[8], i ; 
  
  // Types : 0 : reconstruit les FORCEWALL à partir des 
  
  if (type==0)
  { lstep.nb_idx+=3 ; 
    lstep.idx_col.push_back(IDS("FORCEWALLX")) ; lstep.idx_col.push_back(IDS("FORCEWALLY")) ; lstep.idx_col.push_back(IDS("FORCEWALLZ")) ;     
    lstep.datas.resize(lstep.idx_col.size()) ; 
    idx[0]=lstep.find_idx(IDS("FORCEWALLX")) ; idx[1]=lstep.find_idx(IDS("FORCEWALLY")) ; idx[2]=lstep.find_idx(IDS("FORCEWALLZ")) ; 
    lstep.datas[idx[0]].resize(lstep.nb_atomes, 0) ; lstep.datas[idx[1]].resize(lstep.nb_atomes, 0) ; lstep.datas[idx[2]].resize(lstep.nb_atomes, 0) ; 
  }
  else if (type==1)
  {
    lstep.nb_idx+=3 ; 
    lstep.idx_col.push_back(IDS("FX")) ; lstep.idx_col.push_back(IDS("FY")) ; lstep.idx_col.push_back(IDS("FZ")) ; 
    lstep.datas.resize(lstep.idx_col.size()) ; 
    idx[0]=lstep.find_idx(IDS("FX")) ; idx[1]=lstep.find_idx(IDS("FY")) ; idx[2]=lstep.find_idx(IDS("FZ")) ; 
    lstep.datas[idx[0]].resize(lstep.nb_atomes, 0) ; lstep.datas[idx[1]].resize(lstep.nb_atomes, 0) ; lstep.datas[idx[2]].resize(lstep.nb_atomes, 0) ; 
  }
  else 
    DISP_Warn("Type inconnu LCFForce_by_particle") ; 
  
  idx[3]=find_idx(IDS("CFID1")) ; idx[4]=find_idx(IDS("CFID2")) ; 
  idx[5]=find_idx(IDS("CFFORCEX")) ; idx[6]=find_idx(IDS("CFFORCEY")) ; idx[7]=find_idx(IDS("CFFORCEZ")) ;
  for (i=0 ; i<nb_atomes ; i++)
  {
    lstep.datas[idx[0]][datas[idx[3]][i]-1]+=datas[idx[5]][i] ;
    lstep.datas[idx[1]][datas[idx[3]][i]-1]+=datas[idx[6]][i] ;
    lstep.datas[idx[2]][datas[idx[3]][i]-1]+=datas[idx[7]][i] ;
  
  //if (datas[idx[3]][i]==228073) {printf("[%f %f %f]
  
    lstep.datas[idx[0]][datas[idx[4]][i]-1]+=(-datas[idx[5]][i]) ;
    lstep.datas[idx[1]][datas[idx[4]][i]-1]+=(-datas[idx[6]][i]) ;
    lstep.datas[idx[2]][datas[idx[4]][i]-1]+=(-datas[idx[7]][i]) ;
  }
return 1 ;   
}*/
// Better function LCFforce_by_particle
int Step::LCFforce_by_particle_v2 (Step &lstep, int type)
{
 bool rebuild ; vector <double> ftotx, ftoty, ftotz ; int i, j ; double fatm[3] ; double tankzmin, tankzmax, tankrayon ;
 double dst ; Vector x, ur, pos, cf ; 
 int idxatm[15], idxcf[8] ; int added ; double rayon ; 
 // Type 1 : reconstruit les FX FY FZ des particules à partir des chaines de force fournies (les FX FY FZ ne doivent pas exister).
 // Type 2 : reconstruit les FORCEWALL des particules à partir des chaines de force fournies.
 //          2bis (mais ça reste 2 ...) si FX FY FZ existent, on les utilisent pour la force due au murs
 // Type 3 : reconstruit les chaines de forces dues au mur (donc il faut d'abord reconstruire les FORCEWALL si besoin
 
 // Calcul des forces totales par particules
 
 ftotx.resize(lstep.nb_atomes, 0) ; ftoty.resize(lstep.nb_atomes, 0) ; ftotz.resize(lstep.nb_atomes, 0) ; 
 idxatm[0]=lstep.find_idx(IDS("ID")) ; 
 idxcf[0]=find_idx(IDS("CFID1")) ; idxcf[1]=find_idx(IDS("CFID2")) ; 
 idxcf[2]=find_idx(IDS("CFFORCEX")) ; idxcf[3]=find_idx(IDS("CFFORCEY")) ; idxcf[4]=find_idx(IDS("CFFORCEZ")) ;
 for (i=0 ; i<nb_atomes ; i++)
 {
   if (lstep.datas[idxatm[0]][datas[idxcf[0]][i]-1]!=datas[idxcf[0]][i] || lstep.datas[idxatm[0]][datas[idxcf[1]][i]-1]!=datas[idxcf[1]][i]) 
   {DISP_Err(_("Incohérence dans les ID et les indices de tableau (Step::LCFforce_by_particle_v2)")) ; 
     printf("%d %d\n",idxatm[0], lstep.idx_col[idxatm[0]]) ; }
   
   ftotx[datas[idxcf[0]][i]-1]+=  datas[idxcf[2]][i]  ; ftoty[datas[idxcf[0]][i]-1]+=  datas[idxcf[3]][i]  ; ftotz[datas[idxcf[0]][i]-1]+=  datas[idxcf[4]][i]  ;
   ftotx[datas[idxcf[1]][i]-1]+=(-datas[idxcf[2]][i]) ; ftoty[datas[idxcf[1]][i]-1]+=(-datas[idxcf[3]][i]) ; ftotz[datas[idxcf[1]][i]-1]+=(-datas[idxcf[4]][i]) ;
 }

 switch (type) {
   case 1 : //DISP_Warn("Unchecked function") ; 
	    if (lstep.check_idx(IDS("FX"))!=-1 || lstep.check_idx(IDS("FY"))!=-1 || lstep.check_idx(IDS("FZ"))!=-1) DISP_Warn (_("Les FX/FY/FZ existent déjà, pas de reconstruction\n")) ; 
	    else
	    { lstep.idx_col.push_back(IDS("FX")) ; lstep.idx_col.push_back(IDS("FY")) ; lstep.idx_col.push_back(IDS("FZ")) ; 
	      lstep.nb_idx+=3 ; 
	      idxatm[1]=lstep.check_idx(IDS("FX")) ; idxatm[2]=lstep.check_idx(IDS("FY")) ; idxatm[3]=lstep.check_idx(IDS("FZ")) ; 
	      lstep.datas.resize(lstep.nb_idx) ; 
	      lstep.datas[idxatm[1]].resize(lstep.nb_atomes,0) ; lstep.datas[idxatm[2]].resize(lstep.nb_atomes,0) ; lstep.datas[idxatm[3]].resize(lstep.nb_atomes,0) ;
	      
	      for (i=0 ; i<lstep.nb_atomes ; i++) 
	      {
		lstep.datas[idxatm[1]][i]=ftotx[i] ; 
		lstep.datas[idxatm[2]][i]=ftoty[i] ; 
		lstep.datas[idxatm[3]][i]=ftotz[i]+lstep.gravite(i) ; 
	      }
	    }
     break ; 
     
   case 2 : //DISP_Warn("Unchecked function") ; 
	    if (lstep.check_idx(IDS("FORCEWALLX"))!=-1 || lstep.check_idx(IDS("FORCEWALLY"))!=-1 || lstep.check_idx(IDS("FORCEWALLZ"))!=-1) DISP_Warn (_("Les FORCEWALLX/Y/Z existent déjà, pas de reconstruction\n")) ; 
   	    else
	    { lstep.idx_col.push_back(IDS("FORCEWALLX")) ; lstep.idx_col.push_back(IDS("FORCEWALLY")) ; lstep.idx_col.push_back(IDS("FORCEWALLZ")) ; 
	      lstep.nb_idx+=3 ;
	      idxatm[4]=lstep.check_idx(IDS("FORCEWALLX")) ; idxatm[5]=lstep.check_idx(IDS("FORCEWALLY")) ; idxatm[6]=lstep.check_idx(IDS("FORCEWALLZ")) ; 
	      lstep.datas.resize(lstep.nb_idx) ; 
	      lstep.datas[idxatm[4]].resize(lstep.nb_atomes,0) ; lstep.datas[idxatm[5]].resize(lstep.nb_atomes,0) ; lstep.datas[idxatm[6]].resize(lstep.nb_atomes,0) ;
	      //printf("%d %f ", i, lstep.datas[idxatm[4]].at(10)) ; fflush(stdout) ;
	      idxatm[1]=lstep.check_idx(IDS("FX")) ; idxatm[2]=lstep.check_idx(IDS("FY")) ; idxatm[3]=lstep.check_idx(IDS("FZ")) ; 
	       
	      for (i=0 ; i<lstep.nb_atomes ; i++) 
	      {
		if (idxatm[1]==-1 || idxatm[2]==-1 || idxatm[3]==-1) fatm[0]=fatm[1]=fatm[2]=0 ; 
		else {fatm[0]=lstep.datas[idxatm[1]][i] ; fatm[1]=lstep.datas[idxatm[2]][i] ; fatm[2]=lstep.datas[idxatm[3]][i];}
		lstep.datas[idxatm[4]][i]=fatm[0]-ftotx[i] ; 
		lstep.datas[idxatm[5]][i]=fatm[1]-ftoty[i] ; 
		lstep.datas[idxatm[6]][i]=fatm[2]-ftotz[i]-lstep.gravite(i) ; 
		//printf("%g %g %g %g|", fatm[0], fatm[1], ftotx[i],ftoty[i]) ; 
	      }
	    }
     break ; 
     
   case 3 : // Avant de recréer les chainforces vraiment, il faut calculer les forcewall, même si on les met dans la variable temporaire
	    if (lstep.check_idx(IDS("FORCEWALLX"))!=-1 || lstep.check_idx(IDS("FORCEWALLY"))!=-1 || lstep.check_idx(IDS("FORCEWALLZ"))!=-1) 
	      DISP_Warn (_("Les FORCEWALLX/Y/Z existent, mais ne seront pas utilisées (bizarre)\n")) ; 
	    
	    idxatm[1]=lstep.check_idx(IDS("FX")) ; idxatm[2]=lstep.check_idx(IDS("FY")) ; idxatm[3]=lstep.check_idx(IDS("FZ")) ; 
	    idxatm[7]=lstep.find_idx(IDS("POSX")) ; idxatm[8]=lstep.find_idx(IDS("POSY")) ; idxatm[9]=lstep.find_idx(IDS("POSZ")) ; 
	    idxatm[10]=lstep.check_idx(IDS("RAYON")) ; 
	    
	    added=0 ; 
	    for (i=0 ; i<lstep.nb_atomes ; i++) 
	    {
	      if (idxatm[1]==-1 || idxatm[2]==-1 || idxatm[3]==-1) fatm[0]=fatm[1]=fatm[2]=0 ; 
	      else {fatm[0]=lstep.datas[idxatm[1]][i] ; fatm[1]=lstep.datas[idxatm[2]][i] ; fatm[2]=lstep.datas[idxatm[3]][i];}

	      ftotx[i]=fatm[0]-ftotx[i] ; 
	      ftoty[i]=fatm[1]-ftoty[i] ; 
	      ftotz[i]=fatm[2]-ftotz[i]-lstep.gravite(i) ; 
	      
	      if (ftotx[i]==0 && ftoty[i]==0 && ftotz[i]==0) continue ; 

	      // On va supposer que on ne prend en compte que le mur du tank cylindrique (parce que c'est ce dont j'ai besoin).
	      tankrayon=0.04 ; 
	      dst=sqrt(lstep.datas[idxatm[7]][i]*lstep.datas[idxatm[7]][i]+lstep.datas[idxatm[8]][i]*lstep.datas[idxatm[8]][i]) ;
	      
	      x(1)=lstep.datas[idxatm[7]][i] ; x(2)=lstep.datas[idxatm[8]][i] ; x(3)=lstep.datas[idxatm[9]][i] ;
	      
	      ur=x ; ur(3)=0 ; ur=ur/(ur.norm()) ; 
	      if (idxatm[10]==-1) rayon=actions.Cst["Radius"] ; 
	      else rayon=lstep.datas[idxatm[10]][i] ; 
	      pos=x+rayon*ur ; cf=2*rayon*ur ; 

	      if (dst>tankrayon-2*rayon) // x2 pour être tranquille avec le fait qu'on a des facettes
	      {
	       
	       for (j=0 ; j<nb_idx ; j++)
	 	 {
		 if      (idx_col[j]==IDS("CFID1")) { datas[j].push_back(lstep.datas[idxatm[0]][i]) ; }
		 else if (idx_col[j]==IDS("CFID2")) { datas[j].push_back(-1) ; }// ID d'un mur mis à -1 par soucis de reconnaissance
		 else if (idx_col[j]==IDS("CFPERIOD")) { datas[j].push_back(0) ; }
		 else if (idx_col[j]==IDS("CFFORCEX")) { datas[j].push_back(ftotx[i]) ; }
		 else if (idx_col[j]==IDS("CFFORCEY")) { datas[j].push_back(ftoty[i]) ; }
		 else if (idx_col[j]==IDS("CFFORCEZ")) { datas[j].push_back(ftotz[i]) ; }
		 else if (idx_col[j]==IDS("CFMAG")) { datas[j].push_back(sqrt(ftotx[i]*ftotx[i]+ftoty[i]*ftoty[i]+ftotz[i]*ftotz[i])) ; }
		 else if (idx_col[j]==IDS("CFX")) { datas[j].push_back(cf(1)) ; }
		 else if (idx_col[j]==IDS("CFY")) { datas[j].push_back(cf(2)) ; }
		 else if (idx_col[j]==IDS("CFZ")) { datas[j].push_back(cf(3)) ; }
		 else if (idx_col[j]==IDS("CFPOSX")) { datas[j].push_back(pos(1)) ; }
		 else if (idx_col[j]==IDS("CFPOSY")) { datas[j].push_back(pos(2)) ; }
		 else if (idx_col[j]==IDS("CFPOSZ")) { datas[j].push_back(pos(3)) ; }
		 else {DISP_Err(_("Erreur : lors de l'ajout de wallchainforce, un idxcol n'a pas pu être rempli")) ; }
	 	 }
		added++ ;
	      }
	    }
	    
	    // Check
	    for (j=1 ; j<nb_idx ; j++)
	      if (datas[j].size()!=datas[0].size()) 
		DISP_Err(_("Erreur : l'ajout de wallchainforce ne s'est pas fait de la même manière sur tous les indexes.")) ;
	      
	    printf("Ajoutés : %d %d %d \n", nb_atomes, added, (int) datas[0].size()) ; fflush(stdout) ;
	    if (nb_atomes+added != datas[0].size()) DISP_Err ("Erreur : l'ajout de wallforce a généré un problème ! ") ;
	    nb_atomes+=added ;
	    break ; // car oui, on était toujours dans un case ;-)
     
    default : DISP_Warn(_("Type inconnu dans Step::LCFforce_by_particle_v2\n")) ; printf("%d ", type) ; 
	    break ; 
 }
return 0 ; 
}



//-----------------------------------------------------------------
int Step::LCFpressure_by_particle (Step &lstep)
{
int idx[11], i ; 
static bool message=true, infogive=true ; 
Vector f, r ; 

if (message) {DISP_Warn(_("En phase de testing ...")) ; message=false ; }

lstep.nb_idx+=1 ; 
lstep.idx_col.push_back(IDS("ATMPRESSURE")) ;     
lstep.datas.resize(lstep.idx_col.size()) ; 
idx[0]=lstep.find_idx(IDS("ATMPRESSURE")) ; 
lstep.datas[idx[0]].resize(lstep.nb_atomes, 0) ;

idx[3]=find_idx(IDS("CFID1")) ; idx[4]=find_idx(IDS("CFID2")) ; 
idx[5]=find_idx(IDS("CFFORCEX")) ; idx[6]=find_idx(IDS("CFFORCEY")) ; idx[7]=find_idx(IDS("CFFORCEZ")) ;
idx[8]=find_idx(IDS("CFX")) ; idx[9]=find_idx(IDS("CFY")) ; idx[10]=find_idx(IDS("CFZ")) ;
for (i=0 ; i<nb_atomes ; i++)
{
 r(1)=datas[idx[8]][i] ; r(2)=datas[idx[9]][i] ; r(3)=datas[idx[10]][i] ; 
 f(1)=datas[idx[5]][i] ; f(2)=datas[idx[6]][i] ; f(3)=datas[idx[7 ]][i] ;  
 r.normalise() ; 
 lstep.datas[idx[0]][datas[idx[3]][i]-1]+= -(r.dot(f)) ;
 lstep.datas[idx[0]][datas[idx[4]][i]-1]+= -(r.dot(f)) ; // En fait à la fois le signe de r et de f changent !
}

for (i=0 ; i<lstep.nb_atomes ; i++)
{
    lstep.datas[idx[0]][i]/=(2*M_PI*actions.Cst["Radius"]*actions.Cst["Radius"]) ; 
}
return 1 ;   
}


//-------------------------------------------------------------------
map <string, double> Step::LCFenergy (Step &lstep) 
{
 int idx[14] ; static bool first=true ; 
 map <string, double> res ; 
 double nrj_k=0, nrj_r=0, nrj_kn=0, nrj_ktvirtual=0, nrj_gn=0, nrj_gt=0, nrj_mu=0 ;
 double nrj_springkn=0, nrj_ktgrain=0 ;  

 idx[0]=lstep.find_idx(IDS("VX")) ; idx[1]=lstep.find_idx(IDS("VY")) ; idx[2]=lstep.find_idx(IDS("VZ")) ; 
 if (idx[0]==-1 || idx[1]==-1 || idx[2]==-1) DISP_Err(_("Missing velocity value in the ldump\n")) ; 
 idx[3]=lstep.find_idx(IDS("OMEGAX")) ; idx[4]=lstep.find_idx(IDS("OMEGAY")) ; idx[5]=lstep.find_idx(IDS("OMEGAZ")) ; 
 if (idx[3]==-1 || idx[1]==-1 || idx[2]==-1) DISP_Err(_("Missing angular velocity value in the ldump\n")) ;  
 idx[6]=lstep.find_idx(IDS("MASSE")) ; idx[7]=lstep.find_idx(IDS("RAYON")) ; 
 
 //idx[8]=find_idx(IDS("c_cout[7]")) ;   idx[9]=find_idx(IDS("c_cout[10]"))   ;  // deltas
 //idx[10]=find_idx(IDS("c_cout[8]")) ;  idx[11]=find_idx(IDS("c_cout[9]")) ; // kn gamman
 //idx[12]=find_idx(IDS("c_cout[11]")) ;  idx[13]=find_idx(IDS("c_cout[12]")) ; //kt gammat
 //idx[14]=find_idx(IDS("c_cout[13]"))   ;                                        //mu
 
 /*idx[8]=find_idx(IDS("CFNRJKN")) ;
 idx[9]=find_idx(IDS("CFNRJGN")) ;
 idx[10]=find_idx(IDS("CFNRJktgrain")) ;
 idx[11]=find_idx(IDS("NRJktvirtual")) ;
 idx[12]=find_idx(IDS("NRJgt")) ;
 idx[13]=find_idx(IDS("NRJmu")) ;*/
 
 
 //for (int i=0 ; i<14 ; i++) if (idx[i]==-1) {DISP_Err("Missing data:") ; printf("%d \n", i) ; fflush(stdout) ; return (res ); } 
 
 if (first) {DISP_Info(_("Using moment of inertia of a sphere")) ; first=false ; }
 for (int i=0 ; i<lstep.nb_atomes ; i++)
 {
  if (std::isnan(lstep.datas[idx[6]][i])) continue ; 
  nrj_k += 0.5*lstep.datas[idx[6]][i]*(lstep.datas[idx[0]][i]*lstep.datas[idx[0]][i]+lstep.datas[idx[1]][i]*lstep.datas[idx[1]][i]+lstep.datas[idx[2]][i]*lstep.datas[idx[2]][i]) ; 
  nrj_r += 1/5. * lstep.datas[idx[6]][i] * lstep.datas[idx[7]][i] * lstep.datas[idx[7]][i] * (lstep.datas[idx[3]][i]*lstep.datas[idx[3]][i]+lstep.datas[idx[4]][i]*lstep.datas[idx[4]][i]+lstep.datas[idx[5]][i]*lstep.datas[idx[5]][i]) ;
 }
 
 
 
 /*for (int i=0 ; i<nb_atomes+nb_atomes_supp ; i++)
 {
  nrj_kn += datas[idx[8]][i] ; 
  nrj_ktgrain += datas[idx[10]][i] ; 
  nrj_ktvirtual += datas[idx[11]][i] ; 
  nrj_gn += datas[idx[9]][i] ; 
  nrj_gt += datas[idx[12]][i] ; 
  nrj_mu += datas[idx[13]][i] ; 
 }*/
 
 res["NRJkinetic"] = nrj_k ;
 res["NRJrotation"]=nrj_r ;
 res["NRJkn"]= nrj_kn ; 
 res["NRJktgrain"]= nrj_ktgrain ;
 res["NRJktvirtual"]= nrj_ktvirtual ;
 res["NRJgn"]= nrj_gn ;
 res["NRJgt"]= nrj_gt ;
 res["NRJmu"]= nrj_mu ;
 return res ; 
}

