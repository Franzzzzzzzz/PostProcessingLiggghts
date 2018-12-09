#include "Headers/Step.h"

// =====================================================
// Fonctions de la classe Step, type LStep================
// =====================================================
// =====================================================
// Fonctions orientées LStep ===========================
//======================================================
double Step::gravite (int atm)
{
 int idx ; 
 idx=check_idx(IDS("MASSE")) ;
 if (idx!=-1) return (-actions.Cst["G"]*datas[idx][atm]) ;  
 
 idx=check_idx(IDS("RAYON")) ;
 if (idx!=-1) 
 {
   if (actions["is2D"].set) return (-actions.Cst["G"]*actions.Cst["Rho"]*M_PI*datas[idx][atm]*datas[idx][atm]) ; 
   else return (-actions.Cst["G"]*actions.Cst["Rhograin"]*4./3.*M_PI*datas[idx][atm]*datas[idx][atm]*datas[idx][atm]) ; 
 }
 
 if (actions["is2D"].set) 
   return (-actions.Cst["G"]*actions.Cst["Rhograin"]*actions.Cst["Radius"]*actions.Cst["Radius"]) ; 
 else 
   return (-actions.Cst["G"]*actions.Cst["Rhograin"]*4./3.*M_PI*actions.Cst["Radius"]*actions.Cst["Radius"]*actions.Cst["Radius"]) ; 
}


void Step::Ldisp(void)
{
  cout << "Timestep : " << timestep << "\n" ;
  cout << "Nb atomes : " << nb_atomes << "\n" ;
  cout << "Box : " << box[0][0] << " "<< box[0][1] << " "<< box[1][0] << " "<< box[1][1] << " "<< box[2][0] << " "<< box[2][1] << "\n" ;
  cout << "Nb idx : " << nb_idx << "\n\n" ;
}
//-----------------------------
void Step::del_atm (long int atm)
{
 int i ;
 for (i=0 ; i<nb_idx ; i++)
 {if (idx_col[i]==IDS("POINTX") || idx_col[i]==IDS("POINTY") || idx_col[i]==IDS("POINTZ")) continue ;
  datas[i].erase(datas[i].begin()+atm) ;}
}
//-----------------------------
void Step::swap_atm (long int atm1, long int atm2)
{
 int i ; double tmp ;
 for (i=0 ; i<nb_idx ; i++)
 {/*if (idx_col[i]==IDS("POINTX") || idx_col[i]==IDS("POINTY") || idx_col[i]==IDS("POINTZ")) continue ;*/ // SHOULD NOT HAPPEN
  tmp=datas[i][atm1] ; datas[i][atm1]=datas[i][atm2] ; datas[i][atm2]=tmp ; }
}
//-----------------------------
void Step::crush_atm (long int atm1, long int atm2)
{
 int i ; double tmp ;
 for (i=0 ; i<nb_idx ; i++)
 {datas[i][atm1]=datas[i][atm2] ;}
}
//----------------------------
void Step::del_end_atms (long int nb)
{
 int i ;
 for (i=0 ; i<nb_idx ; i++)
 {if (idx_col[i]==IDS("POINTX") || idx_col[i]==IDS("POINTY") || idx_col[i]==IDS("POINTZ")) continue ;
  datas[i].erase(datas[i].end()-nb, datas[i].end()-1) ;}
}
//-----------------------------
void Step::copy_atm_end (long int nb)
{
 int i ;
 for (i=0 ; i<nb_idx ; i++)
 {if (idx_col[i]==IDS("POINTX") || idx_col[i]==IDS("POINTY") || idx_col[i]==IDS("POINTZ")) continue ;
  datas[i].push_back(datas[i][nb]) ; }
}
//-----------------------------
void Step::Lwrite_asVTK(ofstream &out)
{
int i , idx[3];
out << "# vtk DataFile Version 2.0\n" ;
out << "Generated by AVFF\nASCII\nDATASET POLYDATA\n" ;
out << "POINTS " << nb_atomes << " double\n" ;

idx[0]=find_idx(IDS("POSX")) ;
idx[1]=find_idx(IDS("POSY")) ;
idx[2]=find_idx(IDS("POSZ")) ;

for (i=0 ; i<nb_atomes ; i++)
{
 out << datas[idx[0]][i] << " " << datas[idx[1]][i] << " " << datas[idx[2]][i] << "\n" ;
}

out << "VERTICES " <<nb_atomes << " "<< nb_atomes*2 << "\n" ;

for (i=0 ; i<nb_atomes ; i++)
  {out << "1 " << i <<"\n" ;}

out << "POINT_DATA " << nb_atomes << "\n" ;
out << "SCALARS atom_type int 1\n" ;
out << "LOOKUP_TABLE default\n" ;
idx[0]=find_idx(IDS("TYPE")) ;
if (idx[0]==-1)
  {
  for (i=0 ; i<nb_atomes ; i++) out << "1 " ;
  }
else
  {
  for (i=0 ; i<nb_atomes ; i++) {out << datas[idx[0]][i] << " " ;}
  }

if (actions["w/speed"].set)
 {
 idx[0]=find_idx(IDS("VX")) ; idx[1]=find_idx(IDS("VY")) ; idx[2]=find_idx(IDS("VZ")) ;

 out << "VECTORS vitesse double\n" ;
 for (i=0 ; i<nb_atomes ; i++)
   {
   out <<datas[idx[0]][i] << " " << datas[idx[1]][i] << " " <<datas[idx[2]][i] << "\n" ;
   }
 }


if (actions["w/force"].set)
 {
 idx[0]=find_idx(IDS("FX")) ; idx[1]=find_idx(IDS("FY")) ; idx[2]=find_idx(IDS("FZ")) ;

 out << "VECTORS force double\n" ;
 for (i=0 ; i<nb_atomes ; i++)
   {
   out <<datas[idx[0]][i] << " " << datas[idx[1]][i] << " " <<datas[idx[2]][i] << "\n" ;
   }
 }

if (actions["w/id"].set)
 {
 idx[0]=find_idx(IDS("ID")) ;
 out << "SCALARS atom_id int 1\n" ;
 out << "LOOKUP_TABLE default\n" ;
 for (i=0 ; i<nb_atomes-1 ; i++)
   {
   out << datas[idx[0]][i] << " " ;
   }
 out << datas[idx[0]][nb_atomes-1] <<"\n" ;
 }

if (actions["w/rayon"].set)
 {
 idx[0]=find_idx(IDS("RAYON")) ;
 out << "SCALARS atom_rayon double 1\n" ;
 out << "LOOKUP_TABLE default\n" ;
 for (i=0 ; i<nb_atomes-1 ; i++)
   {
   out << datas[idx[0]][i] << " " ;
   }
 out << datas[idx[0]][nb_atomes-1] <<"\n" ;
 }

if (actions["w/masse"].set)
 {
 idx[0]=find_idx(IDS("MASSE")) ;
 out << "SCALARS atom_masse double 1\n" ;
 out << "LOOKUP_TABLE default\n" ;
 for (i=0 ; i<nb_atomes-1 ; i++)
   {
   out << datas[idx[0]][i] << " " ;
   }
 out << datas[idx[0]][nb_atomes-1] <<"\n" ;
 }

}
//--------------------------------------------------------------------------
void Step::Lwrite_asRESTART (ofstream &out)
{
int i ; int idx[6] ; int maxtype=0 ;
// Nb : ne fonctionne que pur les grains de type liggghts granular !!
DISP_Info("Ceci ne fonctionne que pour les grains de type granular\n") ;
cout << "Ecriture du timestep " << timestep << ".\n";

out << "Donnees atomiques generees par POSTPROCESSING. Permet une insertion d'atomes via la commande liggghts read_data.\n\n" ;
out << "#========================================================\n" ;
out << "# Header\n" ;
out << "#========================================================\n" ;
out << "   " << nb_atomes << " atoms\n" ;
out << "   1 atom types\n" ;
out << box[0][0] << " " << box[0][1] << " xlo xhi\n" ;
out << box[1][0] << " " << box[1][1] << " ylo yhi\n" ;
out << box[2][0] << " " << box[2][1] << " zlo zhi\n\n" ;

out <<"#======================================================\n" ;
out <<"# Atomes\n#======================================================\n\n" ;
out << "Atoms\n\n" ;

idx[0]=find_idx(IDS("ID")) ;
idx[1]=find_idx(IDS("TYPE")) ;
idx[2]=find_idx(IDS("POSX")) ; idx[3]=find_idx(IDS("POSY")) ; idx[4]=find_idx(IDS("POSZ")) ;
idx[5]=find_idx(IDS("RAYON")) ;
DISP_Info("La densité des atomes est fixée à 2500 !\n") ;

for (i=0 ; i<nb_atomes ; i++)
  {
  out << i+1 << " " << datas[idx[1]][i] << " " ;
  //out << datas[idx[0]][i] << " " << datas[idx[1]][i] <<" ";
  out << datas[idx[5]][i]*2 << " 2500 " ;
  out << datas[idx[2]][i] << " " << datas[idx[3]][i] << " " << datas[idx[4]][i] << "\n" ;
  if (datas[idx[1]][i]>maxtype) maxtype=datas[idx[1]][i] ;
  }
if (maxtype>1)
  DISP_Warn("Attention : il y a plus d'1 type d'atome, une erreur surviendra sans doute dans liggghts !") ;

out.close() ;
}
//------------------------------------------------------
int Step::Lmean_forces(double fr[])
{
int idx[3], i ;
fr[0]=fr[1]=fr[2]=0.0 ;

idx[0]=find_idx(IDS("FORCEWALLX")) ; idx[1]=find_idx(IDS("FORCEWALLY")) ; idx[2]=find_idx(IDS("FORCEWALLZ")) ;
if (idx[0]==-1||idx[1]==-1||idx[2]==-1)
   {DISP_Err("Les valeurs FORCEWALLX/FORCEWALLY/FORCEWALLZ sont necessaires avec --w/forcetot et un Ldump") ; return 1 ; }

for (i=0 ; i<nb_atomes ; i++)
   {
   fr[0]+=datas[idx[0]][i] ;
   fr[1]+=datas[idx[1]][i] ;
   fr[2]+=datas[idx[2]][i] ;
   }
return 0;
}
//-----------------------------------
int Step::wall_force(ofstream & out, double ** meanforces, int * meangrains)
{
double dtheta, angle, anglecentre, posx, posz, Fn ; int boites ;
int i, j, idx[6], fen ;
double **forces, factor, sigma ; int *grains ;
static bool info=true ;

// Initialisation des boites
boites=actions["wallforce-by-angle"]["nbbox_theta"] ;
dtheta=2*M_PI/(double)boites ;
posx=actions["wallforce-by-angle"]["xcyl"] ;
posz=actions["wallforce-by-angle"]["zcyl"] ;
// Initialisation des tableaux
forces=(double **)malloc(boites*sizeof(double *)) ;
grains=(int *)malloc(boites*sizeof(int)) ;
for (i=0 ; i<boites ; i++)
   {forces[i]=(double *)malloc(7*sizeof(double)) ; forces[i][0]=forces[i][1]=forces[i][2]=0 ; grains[i]=0 ;
    forces[i][3]=0 ; forces[i][4]=forces[i][5]=forces[i][6]=0 ;}

// Initialisation des index
idx[0]=find_idx(IDS("FORCEWALLX")) ; idx[1]=find_idx(IDS("FORCEWALLY")) ; idx[2]=find_idx(IDS("FORCEWALLZ")) ;
idx[3]=find_idx(IDS("POSX")) ; idx[4]=find_idx(IDS("POSY")) ; idx[5]=find_idx(IDS("POSZ")) ;
if (idx[0]==-1||idx[1]==-1||idx[2]==-1||idx[3]==-1||idx[4]==-1||idx[5]==-1)
   {printf("%d %d %d %d %d %d||", idx[0], idx[1], idx[2], idx[3], idx[4], idx[5] ) ;  
     DISP_Err("Les valeurs FORCEWALLX/FORCEWALLY/FORCEWALLZ sont necessaires avec --w/forcetot et un Ldump") ; return 1 ; }

if (actions["wallforce-by-angle"]["sigma"]==0) {fen=1; if (info) {DISP_Info("Utilisation d'une fenêtre créneau") ; info=false ; } } // Fenêtre créneau
else {sigma=actions["wallforce-by-angle"]["sigma"]/180*M_PI ; fen=2; if(info) {DISP_Info ("Utilisation d'une fenêtre gaussienne") ; info=false ; } } // Fenêtre gaussienne


// Parcours des atomes, ajout des forces aux angles adaptés
for (i=0 ; i<nb_atomes ; i++)
	{
	if (datas[idx[0]][i]!=0 || datas[idx[2]][i]!=0 || datas[idx[1]][i]!=0)
	 {
	 angle=Calcul::arctan(datas[idx[3]][i]-posx, datas[idx[5]][i]-posz) ;
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
	   forces[j][0]+=(factor*datas[idx[0]][i]) ; forces[j][1]+=(factor*datas[idx[1]][i]) ; forces[j][2]+=(factor*datas[idx[2]][i]) ;

	   // Force normale
	   Fn=datas[idx[0]][i]*cos(angle)+datas[idx[2]][i]*sin(angle) ;
	   forces[j][3]+=factor*Fn ;

	   // Force tangentielle
	   forces[j][4]+=factor*(datas[idx[0]][i]-Fn*cos(angle)) ;
	   forces[j][5]+=factor*(datas[idx[1]][i]) ;
	   forces[j][6]+=factor*(datas[idx[2]][i]-Fn*sin(angle)) ;

	   grains[j]+=factor ;
	   }

	 }
	}

// Ecriture des données et éventuellement ajout pour moyennage
for (i=0 ; i<boites ; i++)
   {
   if (actions["mean"].set)
      { for (j=0 ; j<7 ; j++) meanforces[i][j]+=forces[i][j] ;
        meangrains[i]+=grains[i] ; }
   out << (i*dtheta+dtheta/2) << " " << grains[i] << " "<< forces[i][0] << " " << forces[i][1] << " " << forces[i][2] << " " << forces[i][3] << " ";
   out << forces[i][4] << " " << forces[i][5] << " " << forces[i][6] << " " ;
   }

// Cleanning
for (i=0 ; i<boites ; i++)
   {free(forces[i]) ; }
free(forces) ;
free(grains) ;

return 0 ;
}
//-----------------------------------------------------
int Step::xray_projection (int dir, int w, int h, double **img, double * box)
{
int i, j, k ; int idx[3] ; 
double center[2], rayon, dst, jj, kk ; int bound[4] ; 

double dw, dh ; 
dw=(box[1]-box[0])/(double)(w) ; dh=(box[3]-box[2])/(double)(h) ;
if (abs(dw/dh-1)>0.05) {DISP_Info("Pixels de l'image assez éloignés d'être carrés ...") ; printf("dw=%g, dh=%g\n", dw, dh) ; } 

switch(dir) {
  case 0: idx[0]=find_idx(IDS("POSY")) ; idx[1]=find_idx(IDS("POSZ")) ; break ;
  case 1: idx[0]=find_idx(IDS("POSX")) ; idx[1]=find_idx(IDS("POSZ")) ; break ; 
  case 2: idx[0]=find_idx(IDS("POSX")) ; idx[1]=find_idx(IDS("POSY")) ; break ; 
  default: DISP_Err("Direction inconnue") ; return 1 ; 
}
idx[2]=find_idx(IDS("RAYON")) ; 
if (idx[0]==-1 || idx[1]==-1 || idx[2]==-1) DISP_Err("Data missing (maybe the radius?)\n") ; 

for (j=0 ; j<w ; j++) for (k=0 ; k<h ; k++) img[j][k]=0 ; 

for (i=0 ; i<nb_atomes ; i++)
 {
  // Recherche des bornes
  center[0]=round((datas[idx[0]][i]-box[0])/dw)-dw/2 ;
  center[1]=round((datas[idx[1]][i]-box[2])/dh)-dh/2 ; 
  bound[0]=fmax(round(center[0]-rayon/dw) - 1, 0) ; bound[1]=fmin(round(center[0]+rayon/dw) + 1, w-1) ; 
  bound[2]=fmax(round(center[1]-rayon/dh) - 1, 0) ; bound[3]=fmin(round(center[1]+rayon/dh) + 1, h-1) ; 
  
  // On revient en coordonnées réelles
  center[0]=datas[idx[0]][i] ;
  center[1]=datas[idx[1]][i] ; 
  rayon=datas[idx[2]][i] ; 
   
  for (j=bound[0] ; j<=bound[1] ; j++)
    for (k=bound[2] ; k<=bound[3] ; k++)
    {
      jj=j*dw+dw/2+box[0] ; kk=k*dh+dh/2+box[2] ;
      dst=sqrt( (center[0]-jj)*(center[0]-jj) + (center[1]-kk)*(center[1]-kk)) ; 
      if (dst<=rayon) {img[j][k]+=2*sqrt(rayon*rayon-dst*dst) ;}
    }
 }

for (j=0 ; j<w ; j++) for (k=0 ; k<h ; k++) 
{
  img[j][k]/=(box[5]-box[4]) ; 
  img[j][k]=exp(-2*img[j][k]) ;
}
return 0 ; 
}

//---------------------------------------------------------
double Step::GetVoronoi()
{
#ifndef VORONOI
DISP_Err("Program not compiled with voronoi support, no result computed.\n") ; 
return 0 ; 
#else

int id = int(actions["voronoi"]["id"]) ; 
int i ; 
vector <int> idx ; 

idx.push_back(find_idx(IDS("ID"))) ;
idx.push_back(find_idx(IDS("POSX"))) ;
idx.push_back(find_idx(IDS("POSY"))) ;
idx.push_back(find_idx(IDS("POSZ"))) ;
idx.push_back(find_idx(IDS("RAYON"))) ;

for (auto v:idx) if (v==-1) DISP_Warn("WARN: one of the expected index was not found in Step::GetVoronoi()\n") ; 

for (i=0 ; i<nb_atomes ; i++)
    if (datas[idx[0]][i]==id) 
        break ; 

id=i ; // Use the data index instead of the particle id to identify the particle 
double radius=datas[idx[4]][id] ;
double rd = actions.Cst["Radius"] ; 
double x=datas[idx[1]][id], y=datas[idx[2]][id], z=datas[idx[3]][id] ; 
double zmin, zmax ; 

if (actions["is2D"].set)
{
    if (x-radius-2*rd<box[0][0] || x+radius+2*rd>box[0][1] || y-radius-2*rd<box[1][0] || y+radius+2*rd>box[1][1])
        return -1 ; // Too close to a boundary ...
    zmin=-0.5 ; 
    zmax= 0.5 ; 
}
else
{
    if (x-2*radius<box[0][0] || x+2*radius>box[0][1] || y-2*radius<box[1][0] || y+2*radius>box[1][1] || z-2*radius<box[0][0] || z+2*radius>box[0][1])
        return -1 ; // Too close to a boundary ...
    zmin=z-5*radius ; 
    zmax=z+5*radius ; 
    
}

voro::container_poly cont(x-5*radius,x+5*radius,y-5*radius,y+5*radius, zmin, zmax,3,3,3,false,false,false,8);

for (int i=0,n=0 ; i<nb_atomes ; i++,n++)
{
    double dst=(x-datas[idx[1]][i])*(x-datas[idx[1]][i])+(y-datas[idx[2]][i])*(y-datas[idx[2]][i])+(z-datas[idx[3]][i])*(z-datas[idx[3]][i]) ; 
    if (dst<3*3*radius*radius)
    {
        cont.put(n,datas[idx[1]][i],datas[idx[2]][i],datas[idx[3]][i],datas[idx[4]][i]);
    }
}

voro::c_loop_all lp(cont); 
bool test ; 
for (lp.start(), test=true ; test ; test=lp.inc())
{
 if (lp.x()==x && lp.y()==y && lp.z()==z)
 {
     voro::voronoicell a ; 
     cont.compute_cell(a,lp) ;
     return (a.volume()) ; 
 }  
}
return -2 ; 
#endif
}



//------------------------------------------------------
int Step::Latm_rotate(Matrix3x3 & rot, int id)
{
int idx[3] ; 
Vector pt1 ;
idx[0]=find_idx(IDS("POSX")) ; idx[1]=find_idx(IDS("POSY")) ; idx[2]=find_idx(IDS("POSZ")) ;
if (idx[0]!=-1 && idx[1]!=-1 && idx[2]!=-1)
  {
  pt1.set(datas[idx[0]][id], datas[idx[1]][id], datas[idx[2]][id]) ;
  pt1=Geometrie::rotation(pt1, rot) ;
  datas[idx[0]][id]=pt1(1) ; datas[idx[1]][id]=pt1(2) ; datas[idx[2]][id]=pt1(3) ;
  }

idx[0]=find_idx(IDS("VX")) ; idx[1]=find_idx(IDS("VY")) ; idx[2]=find_idx(IDS("VZ")) ;
if (idx[0]!=-1 && idx[1]!=-1 && idx[2]!=-1)
   {
   pt1.set(datas[idx[0]][id], datas[idx[1]][id], datas[idx[2]][id]) ;
   pt1=Geometrie::rotation(pt1, rot) ;
   datas[idx[0]][id]=pt1(1) ; datas[idx[1]][id]=pt1(2) ; datas[idx[2]][id]=pt1(3) ;
   }

idx[0]=find_idx(IDS("FX")) ; idx[1]=find_idx(IDS("FY")) ; idx[2]=find_idx(IDS("FZ")) ;
if (idx[0]!=-1 && idx[1]!=-1 && idx[2]!=-1)
  {
  pt1.set(datas[idx[0]][id], datas[idx[1]][id], datas[idx[2]][id]) ;
  pt1=Geometrie::rotation(pt1, rot) ;
  datas[idx[0]][id]=pt1(1) ; datas[idx[1]][id]=pt1(2) ; datas[idx[2]][id]=pt1(3) ;
  }

return 1 ;
}
