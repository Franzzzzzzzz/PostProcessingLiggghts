#include "Headers/Step.h"

// =====================================================
// Fonctions de la classe Step, type FStep================
// =====================================================
// =====================================================
// Fonctions orientées FStep ===========================
//======================================================
void Step::Fdisp(void)
{
  cout << "Nb triangles : " << nb_triangles << "\n" ;
  cout << "Nb points (si alloué) : " << nb_pts << "\n" ;
  printf( "Nb idx : %d\n", nb_idx) ;
}
//-------------------------------------------------
void Step::Fwrite_asVTK(ofstream &out)
{
int i ;
int idx[4] ;
out << "# vtk DataFile Version 2.0\nLIGGGHTS+AVFF mesh/gran/VTK export\nASCII\nDATASET POLYDATA\n" ;

out.precision(6) ;

if (find_idx(IDS("POINTX"))>-1)
 {
 out << "POINTS " << nb_pts << " double\n" ;
 idx[0]=find_idx(IDS("POINTX")) ; idx[1]=find_idx(IDS("POINTY")) ; idx[2]=find_idx(IDS("POINTZ")) ;
 for (i=0 ; i<nb_pts ; i++)
   {
   out << fixed << epsilon(datas[idx[0]][i]) << " " << epsilon(datas[idx[1]][i]) << " " << epsilon(datas[idx[2]][i]) << "\n" ;
   }
 }

if (find_idx(IDS("POLY1"))>-1)
 {
 out << "POLYGONS " << nb_triangles << " " << nb_triangles*4 << "\n" ;
 idx[0]=find_idx(IDS("POLY1")) ; idx[1]=find_idx(IDS("POLY2")) ; idx[2]=find_idx(IDS("POLY3")) ; idx[3]=find_idx(IDS("POLY4")) ;
 out.precision(0) ;
 for (i=0 ; i<nb_triangles ; i++)
   {
   out << fixed << datas[idx[0]][i] << " " << datas[idx[1]][i] << " " << datas[idx[2]][i] << " " << datas[idx[3]][i] << "\n" ;
   }
 }

out << "CELL_DATA " << nb_triangles << "\n" ;
out.precision(6) ;
if (find_idx(IDS("PRESSURE"))>-1)
 {
 out << "SCALARS pressure double 1\nLOOKUP_TABLE default\n" ;
 idx[0]=find_idx(IDS("PRESSURE")) ;
 for (i=0 ; i<nb_triangles ; i++)
  {out << fixed << epsilon(datas[idx[0]][i]) << "\n" ;}
 }

if (find_idx(IDS("SHEARSTRESS"))>-1)
 {
 out<<"SCALARS shearstress double 1\nLOOKUP_TABLE default\n" ;
 idx[0]=find_idx(IDS("SHEARSTRESS")) ;
 for (i=0 ; i<nb_triangles ; i++)
  {out << fixed << epsilon(datas[idx[0]][i]) << "\n" ;}
 }

if (find_idx(IDS("FORCEX"))>-1)
 {
 out << "VECTORS forcesTri double\n" ;
 idx[0]=find_idx(IDS("FORCEX")) ; idx[1]=find_idx(IDS("FORCEY")) ; idx[2]=find_idx(IDS("FORCEZ")) ;
 for (i=0 ; i<nb_triangles ; i++)
  {out << fixed << epsilon(datas[idx[0]][i]) << " " << epsilon(datas[idx[1]][i]) << " " << epsilon(datas[idx[2]][i]) << "\n" ; }
 }

if (find_idx(IDS("NORMALEX"))>(-1))
 {
out << "VECTORS normales double\n" ;
 idx[0]=find_idx(IDS("NORMALEX")) ; idx[1]=find_idx(IDS("NORMALEY")) ; idx[2]=find_idx(IDS("NORMALEZ")) ;
 for (i=0 ; i<nb_triangles ; i++)
  { out << fixed << epsilon(datas[idx[0]][i]) << " " << epsilon(datas[idx[1]][i]) << " "<< epsilon(datas[idx[2]][i]) << "\n" ; }
 }

if (find_idx(IDS("CENTREX"))>(-1))
 {
out << "VECTORS centres double\n" ;
 idx[0]=find_idx(IDS("CENTREX")) ; idx[1]=find_idx(IDS("CENTREY")) ; idx[2]=find_idx(IDS("CENTREZ")) ;
 for (i=0 ; i<nb_triangles ; i++)
  { out << fixed << epsilon(datas[idx[0]][i]) << " " << epsilon(datas[idx[1]][i]) << " "<< epsilon(datas[idx[2]][i]) << "\n" ; }
 }
}
//----------------------------------------
int Step::Fmean_forces(double fr[])
{
int i ; int idx[3] ;
fr[0]=fr[1]=fr[2]=0 ;
idx[0]=find_idx(IDS("FORCEX")) ; idx[1]=find_idx(IDS("FORCEY")) ; idx[2]=find_idx(IDS("FORCEZ")) ;
if (idx[0]<0 || idx[1]<0 || idx[2]<0) {DISP_Err("Erreur : pas de donnees de FORCE dans le VTK\n") ; fflush(stdout); }

for (i=0 ; i<nb_triangles ; i++)
  {
  fr[0]+=datas[idx[0]][i] ;
  fr[1]+=datas[idx[1]][i] ;
  fr[2]+=datas[idx[2]][i] ;
  }
return 1;
}
//-----------------------------------
int Step::buildtridata (int type)
{
int idx[3],k ;
double aire ;
Vector normale, centres ;

switch (type) {
// Calcul des forces à partir des stresses
case 1 : 	nb_idx+=3 ; datas.resize(nb_idx, vector<double>(1,0)) ;	idx_col.resize(nb_idx,0) ;
			idx_col[nb_idx-3]=IDS("FORCEX") ; datas[nb_idx-3].resize(nb_triangles,0.) ;
			idx_col[nb_idx-2]=IDS("FORCEY") ; datas[nb_idx-2].resize(nb_triangles,0.) ;
			idx_col[nb_idx-1]=IDS("FORCEZ") ; datas[nb_idx-1].resize(nb_triangles,0.) ;

			idx[0]=find_idx(IDS("STRESSX")) ; idx[1]=find_idx(IDS("STRESSY")) ; idx[2]=find_idx(IDS("STRESSZ")) ;
			for (k=0; k<nb_triangles ; k++)
         	 {
        	 aire=get_tri_surface(k) ;
        	 datas[nb_idx-3][k]=datas[idx[0]][k]*aire ;
        	 datas[nb_idx-2][k]=datas[idx[1]][k]*aire ;
        	 datas[nb_idx-1][k]=datas[idx[2]][k]*aire ;
         	 }
			break ;
// Calcul des normales des triangles
case 2 : 	nb_idx+=3 ; datas.resize(nb_idx, vector<double>(1,0)) ;	idx_col.resize(nb_idx,0) ;
		 	idx_col[nb_idx-3]=IDS("NORMALEX") ; datas[nb_idx-3].resize(nb_triangles,0.) ;
		 	idx_col[nb_idx-2]=IDS("NORMALEY") ; datas[nb_idx-2].resize(nb_triangles,0.) ;
		 	idx_col[nb_idx-1]=IDS("NORMALEZ") ; datas[nb_idx-1].resize(nb_triangles,0.) ;

		 	for (k=0; k<nb_triangles ; k++)
		 		{
		 		normale=(get_tri_normal(k)) ;
		 		datas[nb_idx-3][k]=normale(1) ;	datas[nb_idx-2][k]=normale(2) ;	datas[nb_idx-1][k]=normale(3) ;
		 		}
		 	break ;
// Calcul des centres des triangles
case 3 : 	nb_idx+=3 ; datas.resize(nb_idx, vector<double>(1,0)) ;	idx_col.resize(nb_idx,0) ;
		 	idx_col[nb_idx-3]=IDS("CENTREX") ; datas[nb_idx-3].resize(nb_triangles,0.) ;
		 	idx_col[nb_idx-2]=IDS("CENTREY") ; datas[nb_idx-2].resize(nb_triangles,0.) ;
		 	idx_col[nb_idx-1]=IDS("CENTREZ") ; datas[nb_idx-1].resize(nb_triangles,0.) ;

		 	for (k=0; k<nb_triangles ; k++)
		 		{
		 		centres=get_tri_center(k) ;
		 		datas[nb_idx-3][k]=centres(1) ;	datas[nb_idx-2][k]=centres(2) ;	datas[nb_idx-1][k]=centres(3) ;
		 		}
		 	break ;

default : DISP_Err("Erreur : la reconstruction demandée pour le triangle n'est pas possible\n") ; break ;
}
return 1 ;
}
//----------------------------------------
//int Step::couple(Matrix & torque, Matrix c)
int Step::Fcouple(Vector & torque, Vector c)
{
int i, compt, comptb ; int idx[9] ; static int cptold=-1, cptbold=-1 ;
//Matrix g(3,1),force(3,1) ;
Vector g, force, tmp ;
torque=0 ;
idx[0]=find_idx(IDS("POINTX")) ; idx[1]=find_idx(IDS("POINTY")) ; idx[2]=find_idx(IDS("POINTZ")) ;
idx[3]=find_idx(IDS("FORCEX")) ; idx[4]=find_idx(IDS("FORCEY")) ; idx[5]=find_idx(IDS("FORCEZ")) ;
idx[6]=find_idx(IDS("POLY2")) ; idx[7]=find_idx(IDS("POLY3")) ; idx[8]=find_idx(IDS("POLY4")) ;

if (idx[0]<0 || idx[1]<0 || idx[2]<0) {cout << "ERR : pas de donnees de FORCE\n" ; fflush(stdout); }
compt=0 ; comptb=0 ;
for (i=0 ; i<nb_triangles ; i++)
  {
  // Recherche du centre du triangle
  g=get_tri_center (i) ;
  if(actions["iscylperiodic"].set && (g(1)<0 || g(2)<0)) {compt++ ; continue ; }
  if(actions["w/otige"].set && ((g(3)-c(3)-0.0001)>actions.Cst["R"])) {comptb++ ; continue ; }
  force(1)=datas[idx[3]][i] ; force(2)=datas[idx[4]][i] ; force(3)=datas[idx[5]][i] ;
  tmp=g-c ; tmp=tmp.cross(force) ;
  torque=torque+tmp ; //Calcul::cross_product((g-c),force) ;
  }
if (actions["iscylperiodic"].set)
   {
   if (compt!=cptold) {cptold=compt ; printf("\n[%d faces oubliées sur %d faces]", compt, nb_triangles) ;}
   else printf(".") ;
   }
if (actions["w/otige"].set)
   {
   if (comptb!=cptbold) {cptbold=comptb ; printf("\n[%d faces oubliées sur %d faces en w/o tige]", comptb, nb_triangles) ;}
   }
//torque=torque/nb_triangles ;    STUPID !!!!!!!!!!!!!!!!!
return 1 ;
}

//--------------------------------------------
/*Matrix Step::get_tri_center (int polyidx)
{
int idx[9] ;
Matrix pt1(3,1), pt2(3,1), pt3(3,1), res(3,1) ;
idx[0]=find_idx(IDS("POINTX")) ; idx[1]=find_idx(IDS("POINTY")) ; idx[2]=find_idx(IDS("POINTZ")) ;
//idx[3]=find_idx(IDS("FORCEX")) ; idx[4]=find_idx(IDS("FORCEY")) ; idx[5]=find_idx(IDS("FORCEZ")) ;
idx[6]=find_idx(IDS("POLY2")) ; idx[7]=find_idx(IDS("POLY3")) ; idx[8]=find_idx(IDS("POLY4")) ;
pt1=Convert::toVect(datas[idx[0]][datas[idx[6]][polyidx]],datas[idx[1]][datas[idx[6]][polyidx]], datas[idx[2]][datas[idx[6]][polyidx]]) ;
pt2=Convert::toVect(datas[idx[0]][datas[idx[7]][polyidx]],datas[idx[1]][datas[idx[7]][polyidx]], datas[idx[2]][datas[idx[7]][polyidx]]) ;
pt3=Convert::toVect(datas[idx[0]][datas[idx[8]][polyidx]],datas[idx[1]][datas[idx[8]][polyidx]], datas[idx[2]][datas[idx[8]][polyidx]]) ;
res=Geometrie::gravicentre(pt1, pt2, pt3) ;
return res ;
}*/
Vector Step::get_tri_center (int polyidx)
{
int idx[9] ;
Vector pt1, pt2, pt3 ;
Vector *res = new Vector ;
idx[0]=find_idx(IDS("POINTX")) ; idx[1]=find_idx(IDS("POINTY")) ; idx[2]=find_idx(IDS("POINTZ")) ;
idx[6]=find_idx(IDS("POLY2")) ; idx[7]=find_idx(IDS("POLY3")) ; idx[8]=find_idx(IDS("POLY4")) ;
pt1(1)=datas[idx[0]][datas[idx[6]][polyidx]] ; pt1(2)=datas[idx[1]][datas[idx[6]][polyidx]] ; pt1(3)=datas[idx[2]][datas[idx[6]][polyidx]] ;
pt2(1)=datas[idx[0]][datas[idx[7]][polyidx]] ; pt2(2)=datas[idx[1]][datas[idx[7]][polyidx]] ; pt2(3)=datas[idx[2]][datas[idx[7]][polyidx]] ;
pt3(1)=datas[idx[0]][datas[idx[8]][polyidx]] ; pt3(2)=datas[idx[1]][datas[idx[8]][polyidx]] ; pt3(3)=datas[idx[2]][datas[idx[8]][polyidx]] ;
(*res)=Geometrie::gravicentre(pt1, pt2, pt3) ;
return (*res) ;
}
//-------------------------------------------
Vector Step::get_tri_normal (int polyidx)
{
int idx[9] ;
Vector pt1, pt2, pt3, res ;
idx[0]=find_idx(IDS("POINTX")) ; idx[1]=find_idx(IDS("POINTY")) ; idx[2]=find_idx(IDS("POINTZ")) ;
idx[6]=find_idx(IDS("POLY2")) ; idx[7]=find_idx(IDS("POLY3")) ; idx[8]=find_idx(IDS("POLY4")) ;
pt1(1)=datas[idx[0]][datas[idx[6]][polyidx]] ; pt1(2)=datas[idx[1]][datas[idx[6]][polyidx]] ; pt1(3)=datas[idx[2]][datas[idx[6]][polyidx]] ;
pt2(1)=datas[idx[0]][datas[idx[7]][polyidx]] ; pt2(2)=datas[idx[1]][datas[idx[7]][polyidx]] ; pt2(3)=datas[idx[2]][datas[idx[7]][polyidx]] ;
pt3(1)=datas[idx[0]][datas[idx[8]][polyidx]] ; pt3(2)=datas[idx[1]][datas[idx[8]][polyidx]] ; pt3(3)=datas[idx[2]][datas[idx[8]][polyidx]] ;
pt2=pt2-pt1 ;
pt3=pt3-pt1 ;
res=pt2.cross(pt3);
res=res/Calcul::norm(res) ;
return res ;
}
//-------------------------------------------
double Step::get_tri_surface (int polyidx)
{
int idx[9] ;
Vector pt1, pt2, pt3; double res ;
idx[0]=find_idx(IDS("POINTX")) ; idx[1]=find_idx(IDS("POINTY")) ; idx[2]=find_idx(IDS("POINTZ")) ;
//idx[3]=find_idx(IDS("FORCEX")) ; idx[4]=find_idx(IDS("FORCEY")) ; idx[5]=find_idx(IDS("FORCEZ")) ;
idx[6]=find_idx(IDS("POLY2")) ; idx[7]=find_idx(IDS("POLY3")) ; idx[8]=find_idx(IDS("POLY4")) ;
pt1(1)=datas[idx[0]][datas[idx[6]][polyidx]] ; pt1(2)=datas[idx[1]][datas[idx[6]][polyidx]] ; pt1(3)=datas[idx[2]][datas[idx[6]][polyidx]] ;
pt2(1)=datas[idx[0]][datas[idx[7]][polyidx]] ; pt2(2)=datas[idx[1]][datas[idx[7]][polyidx]] ; pt2(3)=datas[idx[2]][datas[idx[7]][polyidx]] ;
pt3(1)=datas[idx[0]][datas[idx[8]][polyidx]] ; pt3(2)=datas[idx[1]][datas[idx[8]][polyidx]] ; pt3(3)=datas[idx[2]][datas[idx[8]][polyidx]] ;
pt2=pt2-pt1 ;
pt3=pt3-pt1 ;
res=Calcul::norm(pt2.cross(pt3))/2.0 ;
return res ;
}
//--------------------------------------------
int Step::Fatm_rotate (Matrix3x3 & rot, int id)
{
int idx[3] ;
Vector pt1 ;
printf("Should check that Step::Fatm_rotate") ; fflush(stdout) ;
if (id>=nb_pts)
  {
  id-=nb_pts ;
  idx[0]=find_idx(IDS("POINTX")) ; idx[1]=find_idx(IDS("POINTY")) ; idx[2]=find_idx(IDS("POINTZ")) ;
  if (idx[0]!=-1 && idx[1]!=-1 && idx[2]!=-1)
   {
   pt1.set(datas[idx[0]][id], datas[idx[1]][id], datas[idx[2]][id]) ;
   pt1=Geometrie::rotation(pt1, rot) ;
   datas[idx[0]][id]=pt1(1) ; datas[idx[1]][id]=pt1(2) ; datas[idx[2]][id]=pt1(3) ;
   }
  }

else
  {
  idx[0]=find_idx(IDS("FORCEX")) ; idx[1]=find_idx(IDS("FORCEY")) ; idx[2]=find_idx(IDS("FORCEZ")) ;
  if (idx[0]!=-1 && idx[1]!=-1 && idx[2]!=-1)
   {
   pt1.set(datas[idx[0]][id], datas[idx[1]][id], datas[idx[2]][id]) ;
   pt1=Geometrie::rotation(pt1, rot) ;
   datas[idx[0]][id]=pt1(1) ; datas[idx[1]][id]=pt1(2) ; datas[idx[2]][id]=pt1(3) ;
   }

  idx[0]=find_idx(IDS("NORMALEX")) ; idx[1]=find_idx(IDS("NORMALEY")) ; idx[2]=find_idx(IDS("NORMALEZ")) ;
  if (idx[0]!=-1 && idx[1]!=-1 && idx[2]!=-1)
   {
   pt1.set(datas[idx[0]][id], datas[idx[1]][id], datas[idx[2]][id]) ;
   pt1=Geometrie::rotation(pt1, rot) ;
   datas[idx[0]][id]=pt1(1) ; datas[idx[1]][id]=pt1(2) ; datas[idx[2]][id]=pt1(3) ;
   }
  }

return 1 ;
}

