#include "Headers/Multisphere.h"

Multisphere::Multisphere() :initialized(false),ngp(0), currentstepinit(false), currentstep(-1), isflux(false) 
  { if (actions["multisphereflux"].set) isflux=true ;  for (int i=0; i<7; i++) symetrie[i]=false ; } 

int Multisphere::init(Step &step)
{
int i, j ; int idx[5] ; int longest ;
Vector null(0) ;

idx[0]=step.find_idx(IDS("POSX")) ; idx[1]=step.find_idx(IDS("POSY")) ;  idx[2]=step.find_idx(IDS("POSZ")) ;  
idx[3]=step.find_idx(IDS("IDMULTISPHERE")) ; idx[4]=step.find_idx(IDS("ID")) ; 

if (actions["symetriser"].set)
  {
    symetrie[0]=true ; double r=actions["symetriser"]["axes"] ;
    if (r>=100) {symetrie[1]=true ; r=r-100 ; }
    if (r>=10)  {symetrie[2]=true ; r=r-10 ; }
    if (r>=1)   {symetrie[3]=true ;}
  }
printf("%d %d %d %d %d %d %d-----------------\n", symetrie[0], symetrie[1],symetrie[2],symetrie[3],symetrie[4],symetrie[5],symetrie[6]) ; 

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
if (ngp==-1) {DISP_Warn("Aucun groupe multisphere trouvé, il y a un problème.") ; printf("Step ID=%d", currentstep) ; fflush(stdout) ;   }   

for (i=1, longest=0 ; i<=ngp ; i++) if (longest<gps[i][0]) longest=gps[i][0] ; 
pts.resize(longest, null) ; 
segments.resize(longest*(longest-1)/2, null) ;  

data.resize(7) ; 
for (i=0;i<7 ; i++) data[i].resize(ngp+1,0) ; 

initialized=true ; 
return 0 ;
}

//=====================================================================================================
int Multisphere::get_orientations_noflux (Step & step)
{
double box[6] ;
int j, k, l, n ; int idx[5] ; 
Vector t ; 

Vector centroid ; int longest ; double maxlen ;
Vector vsph, null(0) ; int idmax ; 
double radius ; 

if (!initialized) init(step) ;   

radius=actions.Cst["Radius"] ; 
idx[0]=step.find_idx(IDS("POSX")) ; idx[1]=step.find_idx(IDS("POSY")) ;  idx[2]=step.find_idx(IDS("POSZ")) ;  
idx[3]=step.find_idx(IDS("IDMULTISPHERE")) ; idx[4]=step.find_idx(IDS("ID")) ; 

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

for (j=1 ; j<=ngp ; j++)
{
  centroid=0 ;   
  if (data[0][j]==GP_LOST) continue ;  
  if (data[0][j]<GP_LOST) data[0][j]=GP_OK ; // A priori les groupes qui ne sont pas perdus ou problematique sont OK. Eventuellement ils seront marqués PBC ou OUT ensuite. 
  for (k=0 ; k<gps[j][0] ; k++)
  {
    if (step.datas[idx[4]][gps[j][k+1]-1]!=gps[j][k+1]) 
    {
      printf("%g %g ", step.datas[idx[4]][gps[j][k+1]-1],gps[j][k+1] ) ; 
      DISP_Err("Probleme in multisphere ici\n") ; data[0][j]=GP_BAD ; 
    } 
    t.set(step.datas[idx[0]][gps[j][k+1]-1], step.datas[idx[1]][gps[j][k+1]-1], step.datas[idx[2]][gps[j][k+1]-1]);   
    pts[k]=t ;
    if (t.isnan()) 
    {
      data[0][j]=GP_LOST ; printf("Le groupe %d a été perdu. Atomes:", j) ; for (l=0 ; l<gps[j][0] ; l++) {printf("%d ", gps[j][l+1]-1) ; } printf("\n") ; break ;  
    }
    centroid=centroid+pts[k] ;
  } 
  if (data[0][j]==GP_LOST) continue ; 
  
  centroid=centroid/gps[j][0] ; 
  if (symetrie[0]) 
     {
       if (symetrie[1]==true && centroid(1)<0) {symetrie[4]=true ; centroid(1)=-centroid(1) ; } else {symetrie[4]=false ; }
       if (symetrie[2]==true && centroid(2)<0) {symetrie[5]=true ; centroid(2)=-centroid(2) ; } else {symetrie[5]=false ; }
       if (symetrie[3]==true && centroid(3)<0) {symetrie[6]=true ; centroid(3)=-centroid(3) ; } else {symetrie[6]=false ; }
     }
  
  if (centroid(1)<box[0] || centroid(1)>box[1] || centroid(2)<box[2] || centroid[2]>box[3] || centroid(3)<box[4] || centroid(3)>box[5]) { data[0][j]=GP_OUT ; continue ; }
  
  
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

 if (symetrie[0])
     {
       if (symetrie[4]) {segments[idmax](1)=-segments[idmax](1) ; }
       if (symetrie[5]) {segments[idmax](2)=-segments[idmax](2) ; }
       if (symetrie[6]) {segments[idmax](3)=-segments[idmax](3) ; }
     }
  
  vsph=Geometrie::cart2sph(segments[idmax]) ;
  if (type==0)
  {
    if (vsph(1)>radius*gps[j][0]*2) 
    {data[0][j]=GP_PBC ; segments[idmax](1)=NAN ; segments[idmax](2)=NAN ; segments[idmax](3)=NAN ;}
  }
  else if (type==1)
  {
    //printf("%g %g \n", vsph(1), radius*(floor(log2((gps[j][0]-1)/3.))*2)) ; 
    //if (vsph(1)>radius*(floor(log2((gps[j][0]-1)/3.))*2)) {
      if (vsph(1)>radius*2.1) {data[0][j]=GP_PBC ; printf("s") ; segments[idmax](1)=NAN ; segments[idmax](2)=NAN ; segments[idmax](3)=NAN ;}
  }        
  data[1][j]=centroid(1) ; 
  data[2][j]=centroid(2) ; 
  data[3][j]=centroid(3) ; 
  data[4][j]=segments[idmax](1) ; 
  data[5][j]=segments[idmax](2) ; 
  data[6][j]=segments[idmax](3) ; 
  if ((isnan(data[4][j]) || isnan(data[5][j]) || isnan(data[6][j])) && data[0][j]==GP_OK) 
  {DISP_Warn("NaN dans le data multisphere avec GP_OK, probleme\n") ; data[0][j]=GP_BAD ;
  }
}

currentstepinit=true ; 

return 0 ; 
}

//================================ get_orientations for flux of multispheres. Nop filter::fill should be applied, because it would be too slow for large fluxes ==============
//Reasonable assumption: ids for spheres and multispheres are never reattibuted if the sphere/multisphere disappears -------
int Multisphere::init_flux (Step &step ) // Should be called each timestep
{
int i, j, k ; int idx[5] ; int longest ;
Vector null(0) ; static bool firstinit = true ; 
bool isnew ; 
  
idx[3]=step.find_idx(IDS("IDMULTISPHERE")) ; idx[4]=step.find_idx(IDS("ID")) ; 

if (actions["symetriser"].set)
  {
    DISP_Warn("symetriser not implemented in multisphere for fluxes\n") ; 
  }

type=actions["multisphereflux"]["type"] ; 

j=0 ; isnew=false ; k=1 ; 
for (i=0 ; i<step.nb_atomes ; i++)
{
  if (step.datas[idx[3]][i]<=0) continue ; 
  
  if (firstinit && i==0) {gps_id.push_back(step.datas[idx[3]][i]) ; gps.resize(gps.size()+1) ; gps[j].push_back(0) ; isnew=true ; }
  
  if (step.datas[idx[3]][i]==gps_id[j])
  {   
    if (isnew) {gps[j][0]++ ; gps[j].push_back(i) ; }
    else {gps[j][k]=i ; k++ ; }
  }
  else // step.datas[idx[3]][i]>gps_id[j]
  {
   j++ ; 
   for (k=j ; k<gps_id.size() ; k++)
     if (step.datas[idx[3]][i]==gps_id[j])
       break ; 
   if (k==gps_id.size()) isnew=true ; 
   gps_id.erase(gps_id.begin()+j, gps_id.begin()+k) ; 
   gps.erase(gps.begin()+j, gps.begin()+k) ; 
   if (isnew)
   {
     gps_id.push_back(step.datas[idx[3]][i]) ;
     gps.resize(gps.size()+1) ; gps[j].push_back(0) ; 
     gps[j][0]++; gps[j].push_back(i) ; 
   }
   else
   {gps[j][1]=i ; k=2 ; }
  }
}
ngp=gps_id.size() ; 
//printf("[%d]\n", ngp) ;

if (firstinit==true) //WARNING will not work if the max size of particle evolves.
{
  for (i=0, longest=0 ; i<ngp ; i++) if (longest<gps[i][0]) longest=gps[i][0] ; 
  pts.resize(longest, null) ; 
  segments.resize(longest*(longest-1)/2, null) ;  
  data.resize(7) ;
  firstinit=false ; 
}
if (data[0].size() < ngp+1) //+1 to be safe, not needed actually
  for (i=0;i<7 ; i++) data[i].resize(ngp+50,0) ; // +50 to avoid too many reallocations

return 0 ;  
}
//---------------------------------------------------------
int Multisphere::get_orientations_flux (Step & step)
{
double box[6] ;
int j, k, l, n ; int idx[5] ; 
Vector t ; 

Vector centroid ; int longest ; double maxlen ;
Vector vsph, null(0) ; int idmax ; 
double radius ; 

init_flux(step) ; // Always initialize 

radius=actions.Cst["Radius"] ; 
idx[0]=step.find_idx(IDS("POSX")) ; idx[1]=step.find_idx(IDS("POSY")) ;  idx[2]=step.find_idx(IDS("POSZ")) ;  
idx[3]=step.find_idx(IDS("IDMULTISPHERE")) ; idx[4]=step.find_idx(IDS("ID")) ; 

if (actions["use-box"].set) 
{box[0] =  actions["use-box"]["box_xmin"] ; box[1] =  actions["use-box"]["box_xmax"] ;
 box[2] =  actions["use-box"]["box_ymin"] ; box[3] =  actions["use-box"]["box_ymax"] ;
 box[4] =  actions["use-box"]["box_zmin"] ; box[5] =  actions["use-box"]["box_zmax"] ;}
else
{box[0]=box[2]=box[4]=-std::numeric_limits <double>::infinity() ;                      
 box[1]=box[3]=box[5]= std::numeric_limits <double>::infinity()  ; }

for (j=0 ; j<ngp ; j++)
{
  centroid=0 ;   
  if (data[0][j]<GP_LOST) data[0][j]=GP_OK ; // A priori les groupes qui ne sont pas perdus ou problematique sont OK. Eventuellement ils seront marqués PBC ou OUT ensuite. 
  
//printf("[%d %d %d]",gps[j][0],gps[j][1],gps[j][2]) ; fflush(stdout) ; 
  for (k=0 ; k<gps[j][0] ; k++)
  {
    t.set(step.datas[idx[0]][gps[j][k+1]], step.datas[idx[1]][gps[j][k+1]], step.datas[idx[2]][gps[j][k+1]]);   
    pts[k]=t ;
    if (t.isnan()) 
    {
      data[0][j]=GP_LOST ; printf("Le groupe %d a été perdu, ne devrait pas arriver pour multisphere flux. Atomes:", j) ; for (l=0 ; l<gps[j][0] ; l++) {printf("%d ", gps[j][l+1]-1) ; } printf("\n") ; break ;  
    }
    centroid=centroid+pts[k] ;
  } 
  if (data[0][j]==GP_LOST) continue ; 
  
  centroid=centroid/gps[j][0] ; 

  if (centroid(1)<box[0] || centroid(1)>box[1] || centroid(2)<box[2] || centroid[2]>box[3] || centroid(3)<box[4] || centroid(3)>box[5]) { data[0][j]=GP_OUT ; fflush(stdout) ; continue ; }
   
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
  
  vsph=Geometrie::cart2sph(segments[idmax]) ;
  if (type==0)
  {
    if (vsph(1)>radius*gps[j][0]*2) 
    {data[0][j]=GP_PBC ; segments[idmax](1)=NAN ; segments[idmax](2)=NAN ; segments[idmax](3)=NAN ;}
  }
  else if (type==1)
  {
    //printf("%g %g \n", vsph(1), radius*(floor(log2((gps[j][0]-1)/3.))*2)) ; 
    //if (vsph(1)>radius*(floor(log2((gps[j][0]-1)/3.))*2)) {
      if (vsph(1)>radius*2.1) {data[0][j]=GP_PBC ; printf("s") ; segments[idmax](1)=NAN ; segments[idmax](2)=NAN ; segments[idmax](3)=NAN ;}
  }        
  data[1][j]=centroid(1) ; 
  data[2][j]=centroid(2) ; 
  data[3][j]=centroid(3) ; 
  data[4][j]=segments[idmax](1) ; 
  data[5][j]=segments[idmax](2) ; 
  data[6][j]=segments[idmax](3) ; 
  if ((isnan(data[4][j]) || isnan(data[5][j]) || isnan(data[6][j])) && data[0][j]==GP_OK) 
  {DISP_Warn("NaN dans le data multisphere avec GP_OK, probleme\n") ; data[0][j]=GP_BAD ;
  }
}
currentstepinit=true ; 

return 0 ; 
}




// -----------------------------------------
Matrix3d Multisphere::compute_K (Step &step)
{
  int Kn ; 
  Matrix3d K, Kvec, Kmatseg ;
  //Map<Vector3d> Ksegment(NULL);
  Vector3d Ksegment ; 
  
  K=Matrix3d::Zero() ; Kn=0 ;
  
  if (! currentstepinit) get_orientations(step) ;
  
  for (int j=1 ; j<=ngp ; j++)
  {
    if (data[0][j]!=GP_OK) continue ; 
     Ksegment[0]=data[4][j] ; Ksegment[1]=data[5][j] ; Ksegment[2]=data[6][j] ; 
     //new (&Ksegment) Map<Vector3d>(&(data[4][j])); // THIS IS NOT AN ALLOCATION (no delete) ; 
     Ksegment=Ksegment/(Ksegment.norm()) ; 
     Kmatseg=Ksegment*(Ksegment.transpose());
     K=K+Kmatseg ; Kn++ ; 
  }
  if (Kn==0) DISP_Err("Kn = 0, problem.") ; 
  K=K/Kn ; 
  return K ;   
}
//--------------------------------------------------
double Multisphere::compute_dzeta (Step & step)
{
  Matrix3d K ;
  K=compute_K(step) ; 
  return (sqrt((3*(K.norm())*(K.norm())-1)/2)) ;
}
//--------------------------------------------------
void Multisphere::compute_eigen(Step &step)
{
 DISP_Err("This function hasn't been tested and probably don't work. Please check the source (Multisphere::compute_eigen()) to check what to do") ; 
 return ; 
 /*    //Calcul::eigen(K, Kval, Kvec) ;
    
    if (type==0)
    {if (Kvec(0,0)<0) {Kvec.col(0)=-Kvec.col(0) ;}}
    else if (type==1)
    {if (Kvec(2,0)<0) {Kvec.col(0)=-Kvec.col(0) ;}}
    else
      DISP_Err("Unknown multisphere type\n") ; 
    
    if (Kvec.determinant()<0) {Kvec.col(2)=-Kvec.col(2) ; }

    //fprintf(out, "%d %g %g %g %g %g %g %g %g %g %g %g %g %g\n", steps[i].timestep, Phi, Kval(0), Kval(1), Kval(2), Kvec(0,0), Kvec(0,1), Kvec(0,2), Kvec(1,0), Kvec(1,1), Kvec(1,2), Kvec(2,0), Kvec(2,1), Kvec(2,2)) ; 
    */
  
}

void Multisphere::check()
{
 for (int j=1; j<ngp ; j++)
 {
   if (isnan(data[4][j])) printf("A%d ", j) ; 
   if (isnan(data[5][j])) printf("B%d ", j) ; 
   if (isnan(data[6][j])) printf("C%d ", j) ; 
 }
}

//-----------------------------------------------
int Multisphere::prepare_Writing (Step & step)  
{

 if (!currentstepinit || isflux) get_orientations(step) ;

 if (isflux)
 {
   /*for (int i=0 ; i<ngp ; i++)
    if (data[0][i] != GP_OK)
    { 
      gps.erase (gps.begin()+i) ; ngp-- ; i-- ; 
      initialized=false ; 
    }*/
 }
 else
 {
  for (int i=1 ; i<=ngp ; i++)
    if (data[0][i] != GP_OK)
    { 
      gps.erase (gps.begin()+i) ; ngp-- ; i-- ; 
      initialized=false ; 
    }
 }
  return 0 ; 
}

//--------------------------------------------
int Multisphere::remove_atoms (Step &step)
{
 vector <int> todelete ;  
 if (!currentstepinit) get_orientations(step) ;

 for (int i=1 ; i<=ngp ; i++)
  if (data[0][i] != GP_OK)
  { 
    for (int j=1 ; j<gps[i][0]+1 ; j++)
      todelete.push_back(gps[i][j]-1) ; 
  }
  std::sort(todelete.begin(), todelete.end()) ; 
 
  for (int i=todelete.size()-1 ; i>=0 ; i--)
  {
    printf("%d ", todelete[i]) ; fflush(stdout) ; 
    step.del_atm(todelete[i]) ; 
  }
  return 0 ; 
}





























