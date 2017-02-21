#include "Headers/Multisphere.h"

int Multisphere::init(Step &step)
{
int i, j ; int idx[5] ; int longest ;
Vector null(0) ;

idx[0]=step.find_idx(IDS("POSX")) ; idx[1]=step.find_idx(IDS("POSY")) ;  idx[2]=step.find_idx(IDS("POSZ")) ;  
idx[3]=step.find_idx(IDS("IDMULTISPHERE")) ; idx[4]=step.find_idx(IDS("ID")) ; 

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

data.resize(7) ; 
for (i=0;i<7 ; i++) data[i].resize(ngp+1,0) ; 

initialized=true ; 
return 0 ;
}

//=====================================================================================================
int Multisphere::get_orientations (Step & step)
{
if (!initialized) init(step) ;   
  
double box[6] ;
int j, k, l, n ; int idx[5] ; 
Vector t ; 

Vector centroid ; int longest ; double maxlen ;
Vector vsph, null(0) ; int idmax ; 
double radius ; 

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
  for (k=0 ; k<gps[j][0] ; k++)
  {
    if (step.datas[idx[4]][gps[j][k+1]-1]!=gps[j][k+1]) 
    {
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
  
  vsph=Geometrie::cart2sph(segments[idmax]) ;
  if (type==0)
  {
    if (vsph(1)>radius*gps[j][0]*2) 
    {data[0][j]=GP_PBC ; segments[idmax](1)=NAN ; segments[idmax](2)=NAN ; segments[idmax](3)=NAN ;}
  }
  else if (type==1)
  {
    if (vsph(1)>radius*(floor(log2((gps[j][0]-1)/3.))*2)) {data[0][j]=GP_PBC ; segments[idmax](1)=NAN ; segments[idmax](2)=NAN ; segments[idmax](3)=NAN ;}
  }        
  data[1][j]=centroid(1) ; 
  data[2][j]=centroid(2) ; 
  data[3][j]=centroid(3) ; 
  data[4][j]=segments[idmax](1) ; 
  data[5][j]=segments[idmax](2) ; 
  data[6][j]=segments[idmax](3) ; 
}

return 0 ; 
}






