#include "Headers/Surface.h"


int Surface::detect_surface (Dump & dump)
{
long int debut, fin, step ; 
double largeur_bloc[2] ;
vector <vector <double> > vv_double_null ;
vector <double> v_double_null ;
int idx[3], bloc[2] ;
int i, j, k ;

if (actions["extract"].set)
   {debut=(int)actions["extract"]["extract_deb"] ; fin=(int)actions["extract"]["extract_fin"] ;
     if (debut<0) {debut=0 ; cout << _("WARNING : le debut de l'extraction est <0. RAZ.\n") ; }
     if (fin>dump.nbsteps) {fin=dump.nbsteps ; cout << _("WARNING : la fin de l'extraction est >nbsteps. Remise à nbsteps.\n") ; }
  }
else
   {debut=0 ; fin=dump.nbsteps;}

if (actions["downsampling"].set)
   { step=(int)actions["downsampling"]["downsampling"] ; 
   if (step<=0 || step>nbsteps) { step=1 ; cout << _("WARNING : downsampling incorrect. Remis à 1.\n") ;}
   }
else
   step=1 ; 

nbsteps=(fin-debut)/step ; 
nbbox[0]=actions["surface"]["nbbox_x"] ;
nbbox[1]=actions["surface"]["nbbox_y"] ;

actions.total=nbsteps ; 
cout << "\nSurface::detect_surface          " ;
actions.disp_progress() ; 

positions.resize(nbsteps, vv_double_null) ;
positionsbas.resize(nbsteps, vv_double_null) ;

for (i=0 ; i<nbsteps ; i++)
{
 positions[i].resize(nbbox[0], v_double_null) ;
 positionsbas[i].resize(nbbox[0], v_double_null) ;
 for (j=0 ;j<nbbox[0] ; j++)
 {
	 positions[i][j].resize(nbbox[1], dump.steps[debut].box[2][0]) ;
	 positionsbas[i][j].resize(nbbox[1], dump.steps[debut].box[2][0]) ;
 }
}
for (i=0 ; i<nbsteps ; i++)
	{for (j=0 ; j<nbbox[0] ; j++)
	    {for (k=0 ; k<nbbox[1] ; k++)
	        {
		    positions[i][j][k]=-1000000 ;
		    positionsbas[i][j][k]=1000000 ;
	        }}}

j=0 ;
for (i=debut ; i<fin ; i+=step)
  {
  actions.valeur=j ; 
  dump.check_timestep(i) ; 
  
  idx[0]=dump.steps[i].find_idx(IDS("POSX")) ; 
  idx[1]=dump.steps[i].find_idx(IDS("POSY")) ;
  idx[2]=dump.steps[i].find_idx(IDS("POSZ")) ;
  largeur_bloc[0]=(dump.steps[i].box[0][1]-dump.steps[i].box[0][0])/nbbox[0] ;
  largeur_bloc[1]=(dump.steps[i].box[1][1]-dump.steps[i].box[1][0])/nbbox[1] ;
  for (k=0 ; k<dump.steps[i].nb_atomes ; k++)
     {
     bloc[0]=floor((dump.steps[i].datas[idx[0]][k]-dump.steps[i].box[0][0])/largeur_bloc[0]) ;
     bloc[0]=(bloc[0]<0?0:bloc[0]) ;
     bloc[0]=(bloc[0]>=nbbox[0]?(nbbox[0]-1):bloc[0]) ;
     
     bloc[1]=floor((dump.steps[i].datas[idx[1]][k]-dump.steps[i].box[1][0])/largeur_bloc[1]) ;
     bloc[1]=(bloc[1]<0?0:bloc[1]) ;
     bloc[1]=(bloc[1]>=nbbox[1]?(nbbox[1]-1):bloc[1]) ;

     if (positions[j][bloc[0]][bloc[1]]<dump.steps[i].datas[idx[2]][k])
        positions[j][bloc[0]][bloc[1]]=dump.steps[i].datas[idx[2]][k] ;
     if (positionsbas[j][bloc[0]][bloc[1]]>dump.steps[i].datas[idx[2]][k])
             positionsbas[j][bloc[0]][bloc[1]]=dump.steps[i].datas[idx[2]][k] ;
      }
  j++ ;
  } 
return 1 ;
}

int Surface::detect_surface2D (Dump & dump)
{
long int debut, fin, step ;
double largeur_bloc ;
vector <double> v_double_null ;
int idx[3], bloc ;
int i, j, k ;

if (actions["extract"].set)
   {debut=(int)actions["extract"]["extract_deb"] ; fin=(int)actions["extract"]["extract_fin"] ;
     if (debut<0) {debut=0 ; cout << _("WARNING : le debut de l'extraction est <0. RAZ.\n") ; }
     if (fin>dump.nbsteps) {fin=dump.nbsteps ; cout << _("WARNING : la fin de l'extraction est >nbsteps. Remise à nbsteps.\n") ; }
  }
else
   {debut=0 ; fin=dump.nbsteps;}

if (actions["downsampling"].set)
   { step=(int)actions["downsampling"]["downsampling"] ;
   if (step<=0 || step>nbsteps) { step=1 ; cout << _("WARNING : downsampling incorrect. Remis à 1.\n") ;}
   }
else
   step=1 ;

nbsteps=(fin-debut)/step ;
nbbox[0]=actions["surface2D"]["nbbox_x"] ;

actions.total=nbsteps ;
cout << "\nSurface::detect_surface2D          " ;
actions.disp_progress() ;

positions2D.resize(nbsteps, v_double_null) ;

for (i=0 ; i<nbsteps ; i++)
{positions2D[i].resize(nbbox[0], dump.steps[debut].box[2][0]) ;}
for (i=0 ; i<nbsteps ; i++)
	{for (j=0 ; j<nbbox[0] ; j++)
	    {
		    positions2D[i][j]=-1000000 ;
	        }}

j=0 ;
for (i=debut ; i<fin ; i+=step)
  {
  actions.valeur=j ;
  dump.check_timestep(i) ;

  idx[0]=dump.steps[i].find_idx(IDS("POSX")) ;
  idx[1]=dump.steps[i].find_idx(IDS("POSY")) ;
  idx[2]=dump.steps[i].find_idx(IDS("POSZ")) ;
  largeur_bloc=(dump.steps[i].box[0][1]-dump.steps[i].box[0][0])/nbbox[0] ;
  for (k=0 ; k<dump.steps[i].nb_atomes ; k++)
     {
     bloc=floor((dump.steps[i].datas[idx[0]][k]-dump.steps[i].box[0][0])/largeur_bloc) ;
     bloc=(bloc<0?0:bloc) ;
     bloc=(bloc>=nbbox[0]?(nbbox[0]-1):bloc) ;

     if (positions2D[j][bloc]<dump.steps[i].datas[idx[1]][k])
        positions2D[j][bloc]=dump.steps[i].datas[idx[1]][k] ;
      }
  j++ ;
  }
return 1 ;
}

//=====================================================
int Surface::write_asmatlab (string chemin)
{
ofstream out ; 
stringstream chemin2 ;
int i, j, k ;

chemin2.str("") ;
try{chemin.erase(chemin.rfind(".gz")) ;}
catch (...) {}
chemin2 << chemin << "-Surface" << ".txt" ;  
out.open(chemin2.str().c_str(), ios::out) ; 

// On écrit la commande matlab de redimensionnement, histoire de l'avoir en dur. On l'écrit aussi à l'écran.
out << _("% Commande matlab de redimensionnement : reshape(X,")<<((int)(nbbox[1]))<<","<<((int)(nbbox[0]))<<","<<((int)(nbsteps))<<");\n" ;
cout<< _("Commande matlab de redimensionnement : reshape(X,")<<((int)(nbbox[1]))<<","<<((int)(nbbox[0]))<<","<<((int)(nbsteps))<<");\n" ;

for (i=0 ; i<nbsteps ; i++)
  {
  for (j=0 ; j<nbbox[0] ; j++)
    {
	for (k=0 ; k<nbbox[1] ; k++)
		{out << positions[i][j][k] << " " ; }
    }
  }
  
out.close() ;   
return 1 ;
}
//----------------------------------------------------
int Surface::write_asmatlab2D (string chemin)
{
ofstream out ;
stringstream chemin2 ;
int i, j, k ;

chemin2.str("") ;
try{chemin.erase(chemin.rfind(".gz")) ; }
catch (...) {}
chemin2 << chemin << "-Surface" << ".txt" ;
out.open(chemin2.str().c_str(), ios::out) ;

// On écrit la commande matlab de redimensionnement, histoire de l'avoir en dur. On l'écrit aussi à l'écran.
for (i=0 ; i<nbsteps ; i++)
  {
  for (j=0 ; j<nbbox[0] ; j++)
    { out << positions2D[i][j] << " " ; }
  out << "\n" ;
  }

out.close() ;
return 1 ;
}
//====================================================
int Surface::write_asmatlab_bottom (string chemin)
{
ofstream out ;
stringstream chemin2 ;
int i, j, k ;

chemin2.str("") ;
chemin.erase(chemin.rfind(".gz")) ; 
chemin2 << chemin << "-SurfaceBas" << ".txt" ;
out.open(chemin2.str().c_str(), ios::out) ;

// On écrit la commande matlab de redimensionnement, histoire de l'avoir en dur. On l'écrit aussi à l'écran.
out << _("% Commande matlab de redimensionnement : reshape(X,")<<((int)(nbbox[1]))<<","<<((int)(nbbox[0]))<<","<<((int)(nbsteps))<<");\n" ;
cout<< _("Commande matlab de redimensionnement : reshape(X,")<<((int)(nbbox[1]))<<","<<((int)(nbbox[0]))<<","<<((int)(nbsteps))<<");\n" ;

for (i=0 ; i<nbsteps ; i++)
  {
  for (j=0 ; j<nbbox[0] ; j++)
    {
	for (k=0 ; k<nbbox[1] ; k++)
		{out << positionsbas[i][j][k] << " " ; }
    }
  }

out.close() ;
return 1 ;
}



