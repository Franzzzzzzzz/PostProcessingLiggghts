#include "Headers/FDump.h"

// =====================================================
// Fonctions de la classe FDump ========================
//=====================================================
int FDump::read (unsigned char read_type, int index)
{
/* read_type==1 : sauver les données lues
   * read_type==2 : ne lire qu'un dump
   * read_type==4 : creer la structure step et incrémenter nb_step */

// lecture d'un dump vtk de forces
string ligne, rien ; int tmp ;
istringstream iligne ;
static bool info=false ; 

Step null_step ;
int retour ; double dtmp ;
int i=0, j, k, stop=0, idx ;

Type=TF ; 

if (read_type & 4) {nbsteps=0 ; printf("Lecture initiale du vtk, peut prendre du temps ...\n") ; }
if (read_type & 2) {idx=index ;}

if (info==false)
DISP_Info("Les tests de réussite de lecture ne sont pas effectués ici. Si des problèmes surviennent, consulter les fonctions LDump::read() pour voir les tests qu'il faudrait effectuer.\n") ; 
info=true ; 

while (!dumpin.eof() && stop<2)
  {
  dumpin.getline(ligne) ;
  retour=identify_ligne_vtk(ligne) ;
  iligne.str(ligne);

  switch(retour)
    {
    case 1 :  if (read_type & 4)
               { nbsteps++ ;
               steps.resize(nbsteps, null_step) ;
               steps[nbsteps-1].posinfile=dumpin.tellg()-(streampos)ligne.length()-(streampos)1 ;
               steps[nbsteps-1].nb_triangles=0 ; steps[nbsteps-1].Type=TF ;
               steps[nbsteps-1].nb_idx=0 ;
               if (index==-1) idx=nbsteps-1 ;
               actions.valeur=dumpin.tellg() ;
               }
              if (read_type & 2) {stop++ ; if (stop>=2) break ;  }
              break ;
    case 2  : if (read_type & 1) {iligne >> rien >> tmp ; steps[idx].nb_triangles=tmp ;}
              break ;
    case 31 : if (read_type & 1)
               {
               dumpin.getline(ligne) ; // ligne de LOOKUP_TABLE
	       steps[idx].nb_idx++ ;
	       steps[idx].idx_col.resize(steps[idx].nb_idx, 0) ;
	       steps[idx].idx_col[steps[idx].nb_idx-1]=IDS("PRESSURE") ;
	       steps[idx].datas.resize(steps[idx].nb_idx, vector<double>(1,0)) ;
	       steps[idx].datas[steps[idx].nb_idx-1].resize(steps[idx].nb_triangles,0.) ;
               for (j=0 ; j<steps[idx].nb_triangles ; j++)
                   {dumpin >> steps[idx].datas[steps[idx].nb_idx-1][j] ; }
               }
               break ;
    case 32 : if (read_type & 1)
               {
               dumpin.getline( ligne) ; // ligne de LOOKUP_TABLE
	       steps[idx].nb_idx++ ;
	       steps[idx].idx_col.resize(steps[idx].nb_idx, 0) ;
	       steps[idx].idx_col[steps[idx].nb_idx-1]=IDS("SHEARSTRESS") ;
	       steps[idx].datas.resize(steps[idx].nb_idx, vector<double>(1,0)) ;
	       steps[idx].datas[steps[idx].nb_idx-1].resize(steps[idx].nb_triangles,0.) ;
               for (j=0 ; j<steps[idx].nb_triangles ; j++)
                   {dumpin >> steps[idx].datas[steps[idx].nb_idx-1][j] ; }
               }
               break ;
    case 41 : if (read_type & 1)
               {
	       steps[idx].nb_idx+=3 ;
	       steps[idx].datas.resize(steps[idx].nb_idx, vector<double>(1,0)) ;
	       steps[idx].idx_col.resize(steps[idx].nb_idx, 0) ;
	       steps[idx].idx_col[steps[idx].nb_idx-3]=IDS("FORCEX") ; steps[idx].datas[steps[idx].nb_idx-3].resize(steps[idx].nb_triangles,0.) ;
	       steps[idx].idx_col[steps[idx].nb_idx-2]=IDS("FORCEY") ; steps[idx].datas[steps[idx].nb_idx-2].resize(steps[idx].nb_triangles,0.) ;
	       steps[idx].idx_col[steps[idx].nb_idx-1]=IDS("FORCEZ") ; steps[idx].datas[steps[idx].nb_idx-1].resize(steps[idx].nb_triangles,0.) ;
               for (j=0 ; j<steps[idx].nb_triangles ; j++)
                 {
                 for (k=0 ; k<3 ; k++)
                  {
                  dumpin >> steps[idx].datas[steps[idx].nb_idx-3+k][j] ;
                  }
                 }
               }
              break ;
    case 42 : if (read_type & 1)
               {
	       steps[idx].nb_idx+=3 ;
	       steps[idx].datas.resize(steps[idx].nb_idx, vector<double>(1,0)) ;
	       steps[idx].idx_col.resize(steps[idx].nb_idx,0) ;
	       steps[idx].idx_col[steps[idx].nb_idx-3]=IDS("NORMALEX") ; steps[idx].datas[steps[idx].nb_idx-3].resize(steps[idx].nb_triangles,0.) ;
	       steps[idx].idx_col[steps[idx].nb_idx-2]=IDS("NORMALEY") ; steps[idx].datas[steps[idx].nb_idx-2].resize(steps[idx].nb_triangles,0.) ;
	       steps[idx].idx_col[steps[idx].nb_idx-1]=IDS("NORMALEZ") ;  steps[idx].datas[steps[idx].nb_idx-1].resize(steps[idx].nb_triangles,0.) ;
               for (j=0 ; j<steps[idx].nb_triangles ; j++)
                 {
                 for (k=0 ; k<3 ; k++)
                  {
                  dumpin >> steps[idx].datas[steps[idx].nb_idx-3+k][j] ;
                  }
                 }
               }
               break ;
    case 43 : if (read_type & 1)
               {
	       steps[idx].nb_idx+=3 ;
	       steps[idx].datas.resize(steps[idx].nb_idx, vector<double>(1,0)) ;
	       steps[idx].idx_col.resize(steps[idx].nb_idx,0) ;
	       steps[idx].idx_col[steps[idx].nb_idx-3]=IDS("STRESSX") ; steps[idx].datas[steps[idx].nb_idx-3].resize(steps[idx].nb_triangles,0.) ;
	       steps[idx].idx_col[steps[idx].nb_idx-2]=IDS("STRESSY") ; steps[idx].datas[steps[idx].nb_idx-2].resize(steps[idx].nb_triangles,0.) ;
	       steps[idx].idx_col[steps[idx].nb_idx-1]=IDS("STRESSZ") ;  steps[idx].datas[steps[idx].nb_idx-1].resize(steps[idx].nb_triangles,0.) ;
               for (j=0 ; j<steps[idx].nb_triangles ; j++)
                 {
                 for (k=0 ; k<3 ; k++)
                  {
                  dumpin >> steps[idx].datas[steps[idx].nb_idx-3+k][j] ;
                  }
                 }
               }
               break ;
    case 5 :   if (read_type & 1)
	        {
                iligne >> rien >> tmp ; steps[idx].nb_pts=tmp ;
		steps[idx].nb_idx+=3 ;
	        steps[idx].datas.resize(steps[idx].nb_idx, vector<double>(1,0)) ;
		steps[idx].idx_col.resize(steps[idx].nb_idx,0) ;
	        steps[idx].idx_col[steps[idx].nb_idx-3]=IDS("POINTX") ;  steps[idx].datas[steps[idx].nb_idx-3].resize(steps[idx].nb_pts,0.) ;
	        steps[idx].idx_col[steps[idx].nb_idx-2]=IDS("POINTY") ;  steps[idx].datas[steps[idx].nb_idx-2].resize(steps[idx].nb_pts,0.) ;
	        steps[idx].idx_col[steps[idx].nb_idx-1]=IDS("POINTZ") ;  steps[idx].datas[steps[idx].nb_idx-1].resize(steps[idx].nb_pts,0.) ;
                for (j=0 ; j<steps[idx].nb_pts ; j++)
                  {
                  for (k=0 ; k<3 ; k++)
                   {
                   dumpin >> steps[idx].datas[steps[idx].nb_idx-3+k][j] ;
                   }
                  }
                }
               break ;
    case 6 :   if (read_type & 1)
                {
                iligne >> rien >> tmp ; steps[idx].nb_triangles=tmp ;
		steps[idx].nb_idx+=4 ;
	        steps[idx].datas.resize(steps[idx].nb_idx, vector<double>(1,0)) ;
	        steps[idx].idx_col.resize(steps[idx].nb_idx,0) ;
		steps[idx].idx_col[steps[idx].nb_idx-4]=IDS("POLY1") ;  steps[idx].datas[steps[idx].nb_idx-4].resize(steps[idx].nb_triangles,0.) ;
	        steps[idx].idx_col[steps[idx].nb_idx-3]=IDS("POLY2") ;  steps[idx].datas[steps[idx].nb_idx-3].resize(steps[idx].nb_triangles,0.) ;
	        steps[idx].idx_col[steps[idx].nb_idx-2]=IDS("POLY3") ;  steps[idx].datas[steps[idx].nb_idx-2].resize(steps[idx].nb_triangles,0.) ;
	        steps[idx].idx_col[steps[idx].nb_idx-1]=IDS("POLY4") ;  steps[idx].datas[steps[idx].nb_idx-1].resize(steps[idx].nb_triangles,0.) ;
                for (j=0 ; j<steps[idx].nb_triangles ; j++)
                  {
                  for (k=0 ; k<4 ; k++)
                   {
                   dumpin >> steps[idx].datas[steps[idx].nb_idx-4+k][j] ;
                   }
                  }
                }
               break ;
    }
  }
dumpin.clear() ;
return 0 ;
}
//----------------------------------------------------
int FDump::identify_ligne_vtk (string ligne)
{
size_t position ; static bool warn=true ; 

position=ligne.find('#') ;
if (position==0)
   return 1 ;
position=ligne.find("CELL_DATA") ;
if (position!=string::npos)
   return 2 ;
position=ligne.find("SCALARS") ;
if (position!=string::npos)
   {
   position=ligne.find("pressure") ;
   if (position!=string::npos)
      return 31 ;
   position=ligne.find("shearstress") ;
   if (position!=string::npos)
      return 32 ;
   DISP_Warn("WARN8 : type de scalaire inconnu dans le vtk\n") ; return 30 ;
   }
position=ligne.find("VECTORS") ;
if (position!=string::npos)
   {
   position=ligne.find("forcesTri") ;
   if (position!=string::npos)
      {if (!actions["v1"].set && warn ){DISP_Warn("Les donnees forcesTri ne devraients pas exister en Lv2. Utiliser --v1 peut aider\n") ; warn=false ;} return 41 ;}
   position=ligne.find("normales") ;
   if (position!=string::npos)
	  {if (!actions["v1"].set && warn){DISP_Warn("Les donnees normales ne devraients pas exister en Lv2. Utiliser --v1 peut aider\n") ; warn=false ; } return 42 ; }
   position=ligne.find("stress") ;
      if (position!=string::npos)
   	  {if (actions["v1"].set && warn){DISP_Warn("Les donnees stress ne devraients pas exister en Lv1.\n") ; warn=false ; } return 43 ; }
   DISP_Warn("WARN9 : type de vecteur inconnu dans le vtk\n") ; return 40 ;
   }
position=ligne.find("POINTS") ;
if (position!=string::npos)
   return 5 ;
// POLYGONS en Lv1 est la même chose (y compris la syntaxe) en Lv2. La lecture (et donc la valeur de renvoie) est donc la même
position=ligne.find("POLYGONS") ;
if (position!=string::npos)
	{if (!actions["v1"].set && warn){DISP_Warn("Les donnees forcesTri ne devraients pas exister en Lv2. Utiliser --v1 peut aider\n") ; warn=false ; } return 6 ;}
position=ligne.find("CELLS") ;
if (position!=string::npos)
	{if (actions["v1"].set && warn){DISP_Warn("Les donnees CELLS devraient être POLYGONS pour la Lv1.\n") ; warn=false ;} return 6 ;}

return 0 ;
}
//---------------------------------------------
void FDump::disp(void)
{
 Dump::disp() ;
 check_timestep(0) ; steps[0].disp() ;
 cout << "------\n" ;
 check_timestep(nbsteps-1) ; steps[nbsteps-1].disp() ;
 cout <<"--------------------\n\n" ;
}
//-----------------------------------------
int FDump::write_forcestot (string chemin)
{
string chem, chem2 ; int i;
int nb_tri=-1 ;
ofstream out, groupout ;
double forces[3] ; long int loop[3] ;

// Variables utilisés en cas de group triangles
vector <int> carte ; vector < Vector > normales_tri ;
vector <vector <double> > datas ;

chem=chemin ;
chem.append("-ForceTotale.txt") ;

out.open(chem.c_str(), ios::out) ;

loopdat(loop) ; 

cout << "\nFDump::write_forcestot          " ;
actions.total=loop[2]-loop[0] ; actions.disp_progress() ;
if (actions["group-triangles"].set)
{
chem2=chemin ; chem2.append("-ForceByElement-group-nomean.txt") ;
groupout.open(chem2.c_str()) ;
}

for (i=loop[0] ; i<loop[2] ; i+=loop[1])
    {
	actions.valeur=i ;
    check_timestep(i) ;
    steps[i].mean_forces(forces) ;
    out<<forces[0] << " " << forces[1] << " " << forces[2] << "\n" ;

    if (actions["group-triangles"].set)
     {
     int idx[3] ; bool counted=false ; double tolerance=0.01 ;
     int k,j ;
     Vector normale ;
     vector <int> nb_group_tri ; //Matrix centrebis(3,1) ;
     idx[0]=steps[i].find_idx(IDS("FORCEX")) ; idx[1]=steps[i].find_idx(IDS("FORCEY")) ; idx[2]=steps[i].find_idx(IDS("FORCEZ")) ;

     // S'il s'agît de la première boucle, on réalise la carte avant tout
     if (i==loop[0])
        {
    	nb_tri=0 ;
        for (k=0 ; k<steps[i].nb_triangles ; k++)
          {
          //normale=steps[i].get_tri_normal(k) ;
          normale[0]=steps[i].datas[steps[i].find_idx(IDS("NORMALEX"))][k] ;
          normale[1]=steps[i].datas[steps[i].find_idx(IDS("NORMALEY"))][k] ;
          normale[2]=steps[i].datas[steps[i].find_idx(IDS("NORMALEZ"))][k] ;
          counted=false ;

          for (j=0 ; j<normales_tri.size() ; j++)
            {
    	    if (Calcul::norm(normales_tri[j]-normale) < tolerance)
    	       {
    	       carte.push_back(j) ; counted=true ;
    	       //datas[0][j]+=steps[i].datas[idx[0]][k] ; datas[1][j]+=steps[i].datas[idx[1]][k] ; datas[2][j]+=steps[i].datas[idx[2]][k] ;
    		   //nb_group_tri[j]++ ; counted=true ;
               //surface[j]+=steps[i].get_tri_surface(k);
    	       }
            }

          if (counted==false)
           {
           nb_tri++ ; normales_tri.push_back(normale) ;
           carte.push_back(nb_tri-1) ;
           //datas[0].push_back(steps[i].datas[idx[0]][k]) ; datas[1].push_back(steps[i].datas[idx[1]][k]) ; datas[2].push_back(steps[i].datas[idx[2]][k]) ;
           //nb_group_tri.push_back(1) ;
           //surface.push_back(steps[i].get_tri_surface(k));
           }
          }
        // On redimensionne la structure qui accueillera les données moyennées
        datas.resize(5) ; datas[0].resize(nb_tri) ; datas[1].resize(nb_tri) ; datas[2].resize(nb_tri) ; datas[3].resize(nb_tri) ; datas[4].resize(nb_tri) ;
        }
     // Dans tous les cas, on effectue la somme et le moyennage en se basant sur la carte créée
     for (j=0 ; j<5 ; j++) { for (k=0; k<nb_tri ; k++) { datas[j][k]=0 ;}}
     for (k=0 ; k<steps[i].nb_triangles ; k++)
        {
    	 datas[0][carte[k]]+=steps[i].datas[idx[0]][k] ; datas[1][carte[k]]+=steps[i].datas[idx[1]][k] ; datas[2][carte[k]]+=steps[i].datas[idx[2]][k] ;
    	 datas[3][carte[k]]+=steps[i].get_tri_surface(k) ;
    	 datas[4][carte[k]]++ ;
        }

     // Writing ;
     double angle ;
     for (k=0 ; k<normales_tri.size() ; k++)
        {
    	angle=Calcul::arctan(normales_tri[k](1),normales_tri[k](3)) ;
    	groupout << normales_tri[k](1) << '\t' << normales_tri[k](2) << '\t' << normales_tri[k](3) << '\t' ;
    	groupout << angle << '\t' << datas[3][k] << '\t' ;
    	groupout << datas[0][k] << '\t' << datas[1][k] << '\t' << datas[2][k] << '\t' ;
    	groupout << datas[4][k] << '\t' ;
        }
    }
   }

if (actions["group-triangles"].set)
{
char chaine[500] ;
groupout.close() ;
sprintf(chaine, "Le fichier %s a ete cree. La commande matlab de redimensionnement est reshape(X,9,%d,%ld).", chem2.c_str(), nb_tri,((loop[2]-loop[0])/loop[1])) ;
DISP_Info(chaine) ;
}

out.close() ;
return 0 ;
}
//------------------------------------------------
int FDump::write_coupletot (string chemin)
{
string chem ; int i;
ofstream out ;
Vector couple ;
Vector c ;
long int loop[3] ; 

chem=chemin ;
chem.append("-CoupleTotal.txt") ;

out.open(chem.c_str(), ios::out) ;

loopdat(loop) ; 

cout << " Entrer les coordonnées du centre du cylindre (x y z) : \n" ;
cin >> c(1) >> c(2) >> c(3) ;
cout << "Coordonnées entrées : " << c(1) << " " << c(2) << " " << c(3) << "\n";
cout << "\nFDump::write_coupletot          " ;
actions.total=loop[2]-loop[0] ; actions.disp_progress() ;

for (i=loop[0] ; i<loop[2] ; i+=loop[1])
    {actions.valeur=i ;
    check_timestep(i) ;
    steps[i].Fcouple(couple, c) ;
    out<<couple(1) << " " << couple(2) << " " << couple(3) << "\n" ;
    }
out.close() ;
return 0 ;
}
//----------------------------------------
int FDump::mean_stress (string chemin)
{
string chem, chem2 ; int i,j,k ;
ofstream out ;
int idx[9] ;
long int compteur ;
char filterdo[30] ; long int loop[3] ; 
Step tmp_forces ; double anglebase ;

chem=chemin ;
chem.append("-ForceByElement.vtk") ;

loopdat(loop) ; 

cout << "\nFDump::mean_stress          " ;
actions.total=loop[2]-loop[0] ; actions.disp_progress() ;
check_timestep(loop[0]) ;
tmp_forces=steps[loop[0]] ;
if (actions["anglebystep"].set)
   {
   anglebase=actions["anglebystep"]["anglebystep"]/180.*M_PI ;
   // Calcul des centres des triangles :
   DISP_Err("NON IMPLEMENTé : il faut reprogrammer anglebystep dans FDump::mean_stress") ;
   /*Matrix mat_null(3,1) ; mat_null=0.0 ;
   centres.resize(tmp_forces.nb_triangles,mat_null) ;
   for (j=0 ; j<tmp_forces.nb_triangles ; j++)
           {centres[j]=steps[debut].get_tri_center(j) ;}*/
  }
double angle ;
compteur=1 ;
for (i=loop[0]+1 ; i<loop[2] ; i+=loop[1])
    {actions.valeur=i-loop[0] ;

    if (actions["anglebystep"].set)
       {
       // Remise en place des positions des points
       // On suppose une rotation par rapport à l'axe vertical passant par x=0,y=0
       angle=-anglebase*(i-loop[0]) ;
       Filter filtre_tmp ;
       // Création du filtre pour le Ldump
       sprintf(filterdo, "null::zrotate::%.10f", angle) ;
       DISP_Err("Please check the type value here (anglebystep in FDump.cpp, meanstress)\n") ; 
       postfiltre=filtre_tmp.parse_arg(filterdo, TF ) ;
       postfiltre[0].valeur=angle ;
       }

    check_timestep(i) ;

    if (actions["anglebystep"].set)
       {
    	   DISP_Err("NON IMPLEMENTé : il faut reprogrammer anglebystep dans FDump::mean_stress") ;
       /*for (j=0 ; j<steps[i].nb_triangles ; j++)
           {
           centrebis=steps[i].get_tri_center(j) ;
           if (fabs(centres[j](1,1)-centrebis(1,1))>0.000001 || fabs(centres[j](2,1)-centrebis(2,1))>0.000001 || fabs(centres[j](3,1)-centrebis(3,1))>0.000001)
             {  printf("%f\n", centres[j](1,1)-centrebis(1,1)) ; fflush(stdout) ; cout << "WARN/ERR : les numéros des triangles ont changé !!\n" ; }
           }*/
       }


    for (j=0 ; j<steps[i].datas.size() ; j++)
      {
      for (k=0 ; k<steps[i].datas[j].size() ; k++)
          {tmp_forces.datas[j][k]=tmp_forces.datas[j][k]+steps[i].datas[j][k] ;}
      }
    compteur++ ;
    }
postfiltre.clear() ;
//check_timestep(debut) ;
for (j=0 ; j<tmp_forces.datas.size() ; j++)
 {
 for (k=0 ; k<tmp_forces.datas[j].size() ; k++)
     {tmp_forces.datas[j][k]=tmp_forces.datas[j][k]/compteur ; }
 }

if (actions["group-triangles"].set)
 {
 int idx[9] ; bool counted=false ; double tolerance=0.01 ;
 vector < Vector > normales_tri ; vector <vector <double> > datas ; vector <double> surface ;
 vector <int> nb_group_tri ;
 datas.resize(3) ;
 idx[0]=tmp_forces.find_idx(IDS("FORCEX"))   ; idx[1]=tmp_forces.find_idx(IDS("FORCEY"))   ; idx[2]=tmp_forces.find_idx(IDS("FORCEZ")) ;
 idx[3]=tmp_forces.find_idx(IDS("CENTREX"))  ; idx[4]=tmp_forces.find_idx(IDS("CENTREY"))  ; idx[5]=tmp_forces.find_idx(IDS("CENTREZ")) ;
 idx[6]=tmp_forces.find_idx(IDS("NORMALEX")) ; idx[7]=tmp_forces.find_idx(IDS("NORMALEY")) ; idx[8]=tmp_forces.find_idx(IDS("NORMALEZ")) ;

 for (i=0 ; i<tmp_forces.nb_triangles ; i++)
   {
	//DISP_Err("NON IMPLEMENTé : il faut reprogrammer group-triangles dans FDump::mean_stress") ;
   //normale=tmp_forces.get_tri_normal(i) ;
   Vector normale (tmp_forces.datas[idx[6]][i], tmp_forces.datas[idx[7]][i], tmp_forces.datas[idx[8]][i]) ;
   counted=false ;

   //centrebis=tmp_forces.get_tri_center(i) ;
   //out << centrebis(1,1) << '\t' << centrebis(2,1) << '\t' << centrebis(3,1) << '\t' ;

   for (j=0 ; j<normales_tri.size() ; j++)
     {
	 if (Calcul::norm(normales_tri[j]-normale) < tolerance)
	    {
		 datas[0][j]+=tmp_forces.datas[idx[0]][i] ; datas[1][j]+=tmp_forces.datas[idx[1]][i] ; datas[2][j]+=tmp_forces.datas[idx[2]][i] ;
		 nb_group_tri[j]++ ; counted=true ;
         surface[j]+=tmp_forces.get_tri_surface(i);
	    }
     }

   if (counted==false)
     {
	 datas[0].push_back(tmp_forces.datas[idx[0]][i]) ; datas[1].push_back(tmp_forces.datas[idx[1]][i]) ; datas[2].push_back(tmp_forces.datas[idx[2]][i]) ;
     normales_tri.push_back(normale) ; nb_group_tri.push_back(1) ;
     surface.push_back(tmp_forces.get_tri_surface(i));
     }
   }

 // Writing ;
 chem2=chemin ; chem2.append("-ForceByElement-group.txt") ;
 out.open(chem2.c_str()) ;
 double angle ;
 for (i=0 ; i<normales_tri.size() ; i++)
    {
	angle=Calcul::arctan(normales_tri[i](1),normales_tri[i](3)) ;
	out << normales_tri[i](1) << '\t' << normales_tri[i](2) << '\t' << normales_tri[i](3) << '\t' ;
	out << angle << '\t' << surface[i] << '\t' ;
	out << datas[0][i] << '\t' << datas[1][i] << '\t' << datas[2][i] << '\t' ;
	out << nb_group_tri[i] << '\n' ;
    }
 out.close() ;
 }

out.open(chem.c_str()) ;
tmp_forces.write_asVTK(out) ;
return 0 ;
}
