#include "Headers/LDump.h"

// =====================================================
// Fonctions de la classe LDump ========================
//======================================================
void LDump::disp(void)
{
Dump::disp() ;
check_timestep(0) ; steps[0].disp() ;
cout << "------\n" ;
check_timestep(nbsteps-1) ; steps[nbsteps-1].disp() ;
cout <<"--------------------\n\n" ;
}
//------------------------------------------------------
int LDump::write_asDUMP (string chemin)
{
long int i ; long int loop[2] ; 
ofstream out ;
stringstream chemin2 ;

loopdat(loop) ;

actions.total=(loop[2]-loop[0]) ;
actions.valeur=0 ;
cout << "\nLDump::write_asDUMP          " ;
actions.disp_progress() ;
chemin2.str("") ;
chemin2 << chemin << ".ds" ;
out.open(chemin2.str().c_str(), ios::out) ;

for (i=loop[0] ; i<loop[2] ; i+=loop[1])
   {actions.valeur=(i-loop[0]) ;
   check_timestep(i) ;
   steps[i].write_asDUMP(out) ;
   }
out.close() ;
return 0 ;
}
//------------------------------------
int LDump::write_asRESTART(string chemin)
{
long int i ;
ofstream out ;
stringstream chemin2 ;

cout << "\nLDump::write_asRESTART          " ;
chemin2.str("") ;
chemin2 << chemin << ".data" ;
out.open(chemin2.str().c_str(), ios::out) ;
check_timestep(actions["dump2restart"]["timestep"]) ;
steps[actions["dump2restart"]["timestep"]].Lwrite_asRESTART(out) ;
out.close() ;
return 0 ;
}

// =====================================================
// Fonctions de la sous classe LucDump =================
//======================================================
int LucDump::read (unsigned char read_type, int index)
{
// read_type==1 : sauver les données lues
 //* read_type==2 : ne lire qu'un dump
 //* read_type==4 : creer la strcuture step et incrémenter nb_step

  string ligne;
  Step null_step ;
  int i, j, type, stop=0, idx ;
  char *ligneC, *actu,*suiv ; 
  bool istriclinic=false ; static bool first_naninf=true ; 

  // Function pointer readline(ifstream dumpin, string ligne)
  // if compress => the function decompress and write the ligne string
  // else: the function is just a wrapper for getline
  // idem 
  
  ligneC=(char*)malloc(5000) ; 
  Type=TL ; 
  
  if (read_type & 4) {nbsteps=0; printf("Lecture initiale du dump, peut prendre du temps ...\n") ; }
  if (read_type & 2) idx=index ;

  while (!dumpin.eof() && stop<2 )
  {
   dumpin.getline(ligne) ; 
   type=isitem(ligne) ; 
   switch (type)
   {
     case 1 : if (read_type & 4)
		{
		 nbsteps++ ; // Create new step and read 1 line to have timestep number.
		 steps.resize(nbsteps, null_step) ;
		 steps[nbsteps-1].multisphere=&multisphere ;
		 steps[nbsteps-1].nb_atomes=0 ; steps[nbsteps-1].Type=TL ;
		 steps[nbsteps-1].posinfile=dumpin.tellg()-(streampos)1 -(streampos)ligne.length() ;
                 if (index==-1) idx=nbsteps-1 ;
                 actions.valeur=dumpin.tellg() ;
		}
	      if (read_type & 2) {stop++ ; if (stop>=2) break ; }
	      //if (read_type & 1) { }
	      dumpin >> steps[idx].timestep ;
	      break ;

     case 2 : if (read_type & 1) {dumpin >> steps[idx].nb_atomes ; }
	      break ;

     case 3 : if (read_type & 1)
	      {
	       // Taking care of triclinic boxes ...
	       istriclinic = (ligne.find("xy") != string::npos) ; 
	       dumpin >> steps[idx].box[0][0] >> steps[idx].box[0][1] ; if (istriclinic) dumpin >> steps[idx].triclinic[0] ; 
	       dumpin >> steps[idx].box[1][0] >> steps[idx].box[1][1] ; if (istriclinic) dumpin >> steps[idx].triclinic[1] ; 
	       dumpin >> steps[idx].box[2][0] >> steps[idx].box[2][1] ; if (istriclinic) dumpin >> steps[idx].triclinic[2] ; 
	      }
	      break ;

     case 4 : if (read_type & 1)
	      {
	       sparselabels(steps[idx], ligne) ; // Nouvelles données atomes. On continue la lecture de la ligne et on écrit la liste d'index
	       // On alloue le tableau de données nécessaire
	       steps[idx].datas.resize(steps[idx].nb_idx, vector<double>(1,0)) ;
	       for (i=0 ; i<steps[idx].nb_idx; i++)
	       {steps[idx].datas[i].resize(steps[idx].nb_atomes, 0) ;}
	       // On lit l'ensemble des données que l'on range
	       // TEST Optimized version
	       for (i=0 ; i<steps[idx].nb_atomes ; i++)
	       {
		dumpin.getline(&ligneC,5000) ;
	        actu=ligneC ; 
		for (j=0 ; j<steps[idx].nb_idx ; j++)
		{
		  steps[idx].datas[j][i]=strtod(actu, &suiv) ;
		  if (actu==suiv) 
		  { 
		    DISP_Err("Conversion to double failed dans LDump! Something really bad happened. Going to next field...\n") ; 
		    while (*actu!=' ' && *actu!=0) actu++ ; 
		  }
		  else 
		    actu=suiv ; 
		  if ( !isfinite(steps[idx].datas[j][i]) && first_naninf)
		  {
		    DISP_Warn("\n\n!!! Une valeur lue dans LDump est inf ou nan. ") ;  
		    DISP_Warn("Elle n'a PAS été supprimée.\n") ;
                    first_naninf=false ; 
		      //steps[idx].datas[j][i]=0 ;
		  }
		}
	       }
	       
	       
	       /*for (i=0 ; i<steps[idx].nb_atomes ; i++)
	       {
		for (j=0 ; j<steps[idx].nb_idx ; j++)
		{
		  dumpin >> steps[idx].datas[j][i]; 
		  if (dumpin.fail())
			{
			 int st=dumpin.rdstate() ;
			 dumpin.clear() ; 
			 DISP_Warn("Une erreur de lecture (L) est apparue à l'octet : ") ; 
			 printf("%d (error number %d, ts %d). ", (int)dumpin.tellg(), st, index) ;
			 DISP_Warn("Elle a été supprimée et la valeur est mise à 0.\n") ;
			 steps[idx].datas[j][i]=0 ; 
			}
		}
               }*/
	     }
     	 break ;
   }
  }
dumpin.clear() ;
return 0 ;
}
//------------------------------------------
bool LucDump::operator == (Dump &dmp)
{
int i ;

if (nbsteps!=dmp.nbsteps) return false ;

for (i=0 ; i<nbsteps ; i++)
  {
  check_timestep(i) ;
  dmp.check_timestep(i) ;
  if (!(steps[i]==dmp.steps[i]))
     return false ;
  }
return true ;
}
//-------------------------------------------------
int LucDump::isitem (string ligne)
{
size_t position ;

position=ligne.find("ITEM: ") ;
if (position!=string::npos)
  {
  position=ligne.find("TIMESTEP") ;
  if (position!=string::npos)
    return 1 ;
  position=ligne.find("NUMBER OF ATOMS") ;
  if (position!=string::npos)
    return 2 ;
  position=ligne.find("BOX BOUNDS") ;
  if (position!=string::npos)
    return 3 ;
  position=ligne.find("ATOMS") ;
  if (position!=string::npos)
   return 4 ;

  cerr << "WARN1 : item inconnu " << ligne << ".\n" ;
  return -1 ;
  }
else
   return 0 ;
}
//-------------------------------------------------------------
int LucDump::sparselabels(Step &step , string ligne)
{
  int compteur=0, res ; // i=12 pour virer le "ITEM: ATOMS " en début de chaîne) ;
  size_t espace ;
  string ligne2, item ;

  ligne2=ligne.substr(12) ;
  while (1)
  {
    espace=ligne2.find_first_of(' ') ;
    if (espace==ligne2.npos) break ;
    item=ligne2.substr(0, espace) ;
    ligne2=ligne2.substr(espace+1) ;

    step.idx_col.resize(compteur+1) ;
    res=IDS(item) ;
    
    if (res==-1) 
    {DISP_Warn ("Un type de donné n'est pas référencé, IDS va l'ajouter\n") ; 

     res=IDS.new_id(item, TL) ; 
     if (res!=-1) step.idx_col[compteur]=res ;
     else 
     {
       DISP_Err("Le type inconnu n'a pas pu être ajouté !!\n"); 
       step.idx_col[compteur]=IDS("UNKNOWN") ; 
     }
    }
    else
     step.idx_col[compteur]=res ;
    
    compteur++ ;
  }
  step.nb_idx=compteur ;
  return compteur ;
}
    
    /*if (item=="id")            step.idx_col[compteur]=ID ;
    else if (item == "x")      step.idx_col[compteur]=POSX ;
    else if (item == "y")      step.idx_col[compteur]=POSY ;
    else if (item == "z")      step.idx_col[compteur]=POSZ ;
    else if (item == "xu")      {DISP_Warn("les item xu sont traites comme des x\n") ;  step.idx_col[compteur]=POSX ;}
    else if (item == "yu")      {DISP_Warn("les item yu sont traites comme des y\n") ; step.idx_col[compteur]=POSY ;}
    else if (item == "zu")      {DISP_Warn("les item zu sont traites comme des z\n") ; step.idx_col[compteur]=POSZ ;}
    else if (item == "vx")      step.idx_col[compteur]=VX ;
    else if (item == "vy")      step.idx_col[compteur]=VY ;
    else if (item == "vz")      step.idx_col[compteur]=VZ ;
    else if (item == "fx")      step.idx_col[compteur]=FX ;
    else if (item == "fy")      step.idx_col[compteur]=FY ;
    else if (item == "fz")      step.idx_col[compteur]=FZ ;
    else if (item == "type")    step.idx_col[compteur]=TYPE ;
    else if (item == "radius")  step.idx_col[compteur]=RAYON ;
    else if (item == "mass")    step.idx_col[compteur]=MASSE ;
    else if (item == "omegax")  step.idx_col[compteur]=OMEGAX ;
    else if (item == "omegay")  step.idx_col[compteur]=OMEGAY ;
    else if (item == "omegaz")  step.idx_col[compteur]=OMEGAZ ;
    else if (item == "sigxx")   step.idx_col[compteur]=SIGMAXX ;
    else if (item == "sigxy")   step.idx_col[compteur]=SIGMAXY ;
    else if (item == "sigxz")   step.idx_col[compteur]=SIGMAXZ ;
    else if (item == "sigyx")   step.idx_col[compteur]=SIGMAYX ;
    else if (item == "sigyy")   step.idx_col[compteur]=SIGMAYY ;
    else if (item == "sigyz")   step.idx_col[compteur]=SIGMAYZ ;
    else if (item == "sigzx")   step.idx_col[compteur]=SIGMAZX ;
    else if (item == "sigzy")   step.idx_col[compteur]=SIGMAZY ;
    else if (item == "sigzz")   step.idx_col[compteur]=SIGMAZZ ;
    else if (item == "f_force_cyl[1]") step.idx_col[compteur]=FORCEWALLX ;
    else if (item == "f_force_cyl[2]") step.idx_col[compteur]=FORCEWALLY ;
    else if (item == "f_force_cyl[3]") step.idx_col[compteur]=FORCEWALLZ ;
    else {cerr << "WARN2 : le type de données" << item <<"est inconnu\n" ; step.idx_col[compteur]=UNKNOWN ; }*/

//---------------------------------------------------
int LucDump::write_forcestot (string chemin)
{
string chem ; int i;
int nb_tri=-1; long int loop[3] ;
ofstream out ;
double forces[3] ;

chem=chemin ;
chem.append("-ForceTotale.txt") ;

out.open(chem.c_str(), ios::out) ;

loopdat(loop) ;

cout << "\nLucDump::write_forcestot          " ;
actions.total=loop[2]-loop[0] ; actions.disp_progress() ;

for (i=loop[0] ; i<loop[2] ; i+=loop[1])
    {
    actions.valeur=i ;
    check_timestep(i) ;
    steps[i].mean_forces(forces) ;
    out<<forces[0] << " " << forces[1] << " " << forces[2] << "\n" ;
   }

out.close() ;
return 0 ;
}
//------------------------------------------------
int LucDump::write_wallforce (string chemin)
{
string chem ;
char infos[500] ;
long int i, loop[3];
double ** meanforces ; int *meangrains ;
ofstream out ;
//double forces[3] ;

chem=chemin ; chem.append("-ForceByAngle.txt") ;
out.open(chem.c_str(), ios::out) ;

loopdat(loop) ;

cout << "\nLucDump::write_wallforce          " ;
actions.total=loop[2] ; actions.disp_progress() ;

// Création des tableaux de données moyennes si nécessaires
if (actions["mean"].set)
  {
  meanforces=(double **)malloc((int)(actions["wallforce-by-angle"]["nbbox_theta"])*sizeof(double *)) ;
  meangrains=(int *)malloc((int)(actions["wallforce-by-angle"]["nbbox_theta"])*sizeof(int)) ;
  for (i=0 ; i<actions["wallforce-by-angle"]["nbbox_theta"] ; i++)
     {meanforces[i]=(double *)malloc(7*sizeof(double)) ; meanforces[i][0]=meanforces[i][1]=meanforces[i][2]=0 ; meangrains[i]=0 ;
     meanforces[i][3]=0 ;meanforces[i][4]=meanforces[i][5]=meanforces[i][6]=0 ;}
  }

// Boucle sur les ts
for (i=loop[0] ; i<loop[2] ; i+=loop[1])
    {
    actions.valeur=i ;
    check_timestep(i) ;
    steps[i].wall_force(out, meanforces, meangrains) ;
    }
out.close() ;
double denom=(loop[2]-loop[0])/(double)loop[1] ;
sprintf(infos, "La largeur en angle des boites est de %f degrés (%f radians).\n", 360.0/actions["wallforce-by-angle"]["nbbox_theta"], 2*M_PI/actions["wallforce-by-angle"]["nbbox_theta"]) ;    ;
DISP_Info(infos) ;
sprintf(infos, "La commande matlab à utiliser est reshape(X, %d, %d, %d).", 9, (int)actions["wallforce-by-angle"]["nbbox_theta"], (int)denom) ;
DISP_Info(infos) ;
// Création du fichier avec moyennage si besoin
if (actions["mean"].set)
 {

 chem=chemin ; chem.append("-ForceByAngle-mean.txt") ;
 out.open(chem.c_str(), ios::out) ;
 for (i=0; i<actions["wallforce-by-angle"]["nbbox_theta"] ; i++ )
    {
	out << (2*i*M_PI)/actions["wallforce-by-angle"]["nbbox_theta"]+M_PI/actions["wallforce-by-angle"]["nbbox_theta"] << " "  << meangrains[i]/denom << " " ;
	out << meanforces[i][0]/denom << " " << meanforces[i][1]/denom << " " << meanforces[i][2]/denom << " " << meanforces[i][3]/denom ;
	out << " " << meanforces[i][4]/denom << " " << meanforces[i][5]/denom << " " << meanforces[i][6]/denom <<"\n" ;
    }
 out.close() ;
 }
return 0 ;
}
//==========================================
int LucDump::write_xray (string chemin)
{
long int loop[3] ; string chem ; char num[50] ; 
int i, dir ; int width, height ; 
double box[6] ; double **img ; 

loopdat(loop) ;
cout << "\nLucDump::write_xray          " ;
actions.set_progress(loop) ; actions.disp_progress() ; 

// Box setting up
check_timestep(loop[0]) ; 
i=loop[0] ; 
if (! actions["use-box"].set)
{
  actions["use-box"].manual_set("box_xmin",-1e15) ; actions["use-box"].manual_set("box_xmax",1e15) ; 
  actions["use-box"].manual_set("box_ymin",-1e15) ; actions["use-box"].manual_set("box_ymax",1e15) ; 
  actions["use-box"].manual_set("box_zmin",-1e15) ; actions["use-box"].manual_set("box_zmax",1e15) ; 
}
// No usebox in the direction of averaging ; 
dir=(int)actions["xray"]["dir"] ;
switch(dir) {
  case 0:     
    box[0]=fmax(actions["use-box"]["box_ymin"], steps[i].box[1][0]) ; 
    box[1]=fmin(actions["use-box"]["box_ymax"], steps[i].box[1][1]) ; 
    box[2]=fmax(actions["use-box"]["box_zmin"], steps[i].box[2][0]) ; 
    box[3]=fmin(actions["use-box"]["box_zmax"], steps[i].box[2][1]) ; 
    box[4]=steps[i].box[0][0] ; box[5]=steps[i].box[0][1] ;
    break ; 
  case 1:     
    box[0]=fmax(actions["use-box"]["box_xmin"], steps[i].box[0][0]) ; 
    box[1]=fmin(actions["use-box"]["box_xmax"], steps[i].box[0][1]) ; 
    box[2]=fmax(actions["use-box"]["box_zmin"], steps[i].box[2][0]) ; 
    box[3]=fmin(actions["use-box"]["box_zmax"], steps[i].box[2][1]) ; 
    box[4]=steps[i].box[1][0] ; box[5]=steps[i].box[1][1] ;
    break ; 
  case 2: 
    box[0]=fmax(actions["use-box"]["box_xmin"], steps[i].box[0][0]) ; 
    box[1]=fmin(actions["use-box"]["box_xmax"], steps[i].box[0][1]) ; 
    box[2]=fmax(actions["use-box"]["box_ymin"], steps[i].box[1][0]) ; 
    box[3]=fmin(actions["use-box"]["box_ymax"], steps[i].box[1][1]) ; 
    box[4]=steps[i].box[2][0] ; box[5]=steps[i].box[2][1] ; 
    break ; 
  default: DISP_Warn("Direction de moyennage xray inconnue\n") ; 
} 
width=actions["xray"]["width"] ; height=actions["xray"]["height"] ; 
img=(double**)malloc(sizeof(double*)*width) ;
for (i=0 ; i<width ; i++) img[i]=(double*)malloc(sizeof(double)*height) ; 

for (i=loop[0] ; i<loop[2] ; i+=loop[1])
  {
  actions.valeur=i ;
  check_timestep(i) ;   
  steps[i].xray_projection(dir, width, height, img, box) ;
  
  sprintf(num, "%04d", (i-loop[0])/loop[1]) ;
#ifdef USETIFF
  chem=chemin ; chem.append("-") ; chem.append(num) ; chem.append(".tif") ; 
  Writing::TIFF_bw_write(chem, width, height, img) ;
#else
  printf("USETIFF not set when compiled => cannot export to tiff") ; 
#endif
  }

for (i=0 ; i<width ; i++) free(img[i]) ; 
free(img) ;
return 0 ; 
}
//---------------------------------------
int LucDump::write_multisphere_dumbell (string chemin)
{
  /*int Meridien=10, Parallel=10, *Numbers; // Hard coded sampling. On va travailler que sur la demi-sphere +x, y, z ; de pôles z
  double dtheta, dphi ; int dt, dp ;
  Numbers=(int*)calloc(Meridien*Parallel,sizeof(int)) ;  
  dtheta=M_PI/Parallel ; 
  dphi=M_PI/Meridien ; */
  Icosahedre Ico ; 
  //if (actions["mean"].set) Ico.subdivide(4) ; 
  
  long int loop[3] ; int idx[5] ; string chem ; Matrix3d K, Kvec, Kmatseg ; Vector3d Kval ; Map<Vector3d> Ksegment(NULL); int Kn ; 
  int i, j, k, l, n ; double maxlen ; int idmax ; FILE *out,*out2 ; 
  int *id1, *id2 ; int gp ; Vector v, vsph, null(0) ; double theta=5, angular[72], angular2[72], nangular, tmptheta ; int tmpidxtheta ;
  int **gps, ngp=-1 ; vector < Vector > pts, segments ; Vector t, dir, vect, centroid ; int longest ; double Phi ;
  int type ; type=actions["multisphere"]["type"] ; double box[6] ; 
  bool symetrie[7]={false,false,false,false,false,false,false} ; 
  
  loopdat(loop) ; 
  cout << "\nLucDump::multisphere_dumbell          \n" ;
  actions.set_progress(loop) ; actions.disp_progress() ; 
  check_timestep(loop[0]) ; 
  idx[0]=steps[loop[0]].find_idx(IDS("POSX")) ; idx[1]=steps[loop[0]].find_idx(IDS("POSY")) ;  idx[2]=steps[loop[0]].find_idx(IDS("POSZ")) ;  
  idx[3]=steps[loop[0]].find_idx(IDS("IDMULTISPHERE")) ; idx[4]=steps[loop[0]].find_idx(IDS("ID")) ; 
  
  chem=chemin ; 
  chem.append("-") ; chem.append("multisphere-tensor") ; chem.append(".txt") ; 
  out=fopen(chem.c_str(), "w") ; 
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
   
  if (actions["symetriser"].set)
  {
    symetrie[0]=true ; double r=actions["symetriser"]["axes"] ;
    if (r>=100) {symetrie[1]=true ; r=r-100 ; }
    if (r>=10)  {symetrie[2]=true ; r=r-10 ; }
    if (r>=1)   {symetrie[3]=true ;}
  }
  printf("%d %d %d %d %d %d %d-----------------\n", symetrie[0], symetrie[1],symetrie[2],symetrie[3],symetrie[4],symetrie[5],symetrie[6]) ; 
  
  for (i=0 ; i<360/theta ; i++) angular[i]=angular2[i]=0 ; 
  theta=theta/180*M_PI ; 
  
  gps=(int **)malloc(sizeof(int*)*0) ;
  for (i=0 ; i<steps[loop[0]].nb_atomes ; i++)
  {
    if (steps[loop[0]].datas[idx[3]][i]>=0)
    {
      if (steps[loop[0]].datas[idx[3]][i]>ngp)
      {
	gps=(int**)realloc(gps,sizeof(int*)*(steps[loop[0]].datas[idx[3]][i]+1)) ; 
	for (j=ngp+1 ; j<=steps[loop[0]].datas[idx[3]][i] ; j++)
	{
	  gps[j]=(int*)malloc(sizeof(int)*1) ; 
	  gps[j][0]=0 ; 
	}
      ngp=steps[loop[0]].datas[idx[3]][i] ; 
      }
      gps[(int)(steps[loop[0]].datas[idx[3]][i])][0]++ ; 
      gps[(int)(steps[loop[0]].datas[idx[3]][i])]=(int*)realloc(gps[(int)(steps[loop[0]].datas[idx[3]][i])], sizeof(int)*(gps[(int)(steps[loop[0]].datas[idx[3]][i])][0]+1)) ;
      gps[(int)(steps[loop[0]].datas[idx[3]][i])][gps[(int)(steps[loop[0]].datas[idx[3]][i])][0]]=(int)steps[loop[0]].datas[idx[4]][i] ;
    }
  }
  if (ngp==-1) {DISP_Warn("Aucun groupe multisphere trouvé, il y a un problème.") ; printf("%d %d %d", loop[0], steps[loop[0]].nb_atomes, idx[3]) ; fflush(stdout) ;   } 
  
  for (i=1, longest=0 ; i<=ngp ; i++) if (longest<gps[i][0]) longest=gps[i][0] ; 
  pts.resize(longest, null) ; 
  segments.resize(longest*(longest-1)/2, null) ; 
  bool * activegp = new bool[ngp] ;
  std::fill_n(activegp, ngp, true); 
  
  int r ;  int count=0 ;
  for (i=loop[0] ; i<loop[2] ; i+=loop[1])
  {
    actions.valeur=i ; 
    check_timestep(i) ;   
    K=Matrix3d::Zero() ; Kn=0 ; 
    for (j=1 ; j<=ngp ; j++)
    {
     centroid=0 ; 
     if (!activegp[j]) continue ; 
     for (k=0 ; k<gps[j][0] ; k++)
     {
       if (steps[i].datas[idx[4]][gps[j][k+1]-1]!=gps[j][k+1]) {printf("%g %g ", steps[i].datas[idx[4]][gps[j][k+1]-1],gps[j][k+1] ) ; DISP_Err("Probleme in multisphere\n") ;} 
       t.set(steps[i].datas[idx[0]][gps[j][k+1]-1], steps[i].datas[idx[1]][gps[j][k+1]-1], steps[i].datas[idx[2]][gps[j][k+1]-1]);
       pts[k]=t ;
       if (t.isnan()) 
       {
	 activegp[j]=false ; printf("Le groupe %d a été perdu. Atomes:", j) ; for (l=0 ; l<gps[j][0] ; l++) {printf("%d ", gps[j][l+1]-1) ; } printf("\n") ; break ;  
       }
       centroid=centroid+pts[k] ;
     }
     if (!activegp[j]) continue ; 
     centroid=centroid/4 ; 
     if (symetrie[0]) 
     {
       if (symetrie[1]==true && centroid(1)<0) {symetrie[4]=true ; centroid(1)=-centroid(1) ; } else {symetrie[4]=false ; }
       if (symetrie[2]==true && centroid(2)<0) {symetrie[5]=true ; centroid(2)=-centroid(2) ; } else {symetrie[5]=false ; }
       if (symetrie[3]==true && centroid(3)<0) {symetrie[6]=true ; centroid(3)=-centroid(3) ; } else {symetrie[6]=false ; }
     }
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
     if (symetrie[0])
     {
       if (symetrie[4]) {segments[idmax](1)=-segments[idmax](1) ; }
       if (symetrie[5]) {segments[idmax](2)=-segments[idmax](2) ; }
       if (symetrie[6]) {segments[idmax](3)=-segments[idmax](3) ; }
     }
     
     
     vsph=Geometrie::cart2sph(segments[idmax]) ; 
     if (type==0)
       {if (vsph(1)>actions.Cst["Radius"]*gps[j][0]*2+0.0001) {continue ;}}
     else if (type==1)
       {if (vsph(1)>actions.Cst["Radius"]*(floor(log2((gps[j][0]-1)/3.))*2)+0.0001) {continue ;}}
       
     tmptheta=Calcul::angle_0_2pi(atan2(segments[idmax](3), segments[idmax](1))) ; 
     if (tmptheta>(2*M_PI-theta/2)) tmptheta-=(2*M_PI) ; 
     tmpidxtheta=(int)(round(tmptheta/theta)) ; 
     
     if (tmpidxtheta>=0 && tmpidxtheta<=71)
     {
     angular[tmpidxtheta]++ ; 
     nangular++ ; 
     }
     else
       printf("!\n") ; 
     
     for (k=0 ; k<72 ; k++)
     {
      vect(1)=cos(2*M_PI/72.*k) ; vect(2)=0 ; vect(3)=sin(2*M_PI/72.*k) ; 
      angular2[k]+=(segments[idmax].dot(vect))*(segments[idmax].dot(vect)) ;  
     }
     new (&Ksegment) Map<Vector3d>(segments[idmax].dat); // THIS IS NOT AN ALLOCATION (no delete) ; 
     Ksegment=Ksegment/(Ksegment.norm()) ; 
     Kmatseg=Ksegment*(Ksegment.transpose());
     K=K+Kmatseg ; Kn++ ; 

     r=Ico.belonging_tri(segments[idmax]) ;  
     if (r!=-1) Ico.data[r]=Ico.data[r]+1 ; else printf("!") ; //printf("%g %g %g\n", segments[idmax](1), segments[idmax](2), segments[idmax](3)) ;
     r=Ico.belonging_tri(-segments[idmax]) ; 
     if (r!=-1) Ico.data[r]=Ico.data[r]+1 ; else printf("!") ;
    }
    K=K/Kn ; 
    Phi=sqrt((3*(K.norm())*(K.norm())-1)/2) ;
    //Calcul::eigen(K, Kval, Kvec) ;
    
    if (type==0)
    {if (Kvec(0,0)<0) {Kvec.col(0)=-Kvec.col(0) ;}}
    else if (type==1)
    {if (Kvec(2,0)<0) {Kvec.col(0)=-Kvec.col(0) ;}}
    else
      DISP_Err("Unknown multisphere type\n") ; 
    
    if (Kvec.determinant()<0) {Kvec.col(2)=-Kvec.col(2) ; }

    //fprintf(out, "%d %g %g %g %g %g %g %g %g %g %g %g %g %g\n", steps[i].timestep, Phi, Kval(0), Kval(1), Kval(2), Kvec(0,0), Kvec(0,1), Kvec(0,2), Kvec(1,0), Kvec(1,1), Kvec(1,2), Kvec(2,0), Kvec(2,1), Kvec(2,2)) ; 
    fprintf(out, "%d %g %g %g %g %g %g %g %g %g\n",steps[i].timestep, K(0,0), K(0,1), K(0,2),K(1,0), K(1,1), K(1,2),K(2,0), K(2,1), K(2,2)) ; 
    
  }
  
  if (actions["mean"].set)
  {
    chem=chemin ; 
    chem.append("-") ; chem.append("multisphere-2Dorient") ; chem.append(".txt") ; 
    out2=fopen(chem.c_str(), "w") ; 
    //theta=theta*180/M_PI ; 
    fprintf(out2, "# Informations on the average - loop parameters : %d:%d:%d\n", loop[0],loop[1],loop[2]) ; 
    double sumangular2=0 ; 
    for (i=0 ; i<72 ; i++) sumangular2+=angular2[i] ; 
    for (i=0 ; i<72 ; i++) fprintf(out2, "%g %g %g %d\n", i*theta, angular[i]/nangular/(theta), angular2[i]/sumangular2/theta, nangular) ; 
    fprintf(out2, "%g %g %g %d\n", 0*theta, angular[0]/nangular/(theta), angular2[0]/sumangular2/theta, nangular) ; 
    fclose(out2) ; 
    
    chem=chemin ; 
    chem.append("-") ; chem.append("multisphere") ; chem.append(".vtk") ; 
  
    int tot=0 ; 
    for (i=0 ; i<Ico.nfaces ; i++) tot+=Ico.data[i] ;
    for (i=0 ; i<Ico.nfaces ; i++) Ico.data[i]/=tot ;
    Ico.deform() ; 
    Ico.write_vtk(chem) ;
  }
  
}


// =====================================================
// Fonctions de la classe LcpDump ========================
//=====================================================
int LcpDump::read (unsigned char read_type, int index)
{
  /* read_type==1 : sauver les données lues
   * read_type==2 : ne lire qu'un dump
   * read_type==4 : creer la structure step et incrémenter nb_step */
// lecture d'un dump compressé de forces
char datas[500] ;  unsigned char UC ; unsigned short int USI ; unsigned int UI ; unsigned long int ULI ;
double D ; float F ; int I ; long int LI ;

Step null_step ;
int retour ; size_t res ; 
int i=0, j, k,l, idx ;
FILE *in ; in=dumpinc ;

Type=TL ; 

if (!checkfile()) {return -1 ;}

if (read_type & 4) {nbsteps=0 ; printf("Lecture initiale du dump compressé, peut prendre du temps ...\n") ; }
if (read_type & 2) {idx=index ;}

if (isfirst == true) // Lecture de l'header du fichier
   {
   res=fread(datas, 1, 26, in) ; datas[26]=0 ;  // Extraction de l'header AVFF
   if (!strcmp(datas,"AVFFLIGGGHTSDUMPCOMPRESSED")) cout << "Compression avec la version 1.0 du compresseur. Cela ne devrait pas poser de problème. Utiliser sinon une version plus ancienne de PostProcessing.\n" ;
   else if (!strcmp(datas,"AVFFLIGGGHTSDUMPCOMPRES2.0")) cout << "Compression avec la version 2.0, ne devrait pas poser de probl�me \n" ;
   else if (!strcmp(datas,"AVFFLIGGGHTSDUMPCOMPRES2.1")) cout << "Compression avec la version 2.1, pourrait causer des erreurs\n" ;
   else if (strcmp(datas,"AVFFLIGGGHTSDUMPCOMPRES2.3")) cout << "WARN : le fichier compressé ne semble pas avoir l'header correct\n" ;

   res=fread(&USI, 2, 1, in) ;
   res=fread(datas,1, USI, in) ; datas[USI]=0 ;
   strcpy(nomoriginal, datas) ;
   res=fread(&nbsteps, 4, 1, in) ;
   steps.resize(nbsteps) ;

   // Création des colonnes, extraction des formats et dénombrables
   res=fread(&steps[0].nb_idx, 1, 1, in) ;
   cp_dat.formats.resize(steps[0].nb_idx) ;
   cp_dat.same.resize(steps[0].nb_idx) ;
   cp_dat.champs.denom.resize(steps[0].nb_idx) ;
   cp_dat.champs.denombrables.resize(steps[0].nb_idx) ;
   for (i=0 ; i<nbsteps ; i++)
     {
     steps[i].multisphere=&multisphere ;   
     steps[i].idx_col.resize(steps[0].nb_idx) ;
     steps[i].nb_idx=steps[0].nb_idx ;
     }
   for (i=0 ; i<steps[0].nb_idx; i++)
     {
     res=fread(&steps[0].idx_col[i], 1, 1, in) ;
     cp_dat.champs.denombrables[i].resize(MAX_DENOM) ;
     res=fread(&UC, 1,1, in) ; cp_dat.formats[i]=UC&(~MASK_ALWAYS_THE_SAME) ;
     if ((UC & MASK_ALWAYS_THE_SAME) >> 7) cp_dat.same[i]=true ;
     else cp_dat.same[i]=false ;
     cp_dat.champs.denom[i]=0 ;
     if (cp_dat.formats[i]>=17) // C'est un dénombrable !
       {
       unsigned char nombre, format ;
       nombre=cp_dat.formats[i]>>4 ;
       switch (nombre)
         {
	 case 1 : nombre=1 ; break ;
	 case 2 : nombre=2 ; break ;
	 case 3 : nombre=4 ; break ;
	 case 4 : nombre=16 ; break ;
	 }
       cp_dat.champs.denom[i]=nombre ;
       format=cp_dat.formats[i]&15 ;
       for (j=0 ; j<nombre ; j++)
         {
	 if ((format==1) || (format==2))
	    {res=fread(&UC, 1, 1, in); cp_dat.champs.denombrables[i][j]=UC ;}     // TODO Séparer signed et unsigned
	 else if ((format ==3) || (format==4))
	    {res=fread(&USI, 2, 1, in) ; cp_dat.champs.denombrables[i][j]=USI ; }
	 else if ((format==5) || (format==6))
	    {res=fread(&UI, 4, 1,in)  ; cp_dat.champs.denombrables[i][j]=UI ; }
	 else if (format==9)
            {res=fread(&F,4,1,in) ; cp_dat.champs.denombrables[i][j]=F ; }
         else if ((format==7) || (format==8))
	    {res=fread(&LI,8,1,in) ; cp_dat.champs.denombrables[i][j]=LI ; }
	 else if (format==10)
	    {res=fread(&D,8,1,in) ; cp_dat.champs.denombrables[i][j]=D ; }
         else
            {cout << "AFormat inconnu dans la sous fonction convertDBL"  ; std::exit(EXIT_FAILURE) ;}
	 }
       }
     }
   for (i=1 ; i<nbsteps ; i++)
     {
     steps[i].idx_col=steps[0].idx_col ;
     if (steps[0].idx_col[0]>=64 && steps[0].idx_col[0]<128) {steps[i].Type=TCF ;}
     else if (steps[0].idx_col[0]>=128 && steps[0].idx_col[0]<255) {steps[i].Type=TF ;}
     else {steps[i].Type=TL ;}
     }
  steps[0].Type=steps[1].Type ;
  isfirst=false ;
  }

int beggining=0, endding ; endding=nbsteps ;
if (read_type & 2) {beggining=idx ; endding=idx+1 ;  fseek(in, steps[idx].posinfile,SEEK_SET) ;}
for (i=beggining ; i<endding ; i++)
    {
    if (read_type & 4) {steps[i].posinfile=ftell(in) ; }
    if (steps[i].Type==TL || steps[i].Type==TCF)
       {res=fread(&steps[i].timestep,4,1,in) ;
       res=fread(&steps[i].nb_atomes,4,1,in) ; }
    else if (steps[i].Type==TF)
       {res=fread(&steps[i].nb_pts,4,1,in) ;
       res=fread(&steps[i].nb_triangles,4,1,in) ; }

    res=fread(&F,4,1,in) ; steps[i].box[0][0]=F; res=fread(&F,4,1,in) ; steps[i].box[0][1]=F;
    res=fread(&F,4,1,in) ; steps[i].box[1][0]=F; res=fread(&F,4,1,in) ; steps[i].box[1][1]=F;
    res=fread(&F,4,1,in) ; steps[i].box[2][0]=F; res=fread(&F,4,1,in) ; steps[i].box[2][1]=F;
    steps[i].datas.resize(steps[i].nb_idx) ;

    for (j=0 ; j<steps[i].nb_idx ; j++)
        {
        int length_dat ;
        //printf("[%d]", steps[i].idx_col[j]) ;
        if (steps[i].idx_col[j]<128) length_dat=steps[i].nb_atomes ;
        else if (steps[i].idx_col[j]==IDS("POINTX") || steps[i].idx_col[j]==IDS("POINTY") || steps[i].idx_col[j]==IDS("POINTZ")) length_dat=steps[i].nb_pts ;
        else length_dat=steps[i].nb_triangles ;

	if (read_type & 1)   // Récupération et écriture des données
	  {
          steps[i].datas[j].resize(length_dat,0.0) ;
          if (cp_dat.same[j]==true && i>0)
             {
             steps[i].datas[j]=steps[0].datas[j] ;
             continue ;
             }

          if (cp_dat.formats[j]<17) // Ce n'est pas un dénombrable
           {
           unsigned char format ;
           format=(cp_dat.formats[j])&15 ;
           for (k=0 ; k<length_dat ; k++)
            {
            if ((format==1) || (format==2))
              {res=fread(&UC, 1, 1, in); steps[i].datas[j][k]=UC ;}// TODO Séparer signed et unsigned
            else if ((format==3) || (format==4))
              {res=fread(&USI, 2, 1, in) ; steps[i].datas[j][k]=USI ; }
             else if ((format==5) || (format==6))
              {res=fread(&UI, 4, 1,in)  ; steps[i].datas[j][k]=UI ; }
            else if (format==9)
              {res=fread(&F,4,1,in) ; steps[i].datas[j][k]=F ; }
            else if ((format==7) || (format==8))
              {res=fread(&LI,8,1,in) ; steps[i].datas[j][k]=LI ; }
            else if (format==10)
              {res=fread(&D,8,1,in) ; steps[i].datas[j][k]=D ; }
            else
              {cout << "BFormat inconnu dans la sous fonction convertDBL" ; std::exit(EXIT_FAILURE) ;}
	    }
	   }
	  else if (cp_dat.formats[j]>=cp_dat.CHAR_CST && cp_dat.formats[j]<cp_dat.CHAR_DENOM_2)
           {for (k=0 ; k<length_dat ; k++){steps[i].datas[j][k]=cp_dat.champs.denombrables[j][0] ;}}
          else if (cp_dat.formats[j]>=cp_dat.CHAR_DENOM_2 && cp_dat.formats[j]<=cp_dat.DOUBLE_DENOM_16)
           {
	   unsigned char bits, mask, octet ; int longueur ;
	   bits=(cp_dat.formats[j]>>4)&15 ;
	   switch(bits)
	      {
	      case 2 : bits=1 ; mask=1 ; break ;
              case 3 : bits=2 ; mask=3 ; break ;
              case 4 : bits=4 ; mask=7 ; break ;
	      }
	   longueur=ceil(length_dat*bits/8.) ;
	   for (k=0 ; k<longueur; k++)
	      {
	      res=fread(&octet, 1,1, in) ; octet=cp_dat.binary_swap(octet, bits) ;
		for (l=0 ; l<(8/bits) ; l++)
	          {steps[i].datas[j][k*(8/bits)+l]=cp_dat.champs.denombrables[j][octet&mask] ; octet=octet>>bits ;  }
              }
	   }
          else
           {cout << "ERR : format inconnu lors de la décompression des steps." ; }
          }

	else 	// On ne recherche que les positions dans le fichier qui pointent vers les données de chaque ts.
	  {
          if (cp_dat.same[j]==true && i>0)
             {continue ; }
	  if (cp_dat.formats[j]<17) // Ce n'est pas un dénombrable
           {
           unsigned char longueur ; unsigned char format ; format=(cp_dat.formats[j])&15 ;

            if ((format==1) || (format==2)) longueur=1 ;
            else if ((format==3) || (format==4)) longueur=2 ;
            else if ((format==5) || (format==6) || (format==9)) longueur=4 ;
            else if ((format==7) || (format==8) || (format==10)) longueur=8 ;
            else {cout << "Format inconnu dans la sous fonction convertDBL" ; std::exit(EXIT_FAILURE) ;}

	   fseek(in, length_dat*longueur, SEEK_CUR) ;
	   }
	  else if (cp_dat.formats[j]>=cp_dat.CHAR_CST && cp_dat.formats[j]<cp_dat.CHAR_DENOM_2) {/*do nothing*/}
          else if (cp_dat.formats[j]>=cp_dat.CHAR_DENOM_2 && cp_dat.formats[j]<=cp_dat.DOUBLE_DENOM_16)
           {
	   unsigned char bits, mask, octet ; int longueur ; bits=(cp_dat.formats[j]>>4)&15 ;
	   switch(bits)
	      {case 2 : bits=1 ; mask=1 ; break ;
               case 3 : bits=2 ; mask=3 ; break ;
               case 4 : bits=4 ; mask=7 ; break ; }
	   longueur=ceil(length_dat*bits/8.) ;
	   fseek(in, longueur, SEEK_CUR) ;
	   }
          else
           {cout << "ERR : format inconnu lors de la décompression des steps." ; }
	  }
	}
    }
rewind(in) ;
return 0 ;
}
//-------------------------------------------------------
bool LcpDump::checkfile (void)
{
if (dumpinc==NULL)
  {
  cout << "ERR ! Le fichier de dump doit être ouvert au format C (FILE*) pour être utilisé dans un dump compressé. Des erreurs surviendront." ;
  return false ;
  }
return true ;
}
//---------------------------------------------------------
int LcpDump::uncompress()
{
char chemintempor[500]; string chemintmp ;
sprintf(chemintempor,"%s.uc",nomoriginal) ;
chemintmp=chemintempor ;
if (steps[0].Type==TL)
   {cout << "Décompression d'un dump atomique\n" ; write_asDUMP(chemintmp) ; }
else if (steps[0].Type==TCF)
   {cout << "Décompression d'un dump de chainforce\n" ; write_asDUMP(chemintmp) ; }
else
   {cout << "Décompression d'un dump de stress\n" ; write_asOneVTK(chemintmp) ; }
return 0 ;
}
//-------------------------------------------------------
int LcpDump::free_timestep (long int inscrit)
{
int i ;
bool exist_same = false ;
// Do not free s'il existe un champ same
for (i=0 ; i<steps[inscrit].nb_idx ; i++)
    {if (cp_dat.same[i]==true) exist_same=true ; }
if (exist_same && inscrit==0)
   {inscrit=-1 ; return inscrit ; }
if (steps[inscrit].nb_atomes>0)
  {
  for (i=0 ; i<steps[inscrit].nb_idx ; i++)
  {steps[inscrit].datas[i].clear() ; }
  }
/* Créé des erreurs de segmentation dans certains cas pour une raison que je ne m'explique pas. On commente (c'est pas pour la place que ça prend ...)  mais c'est bizarre
steps[inscrit].datas.clear() ; */
steps[inscrit].filtered=false ; 
inscrit=-1 ;
return inscrit ;
}
//----------------------------------------------------
void LcpDump::disp(void) {
if (steps[0].Type==TL || steps[0].Type==TCF)
   {LDump::disp() ; fflush(stdout) ; }
else if (steps[0].Type==TF)
   {FDump::disp() ; fflush(stdout) ; }
else
  cout<<"WARN : le type de steps n'est pas connu pour les dump compressé, disp impossible" ;
}
