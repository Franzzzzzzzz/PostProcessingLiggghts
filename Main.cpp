#include "Headers/Main.h"


Actions actions ;
IDSCLASS IDS ; 

int main (int argc, char *argv[])
{
int i, res ;
   
/*//BEGIN TEST

Matrix3d A ; A=Matrix3d::Random() ; 

Matrix3d Vec ;
Vector3d Lambda;
Calcul::eigen(A, Lambda, Vec) ; 


return 0 ; 
//END TEST*/

// Création des arguments possibles
actions.initialize() ; 
// lecture des arguments de ligne de commande
actions.parse_args(argc, argv) ; 
//actions.parse_args_fromfile() ; 
actions.write_history(0, argc, argv) ;

// Traitement des arguments requierant une action immédiate particulière
if (actions["v1"].set) actions["noperiodicity"].set=true ;
if (actions["mkdir"].set) actions.ecrire.mkdir=true ; 
if (actions["clean"].set) {actions.clean(argv[argc-1]) ; }
if (actions["help"].set) {actions.aide() ; std::exit(EXIT_SUCCESS) ; }
if (actions["version"].set) {cout << "Date de compilation : " << __DATE__ << "\n" ; std::exit(EXIT_SUCCESS) ; }
if (actions["restart2vtk"].set) {actions.regenerate_restart2vtk_script() ; std::exit(EXIT_SUCCESS) ; }
actions.Cst.set_constantes () ;



if (actions["vtk"].set)
     {
     FDump fdump ; 
     fdump.open(actions.dumpnames[0]) ;

     fdump.disp() ;

     if (actions["separe"].set)
        fdump.write_asVTK(actions.dumpnames[0]) ;
     if (actions["recode"].set)
        fdump.write_asOneVTK(actions.dumpnames[0]) ;
     if (actions["w/forcetot"].set)
        fdump.write_forcestot(actions.dumpnames[0]) ;
     if (actions["w/coupletot"].set)
        fdump.write_coupletot(actions.dumpnames[0]) ;
     if (actions["mean"].set)
        fdump.mean_stress(actions.dumpnames[0]) ; 
     if (actions["compress"].set)
        {
	FILE * out ; 
        string temp ; 
        temp=actions.dumpnames[0]+".cp" ; 
        out=fopen(temp.c_str(), "wb") ; 
        Compresser tmp ; 
        tmp.compress (fdump, actions.dumpnames[0], out) ;
        fclose(out) ; 
        }
     }

else if (actions["chainforce"].set)
     {
     LcfDump cfdump ; 
     LucDump dump, walldump ;
     Filter filtre_tmp ; 

     char filterdo[50] ;

     if (!actions["compress"].set) 
        {
     	// Création du filtre pour le Ldump de positions atomiques : rangement des atomes par id
     	if (actions["iscylperiodic"].set) sprintf(filterdo, "id::sort::null;null::quartx+y+::null") ;
     	else sprintf(filterdo, "id::sort::null") ;
     	dump.prefiltre=filtre_tmp.parse_arg(filterdo,TL) ;

     	// Création des filtres pour le CFdump : pas de chaines à travers les pbc, filtre de création des CFR CFMAG etc.
     	if (!actions["noperiodicity"].set) sprintf(filterdo, "null::.cfpp.::null;null::.nopbc.::null") ;
     	else sprintf(filterdo, "null::.cfpp.::null") ;
     	cfdump.prefiltre=filtre_tmp.parse_arg(filterdo,TCF) ;

     	// Création du lien avec les autres dump dans le CFdump, et ajout d'un postfiltre si nécessaire
     	cfdump.prefiltre[0].alter_dump=&(dump) ; cfdump.prefiltre[1].alter_dump=&(dump) ;
     	cfdump.prefiltre[0].wall_dump=&(walldump) ; cfdump.prefiltre[1].wall_dump=&(walldump) ;
     	if (actions["iscylperiodic"].set)
           { sprintf(filterdo, "null::quartx+y+::null;null::.cpplink.::null") ;
    	cfdump.postfiltre=filtre_tmp.parse_arg(filterdo,TCF) ; }

     	// Ouverture réelle des dumps annexes
     	dump.open(actions.dumpnames[1]) ;
     	if (actions["wallchainforce"].set)
	   walldump.open(actions.dumpnames[2]) ;
        }

     cfdump.open(actions.dumpnames[0]) ; 
     
     //dump.disp () ; 
     //cfdump.disp() ; 

     if (actions["coarse-graining"].set)
        {
        CoarseDump dcor ; 
        if (dump.nbsteps != cfdump.nbsteps)
	   DISP_Warn("ATTENTION : dump atomique et de chainforce devraient être au même échantillonnage (ne semble pas le cas pour l'instant ...) (effectuer d'abord un /PostProcessing --reechantillonner --downsampling xx --extract xx yy dump.atm/ si besoin)\n") ; 

        cfdump.check_timestep(0) ;

        /*for (i=0 ; i<dump.nbsteps ; i++) // NOT NEEDED I THINK SINCE .cfpp. filter does it, but to check.
          {
          dump.steps[i].Type=TL ; if (actions["wallchainforce"].set) {walldump.steps[i].Type=TL ;}
          cfdump.steps[i].box[0][0]=dump.steps[0].box[0][0] ; cfdump.steps[i].box[0][1]=dump.steps[0].box[0][1] ; cfdump.steps[i].box[1][0]=dump.steps[0].box[1][0] ;
          cfdump.steps[i].box[1][1]=dump.steps[0].box[1][1] ; cfdump.steps[i].box[2][0]=dump.steps[0].box[2][0] ; cfdump.steps[i].box[2][1]=dump.steps[0].box[2][1] ;
          }*/
        //dcor.do_coarse(dump) ;
        if (actions["noatmcoarse"].set) {LucDump nulldump ; dcor.do_coarse(cfdump, nulldump, 0) ;}
        else dcor.do_coarse(cfdump, dump, 2) ;
        //dcor.write(chemin);
        //dcor.write_asMatlabSpatioTemp (actions.dumpnames[0]) ;
	
        if (actions["mean"].set)
          {
          dcor.mean() ;
          dcor.write(actions.dumpnames[0]) ;
          }
        else
          {
          dcor.write(actions.dumpnames[0]) ;
          }
        }
     else if (actions["grainforce-by-angle"].set || actions["grainforce"].set)
    	cfdump.write_grainforce(actions.dumpnames[0], dump) ;
     //else if (actions["grainforce"].set)
    	//cfdump.write_grainforcetot(actions.dumpnames[0]) ;
     else if (actions["grain-rayon-around"].set)
    	cfdump.write_radiuscontact() ;
     else if (actions["forcetank"].set)
        cfdump.forcetank(actions.dumpnames[0], dump) ; 
     else if (actions["w/coupletot"].set)
     {
	cfdump.coupletot(actions.dumpnames[0], dump) ; 
	//cfdump.write_asVTK (actions.dumpnames[0], dump);
     }
     else if (actions["compress"].set)
     	 {
    	 FILE * out ;
    	 string temp ;
    	 temp=actions.dumpnames[0]+".cp" ;
    	 out=fopen(temp.c_str(), "wb") ;
    	 Compresser tmp ;
    	 tmp.compress (cfdump, actions.dumpnames[0], out) ;
     	 }
     else if (actions["justtmp"].set) {} // do nothing	 
     else
        cfdump.write_asVTK (actions.dumpnames[0], dump);
     }

else 									//Atom dump (chainforce not set)
     {
     
     if (actions["cp"].set)   // C'est un dump compressé // TODO remove that ... gunzip est bien mieux ...
       {
       LcpDump dump ;
       // if (filtre_tmp.lst_op.size()>0) {dump.filtre=filtre_tmp.lst_op ; } 
       dump.dumpinc=fopen(actions.dumpnames[0].c_str(), "rb") ;
       res=dump.open(actions.dumpnames[0]) ;
       dump.disp() ; 
       if (actions["uncompress"].set)
        {
        dump.uncompress() ; 
        std::exit(EXIT_SUCCESS) ; 
	    }
       if (actions["dump2vtk"].set)
        { dump.write_asVTK(actions.dumpnames[0]) ; }
       if (actions["coarse-graining"].set)
        {
	CoarseDump dcor ; LcfDump nulldump ; 
        dcor.do_coarse(dump, nulldump, 1) ; 
       // dcor.write(chemin) ; 
        if (actions["mean"].set)
           {
	   dcor.mean() ;
           dcor.write(actions.dumpnames[0]) ;
           }
         //dcor.write_asMatlab(chemin) ;
         } 
       if (actions["surface"].set)
        {
        Surface surf ; 
        surf.detect_surface (dump) ; 
        surf.write_asmatlab (actions.dumpnames[0])  ; 
        }
       if (actions["justtmp"].set) {} // Do nothing 
       }

     else			// C'est un dump décompressé
       {
       LucDump dump; LcfDump walldump ;
       
       if (actions["multisphere"].set)
       {
     	 Filter filtre_tmp ; char filterdo[50] ;
	 sprintf(filterdo, "id::sort::null;id::fill::null") ;
       	 dump.prefiltre=filtre_tmp.parse_arg(filterdo,TL) ;
       }
     
       if (actions["wallchainforce"].set) // Wallchainforce sans chainforce. If faut filtrer le dump de test avec le wallchainforce
       {
         Filter filtre_tmp ; char filterdo[50] ;
         sprintf(filterdo, "null::.wallforceatm.::null") ;
         dump.prefiltre=filtre_tmp.parse_arg(filterdo,TL) ;
         dump.prefiltre[0].wall_dump=&(walldump) ;
         walldump.open(actions.dumpnames[1]) ;
       }
     
       res=dump.open(actions.dumpnames[0]) ;
       if (res==-1) {std::exit(EXIT_FAILURE) ; }
       //dump.disp();
       //Stats stat ; stat.compute_step(dump, dump.nbsteps-10) ; stat.disp() ;  
       if (actions["dstminmax"].set) 
       	   {
    	   Stats stat ;
    	   stat.minmaxdst(dump, (int) actions["dstminmax"]["timestep"]) ;
       	   }
       if (actions["compress"].set)
        {
    	FILE * out ;
        string temp ; 
        temp=actions.dumpnames[0]+".cp" ; 
        out=fopen(temp.c_str(), "wb") ; 
        Compresser tmp ; 
        tmp.compress (dump, actions.dumpnames[0], out) ;
        }
       if (actions["dump2vtk"].set)
        {
	  if (actions["multisphere"].set) 
	  {
	    Filter filtre_tmp ; char filterdo[50] ;
	    sprintf(filterdo, "id::sort::null;id::fill::null") ;
	    dump.prefiltre=filtre_tmp.parse_arg(filterdo,TL) ;
	  }
	  dump.write_asVTK(actions.dumpnames[0]) ;
	}
       if (actions["dump2restart"].set)
	    {dump.write_asRESTART(actions.dumpnames[0]) ; }
       if (actions["coarse-graining"].set)
        {
	LcfDump nulldump ; 
	CoarseDump dcor ;
        dcor.do_coarse(nulldump, dump, 1) ; 
        if (actions["mean"].set)
          { 
	  dcor.mean() ;
          dcor.write(actions.dumpnames[0]) ;
          }
        else
          dcor.write(actions.dumpnames[0]) ; 
         //dcor.write_asMatlab(chemin) ; 
        }
       if (actions["coarse-graining-basic"].set)
        {
	CoarseDump dcor ;
        dcor.do_coarse_basic(dump, actions.dumpnames[0]) ;
        }
       if (actions["surface"].set || actions["surfaces"].set)
        {
        Surface surf ;
        if (actions["surfaces"].set) {actions.copyarg("surfaces", "surface") ; }
        surf.detect_surface (dump) ; 
        surf.write_asmatlab (actions.dumpnames[0])  ; 
        if (actions["surfaces"].set) surf.write_asmatlab_bottom(actions.dumpnames[0]) ;
        }
       if (actions["surface2D"].set)
               {
               Surface surf ;
               surf.detect_surface2D (dump) ;
               surf.write_asmatlab2D (actions.dumpnames[0])  ;
               }
       //if (actions["downsampling"].set)
       //  dump.write_asDUMP(actions.dumpnames[0]) ;
       if (actions["w/forcetot"].set)
	 dump.write_forcestot(actions.dumpnames[0]) ; 
       if (actions["wallforce-by-angle"].set)
         dump.write_wallforce(actions.dumpnames[0]) ; 
       if (actions["xray"].set)
	 dump.write_xray(actions.dumpnames[0]) ; 
       if (actions["multisphere"].set && !actions["dump2vtk"].set && !actions["coarse-graining"].set)
       {     
	 Filter filtre_tmp ; char filterdo[50] ;
	 sprintf(filterdo, "id::sort::null;id::fill::null") ;
     	 dump.prefiltre=filtre_tmp.parse_arg(filterdo,TL) ;
	 dump.write_multisphere_dumbell(actions.dumpnames[0]) ; 
       }
       if (actions["justtmp"].set) {} //do nothing
       }
     }

actions.write_history(1, argc, argv) ;
printf("\n") ; // Plus joli à la fin...
std::exit(EXIT_SUCCESS) ;
}

//----------------------------
Actions::Actions() : numdump(1), using_loop(false), Cst(*this)
{
  Cst.add("Radius", 0.00075) ; 
  Cst.add("Rhograin", 2500) ; 
  Cst.add("G", 9.81) ; 
  Cst.add("D", 0.005) ; 
  Cst.add("R", 0.0025) ; 
  Cst.add("L", 0.03) ; 
  Cst.add("LSimu", 0.05) ; 
  Cst.add("DSimu", 0.015) ; 
}

double Actions::Constantes::operator[] (string n) 
{ 
 for (int i =0 ; i<count ; i++) 
 {
   if (nom[i]==n) 
   {
     if (!info[i])
     { char message[500] ; 
       sprintf(message, "Utilisation d'une constante (1e fois) : %s=%f.\n", nom[i].c_str(), valeur[i]) ;
       DISP_Info(message) ; info[i]=true ;
     }
     return (valeur[i]) ;
   }
 }
DISP_Err("Actions : Constante inconnue") ; printf("%s\n",n.c_str()) ;   throw 1 ; 
}
//---------------------------
void Actions::initialize(void) 
{
int dep[10] ; 
char * args [10] ; 

actions.new_arg ("writing", "Paramètres pour l'écriture des fichiers de sorti (cf. la documentation)", 0, 0) ; 
actions.new_arg ("W", "Raccourci du précédent", 0, 0) ; 

actions.new_arg ("alias", "génère / appelle / liste / décrit un alias, selon l'argument donné sous la forme : set:alias, call:alias,list, info:alias", 0, 0) ; 

actions.new_arg ("vtk", "le dump est au format vtk", 0, 0) ; 
dep[0]=actions["vtk"].id ; 
actions.new_arg ("separe", "coupe un vtk à chaque pas de temps", 0, 1, dep) ;
actions.new_arg ("w/forcetot", "écrit un fichier dump.forces.txt contenant la force totale à chaque pas de temps", 0, 0) ;

actions.new_arg ("w/coupletot", "écrit un fichier dump.forces.txt contenant le couple total à chaque pas de temps. \
Avec vtk : a partir des données dans un vtk. \
Avec Chainforce : en comparant les écart de force propre des grains dans le ldump avec la somme des chaines de forces par grains. ", 0, 0) ;

actions.new_arg ("w/otige", "Retire la tige de support dans le calcul du couple pour le cas où on l'a simulée", 0, 1, dep) ;
actions.new_arg ("group-triangles" , "Fait la moyenne des triangles de même normale dans un mesh", 0, 1, dep) ;
args[0]=(char *) "anglebystep" ; 
actions.new_arg ("anglebystep", "si le mesh est en rotation autour du centre indiqué, angle de rotation effectué en 1 pas de step (en degrés)", 1, args, 0) ;
actions.new_arg ("recode", "recode le vtk en enlevant des signes - en vadrouille pour comparer avant/après compression", 0, 1, dep) ;
actions.new_arg ("norebuild", "Ne tente pas de reconstruire des données manquantes. Plutôt utilisé pour les vtk dans la cadre de Liggghts 1 (--v1)", 0, 0) ;

actions.new_arg ("cp", "identifie un dump compressé", 0, 0) ;
dep[0]=actions["cp"].id ;  
actions.new_arg ("uncompress", "décompresse un dump compressé", 0, 1, dep) ;

args[0]=(char *) "nbbox_x" ; args[1]=(char *) "nbbox_y" ; args[2]=(char *) "nbbox_z" ; 
actions.new_arg ("coarse-graining", "extrait des vtk avec coarse graining en vitesse d'un dump. Les params sont le nb de boites en x, y, z.",3, args, 0) ;
actions.new_arg ("coarse-graining-basic", "Extrait les vitesses coarsed moyennes (fonctionne pour type non granulaires aussi). Les params sont le nb de boites en x, y, z.",3, args, 0) ;
dep[0]=actions["coarse-graining"].id ; args[0]=(char *) "soustraire" ; 
actions.new_arg ("substract", "soustrait la valeur donnée à la valeur x du champ de vitesse", 1,args, 1, dep) ;
args[0]=(char *) "xcyl" ; args[1]=(char *) "ycyl" ; args[2]=(char *) "zcyl" ; 
args[3]=(char *) "dirxcyl" ; args[4]=(char *) "dirycyl" ; args[5]=(char *) "dirzcyl" ; 
args[6]=(char *) "rayon" ; args[7]=(char *) "longueur" ; 
actions.new_arg ("delcyl", "Supprime le cylindre du calcul de coarse graining dans le cas des winspheres, pour un meilleur PHI près du cylindre en particulier.", 8, args, 1, dep) ;
args[0]=(char *) "rayon" ;
actions.new_arg ("deltank", "Supprime le tank du calcul de coarse graining dans le cas des chaines de forces, pour un meilleur P près du tank en particulier.", 1, args, 1, dep) ;
actions.new_arg ("set-vy-zero", "annule la composante y de la vitesse (pour tracer les lignes de courant ...)", 0, 1, dep) ;
args[0]=(char *) "box_xmin" ; args[1]=(char *) "box_xmax" ; args[2]=(char *) "box_ymin" ; args[3]=(char *) "box_ymax" ; args[4]=(char *) "box_zmin" ; args[5]=(char *) "box_zmax" ;
actions.new_arg ("use-box", "ne pas utiliser la boîte de la simu mais une boîte perso", 6,args, 0) ;
args[0]=(char *) "valeur";
actions.new_arg ("setdeltaboite", "indique la valeur de deltaboite", 1, args, 0) ; 
args[0]=(char *) "sig-x" ; args[1]=(char *) "sig-y" ; args[2]=(char *) "sig-z" ;
actions.new_arg ("sig-boites", "taille des boites de coarse graining en x y z",3,args, 1, dep) ;
actions.new_arg ("winboxyper", "utiliser un fenêtrage box precis grace a une périodicité en y", 0, 1, dep) ;
args[0]=(char *) "sigma" ; 
actions.new_arg ("wingauss", "utiliser un fenêtrage gaussien (sigma=-1 : le calculer automatiquement)", 1,args, 1, dep) ;
args[0]=(char *) "sigmax" ; args[1]=(char *) "sigmay" ; args[2]=(char *) "sigmaz" ; 
actions.new_arg ("wingauss3D", "utiliser un fenêtrage gaussien avec 3 sigma differents (sigmax|y|z=-1 : le calculer automatiquement)", 3,args, 1, dep) ;
actions.new_arg("w/kinetic", "calcul aussi la partie cinétique du tenseur des contraintes", 0, 1, dep) ;
actions.new_arg("pressureFgrain", "Calculer un vecteur de pression à partir de la force sur chaque grain", 0, 1, dep) ;
args[0]=(char *) "nbbox_x" ; args[1]=(char *) "nbbox_y" ;
actions.new_arg ("surface", "trouve la surface libre haute avec l'échantillonnage donné en x et y", 2, args, 0) ;
actions.new_arg ("surface2D", "trouve la surface libre haute en yavec l'échantillonnage donné en x ", 1, args, 0) ;
actions.new_arg ("surfaces", "trouve les surfaces libres haute & basse avec l'échantillonnage donné en x et y", 2, args, 0) ;

args[0]=(char *) "dir" ; args[1]=(char *) "width" ; args[2]=(char *) "height" ; 
actions.new_arg ("xray" , "effectue une projection rayons y (image de la densité locale) selon une direction donnée. (direction:x=0, y=1, z=2) ", 3, args, 0) ; 

actions.new_arg ("dumpall", "charger un dump en mémoire en entier", 0, 0) ;
args[0]=(char *) "downsampling" ; 
actions.new_arg ("downsampling", "n'extraire qu'un timestep sur n", 1 ,args, 0) ;
args[0]=(char *) "timestep" ;
actions.new_arg ("dstminmax", "Calculer la distance minimale et maximale entre atome au timestep ts", 1, args, 0) ; 
args[0]=(char *) "extract_deb" ; args[1]=(char *) "extract_fin" ;
actions.new_arg ("extract", "n'extraire que les timestep entre arg1 et arg2", 2 ,args, 0) ;
actions.new_arg ("mean", "ne pas faire de dump par pas de tps mais seulment des moyenne", 0, 0) ;
args[0]=(char *) "xcyl" ; args[1]=(char *) "zcyl" ; args[2]=(char *) "nbbox_theta" ; args[3]=(char *) "sigma" ;
actions.new_arg ("wallforce-by-angle", "Pour des interactions ycylindre-grains,calcul les forces en fonction de l'angle. arg1=posx, arg2=posz, arg3=discretisation en angle, arg4=sigma en degré pour moyennage gaussien (si =0 moyenne créneau) ", 4, args, 0) ;

args[0]=(char *) "timestep" ; 
actions.new_arg ("dump2restart", "Ecrit un fichier restart pour charger des atomes dans liggghts", 1, args, 0) ; 
actions.new_arg ("dump2vtk", "transforme un dump en vtk séparés", 0, 0) ;
dep[0]=actions["dump2vtk"].id ;
actions.new_arg ("w/speed", "ajoute les données de vitesse au dump2vtk", 0, 1, dep) ;
actions.new_arg ("w/force", "ajoute les forces au dump2vtk", 0, 1, dep) ;
actions.new_arg ("w/id", "ajoute les numéros d'identification au dump2vtk", 0, 1, dep) ;
actions.new_arg ("w/rayon", "ajoute les rayons au dump2vtk", 0, 1, dep) ;
actions.new_arg ("w/masse", "ajoute les masses au dump2vtk", 0, 1, dep) ;
args[0]=(char *) "type" ; 
actions.new_arg ("multisphere", "les particules sont en fait des multispheres... Type: 0=long, 1=flat", 1,args,0) ; // dep va surement sauter bientôt, il n'y a pas de raison d'en rester au dump2vtk...
args[0]=(char *) "axes" ; 
actions.new_arg ("symetriser", "les particules seront symetrisées selon les directions choisies. Ex: 100: symetrie axe x, 101: symetrie x & z", 1, args, 0) ; 
actions.new_arg ("eigen", "Prend les eigenvalues / eigenvectors des tenseurs. Le writing du tenseur génère alors les eigenvectors (columns), et demander lambda génère les eigenvalues",0,0) ; 
actions.new_arg ("compress", "compresse un dump ou un vtk ou un chainforce dump (dans ce cas, mettre un null en avant-dernier argument !!!)", 0, 0) ; 
actions.new_arg ("reechantillonner", "Réécris le dump atomique en modifiant si besoin l'échantillonnage et les timestep de départ et de fin.", 0,0)  ;
// Pour filter : cout << "\t\t\t Ex : \"x<3.25;vx<vz;id::sort::null\n\n" ;

actions.new_arg ("chainforce", "s'occupe des chaines de force. ATTENTION : l'avant dernier argument doit être le dump des positions, le dernier le dump des chaines de forces.", 0, 0) ;
dep[0]=actions["chainforce"].id ;
actions.new_arg("donotusecfpos", "ne pas utiliser le cfdump pour reconstruire les positions des particules",0,1,dep) ; 
args[0]=(char *) "slices" ;
actions.new_arg("forcetank", "force exercée par le tank sur une tranche de grains", 1, args, 1, dep) ; 
actions.new_arg ("filter", "filtre les données de LDUMP avec l'argument suivant entre guillemets. Chaque opération est séparé par un point-virgule (cf. examples ailleurs TODO)", 0, 0) ;
actions.new_arg ("F", "Raccourci du précédent",0,0) ; 
actions.new_arg ("filterCF", "filtre les données de CFDUMP avec l'argument suivant entre guillemets. Chaque opération est séparé par un point-virgule (cf. examples ailleurs TODO)", 0, 1, dep) ;
actions.new_arg ("FCF", "Raccourci du précédent",0,1, dep) ; 
args[0]=(char *) "xcyl" ; args[1]=(char *) "zcyl" ; args[2]=(char *) "rayon" ;
actions.new_arg ("wallchainforce", "utilise un dump externe pour les liaisons grains-murs. 2 manières de l'utiliser:\n 1) avec --chainforce, 3 dump: wallforce, atoms, chainforces.\n 2/ sans --chainforce, 2 dumps: wallforce suivi de test.", 3, args , 0, dep) ;
actions.new_arg ("wallchainforcenodump", "calcul les liaisons grains-murs par équilibre des forces", 0, 1, dep) ; 
args[0]=(char *) "cutoff" ;
actions.new_arg("clone-periodic-chain", "copie les chaines de force traversant les parois periodiques (pas pour cylperiodic !)", 0, 1, dep) ;
actions.new_arg("cutoff", "ne garde que les n percent de chaînes de force de plus grande magnitude", 1, args, 1, dep) ;
actions.new_arg("noperiodicity", "les chaines de forces n'indiquent pas la periodicité (LIGGHTS version < 1.4.6)", 0, 1, dep) ;
actions.new_arg("v1", "modifie les parametres de PostProcessing pour la compatibilité avec Liggghts v1.x.x. -> autoset --noperiodicity", 0,0) ;
args[0]=(char *) "id" ; args[1]=(char *) "nbbox_theta" ; args[2]=(char *) "sigma" ;
actions.new_arg("grainforce-by-angle", "répartition angulaire des forces sur un grain", 3, args, 1, dep) ;
actions.new_arg("grainforce", "répartition angulaire des forces sur un grain", 1, args, 1, dep) ;
actions.new_arg("grain-rayon-around", "répartition des rayons autour d'un grain", 1, args, 1, dep) ;
dep[0]=actions["coarse-graining"].id ; 
actions.new_arg("is2D", "indique qu'il s'agit d'un coarse 2D", 0, 0) ;
actions.new_arg("break", "Supposer qu'il y a des données de break", 0,1,dep) ; 
dep[1]=actions["chainforce"].id ; 
actions.new_arg("noatmcoarse", "ne pas effectuer le coarse atomique", 0, 2, dep) ;
actions.new_arg("wincreneau", "Moyenner sur une fenetre creneau", 0, 1, dep) ;
actions.new_arg("iscylperiodic", "Dans le cas où on ne simule qu'1/4 de la boîte (x+, y+, z) ; utilisé en conjonction avec w/coupletot et coarse-graining",0,0) ;
actions.new_arg("mkdir", "Créer un nouveau dossier pour stocker les fichiers ecrits. Si le dossier existe, il est utilisé", 0, 0) ;
actions.new_arg("justtmp", "Créer seulement les fichiers tmp pour les dumps, sans faire aucun calcul", 0, 0) ;
actions.new_arg("forcetmp", "Force l'écriture du tmp indépendamment de se taille relative avec le dump",0,0) ; 
actions.new_arg("nofileerror", "Rend silencieuses les erreurs de lecture de fichier", 0, 0);
actions.new_arg("clean", "Nettoie des fichiers vtk et tmp avant tout calcul.", 0, 0);
actions.new_arg("restart2vtk", "Genere le script restart2vtk pour transformer les restart en vtk", 0, 0);
actions.new_arg("help", "Affiche tous les arguments possibles et sort.", 0, 0);
actions.new_arg("bashcomplete", "Genère le script pour la completion bash", 0, 0);
actions.new_arg("version", "Indique la date de compilation.", 0, 0);
actions.new_arg("teststress" , "Génére des forces par triangle idéales pour tester",0,0) ; 
}

//---------------------------
void Actions::set_progress(long int *lp) 
{
  loop[0]=lp[0] ; loop[1]=lp[1] ; loop[2]=lp[2] ; valeur=loop[0] ; total=(loop[2]-loop[0])/loop[1] ; using_loop=true ; return ; 
}

void Actions::disp_progress(void)
{
static pthread_t progress ;
static bool alreadyone=false ; 
cout << "\n\n" ; fflush(stdout) ;

if (alreadyone)
{
 pthread_cancel (progress) ; 
 pthread_join (progress, NULL) ;  
}

alreadyone=true ; 
pthread_create (&progress, NULL, progression, NULL);
}

void * progression (void * rien)
{
  double v ; 
v=0 ;   
while (v<99)
  {
    if (actions.using_loop)
     {
     v=(actions.valeur-actions.loop[0])/((double)(actions.loop[2]-actions.loop[0]))*100. ; 
     if (actions.valeur > actions.loop[2]+actions.loop[1] ) {DISP_Warn ("Boucle de progression incohérente") ;}
     if (actions.valeur < actions.loop[0]-actions.loop[1] ) {printf("%g %d %d %d\n", actions.valeur, actions.loop[0], actions.loop[1], actions.loop[2]); DISP_Warn ("Boucle de progression incohérente") ; v=0 ;}
     } 
    else
     v=actions.valeur/actions.total*100. ;  

  cout << setw(5) << setiosflags(ios::fixed) <<setprecision(1) << v << "%\r" ; fflush(stdout) ; 
  sleep(1) ; 
  }
return NULL ;
}

void Actions::clean(string chemin)
{
char reponse ; int res ; 
char commande[1000] ; 
size_t pos ; 

cout << "!! Attention : les fichiers .vtk et .tmp du répertoire où se trouve le dump seront supprimés avant toute opération.\nÊtes-vous sûr de vouloir cela ([Yy | Nn]) ? \n" ;
cin  >> reponse ; 

if (reponse=='Y' || reponse=='y')
  {
  pos=chemin.find_last_of("/\\");
  if (pos==string::npos)
     {
     res=system("rm *.tmp ; rm *.vtk ;") ; 
     }
  else
     {
     chemin=chemin.substr(0,pos+1) ; 
     sprintf(commande,"rm %s*.vtk ; rm %s*.tmp", chemin.c_str(),chemin.c_str()) ;
     res=system(commande) ; 
     }
  }

}

//===========================================================
int Actions::new_arg (string commande, string description, int nb_args, int nb_dep)
{ if (nb_args==0 && nb_dep==0)
     return (new_arg (commande, description, 0, (char **) NULL, 0, (int *) NULL)) ;
  else
     cout << "ERR : il manque des arguments au nouvel argument " << commande ; 
return 1 ;
}
int Actions::new_arg (string commande, string description, int nb_args, char * noms_args[], int nb_dep)
{
if (nb_dep==0)
   return (new_arg(commande, description, nb_args, noms_args, 0, (int *)NULL)) ;
else
   cout << "ERR : il manque des arguments au nouvel argument " << commande ;
return 1 ;
}
int Actions::new_arg (string commande, string description, int nb_args, int nb_dep, int id_dep[])
{
if (nb_args==0)
   return (new_arg (commande, description, 0, (char **) NULL, nb_dep, id_dep)) ;
else
   cout << "ERR : il manque des arguments au nouvel argument " << commande ;
return 1 ;
}
//----------------------------------------------------------
int Actions::new_arg (string commande, string description, int nb_args, char * noms_args[], int nb_dep, int id_dep[])
{
int i ; 
Arguments tmp ; 

tmp.set=false ; 
tmp.commande=commande ; 
tmp.description=description ;
if (nb_args>0)
   {
   tmp.params.resize(nb_args, 0.0) ; 
   for (i=0 ; i<nb_args ; i++)
    {
    tmp.p_names.push_back(noms_args[i]) ; 
    }
   }
if (nb_dep>0)
   {
   for (i=0 ; i<nb_dep ; i++)
       {tmp.dependances.push_back(id_dep[i]) ;}
   }
actions.arguments.push_back(tmp) ; 
actions.arguments[actions.arguments.size()-1].id=actions.arguments.size() ; 
return (actions.arguments[actions.arguments.size()-1].id) ; 
}
//---------------------------------------------------------
Arguments & Actions::operator[] (string name)
{
unsigned int i ;
for (i=0 ; i<arguments.size() ; i++)
    {
    if(arguments[i].commande==name) return arguments[i] ; 
    } 
cout << "WARN : Impossible de trouver un argument correspondant à la string " << name ;
Arguments * null_arg ;
null_arg=new Arguments ; 
return (*null_arg); 
}
//-------------------------------------------------------
double Arguments::operator[] (string name)
{
unsigned int i ;
for (i=0 ; i<p_names.size() ; i++)
    {
    if (p_names[i]==name) return params[i] ; 
    }
cout << "WARN : Impossible de trouver un paramètre correspondant à la string " << name ;
return -1 ; 
}
//---------------------------------------------------------
int Arguments::manual_set(string name, double valeur)
{
unsigned int i ;
for (i=0 ; i<p_names.size() ; i++)
    {
    if (p_names[i]==name) {params[i]=valeur ; return 0 ;} 
    }
cout << "WARN : Impossible de trouver un paramètre correspondant à la string // not set\n" << name ;
return -1 ; 
}
//------------------------------------------------------
int Actions::parse_args (int argc, char * argv[])
{
int i, j ; string commande ;
for (i=1; i<argc ; i++)
    {
    commande=argv[i] ; 
    if (commande[0]!='-' || commande[1]!='-')
    {actions.dumpnames.insert(actions.dumpnames.begin(),commande) ; continue ;}
    
    if (commande[1]=='-') commande.erase(0,2) ; 
    else commande.erase(0,1) ; 
 
    if (actions[commande].id==-1)  {DISP_Warn("WARN: argument de ligne de commande inconnu (verifier les deprecated argument peut-être?) : ") ; printf("%s\n", commande.c_str()) ; continue ; }
    actions[commande].set=true ; 

    if (commande ==  "filter" || commande == "F") 
       {
       Filter filtre_tmp ; vector <struct Op> tmpop ; 
       DISP_Warn ("NB : Filter est une commande fournie sans garantie d'être correcte.\n") ;
       tmpop=filtre_tmp.parse_arg(argv[i+1],TL) ; 
       actions.filtre_global.insert(actions.filtre_global.end(), tmpop.begin(), tmpop.end()) ; 
       i++ ; continue ; 
       }
    if (commande ==  "filterCF" || commande == "FCF") 
       {
       Filter filtre_tmp ; vector <struct Op> tmpop ; 
       DISP_Warn ("NB : Filter est une commande fournie sans garantie d'être correcte.\n") ;
       tmpop=filtre_tmp.parse_arg(argv[i+1],TCF) ;
       actions.filtre_global.insert(actions.filtre_global.end(), tmpop.begin(), tmpop.end()) ; 
       i++ ; continue ; 
       }   
    if (commande ==  "writing" || commande ==  "W") 
       {
       DISP_Info ("L'argument writing est encore en phase de test.\n") ;
       actions.ecrire.parse(argv[i+1]) ; 
       actions.ecrire.init=true ; 
       i++ ; continue ; 
       }
    if (commande ==  "alias") 
       {
       int res ; 
	 // TODO
	 //DISP_Err("L'argument alias n'est pas encore implémenté") ; 
       DISP_Info ("L'argument alias est encore en phase de test.\n") ;
       res=actions.parse_alias(argv[i+1], argc, argv) ;
       if (res>0) argc=res ; 
       i++ ; continue ; 
       }      
    if (commande ==  "bashcomplete") 
       {genere_bash_complete() ; continue ; }

    for (j=0 ; j<(int)actions[commande].params.size() ; j++)
        {sscanf(argv[i+j+1], "%lf", &(actions[commande].params[j])) ; fflush(stdout) ;  }
    i+=actions[commande].params.size() ; 

    if (actions[commande].commande=="chainforce")     actions.numdump++ ; 
    if (actions[commande].commande=="wallchainforce") actions.numdump++ ;
    }
    
if (actions.filtre_global.size()>0)
  Filter::disp(actions.filtre_global) ;      

if (check_dependences()>0)
   std::exit(EXIT_FAILURE) ;
return 1 ;
}
//----------------------------------------------------
int Actions::parse_args_fromfile (void)
{
char buffer[5000] ; int i=0, j, p1 ;   
char ** argv ; int argc=0 ; int change ; 
// Si cin n'est pas un terminal, on le considère comme un fichier d'argument
if (!isatty(fileno(stdin)))
{
  DISP_Info ("Lecture des arguments depuis l'entrée standard ...") ; 
  while (!cin.eof())
    {
     cin.read(buffer+i, 1);
     if (cin.eof())
          *(buffer+i)=0 ; 
    i++ ; if (i>5000) DISP_Err("Buffer trop petit pour l'entree standard") ; 
    }
  // Splitting
  p1=0 ; change=0 ; 
  argc=1 ; // Le premier argument n'est pas utilisé (dans la ligne de commande, c'est le nom du programme)
  argv=(char **) malloc(argc*sizeof(char*)) ; //argv[0]=NULL ; 
  for (j=0 ; j<i ; j++)
  {
    if (buffer[j]==' ' || buffer[j]=='\n' || buffer[j]==0 || buffer[j]=='\t')
    {
      if (change==1) 
      {
	argc++ ; 
	argv=(char **) realloc(argv, argc*sizeof(char*)) ;
	argv[argc-1]=buffer+p1 ; 
	buffer[j]=0 ;
	change=0 ; 
      }
    }
    else
    {
      if (change==0) {change=1 ; p1=j ; }
    }
  }
  parse_args(argc, argv) ;  
}
return 1 ; 
}
//-------------------------------------------------
int Actions::parse_alias (char com[], int argc, char * argv[])
{
 int len=0, column, num, i,j,k; int res; 
 char *action, *nom ;
 vector <string> aliasname ; vector <string> aliascommand ; 
 FILE *aliasfile ; char aliaspath[1000] ; 
 bool sw_nom ; 
 action=com ; 
 while (com[len]!=0)
 {
   if (com[len]==':') {com[len]=0 ; nom=com+len+1 ; }
   len++ ;
 }

 char *temp  = getenv("HOME");
 if (temp==NULL) {DISP_Warn("La variable environnement $HOME n'est pas accessible, impossible d'enregistrer/lire un fichier d'alias de Post Processing\n") ; return -1;}
 aliaspath[0]=0 ; 
 strcat(aliaspath, temp) ; 
 strcat(aliaspath, "/bin/PostProcessing.alias") ; 
 aliasfile=fopen(aliaspath, "r") ; 
 if (aliasfile!=NULL)
 { // Reading alias file
   char ligne[10000] ; 
   num=0 ; 
   sw_nom=true ; 
   while (!feof(aliasfile))
   {
    res=fscanf(aliasfile,"%[^\n]", ligne) ; res=fscanf(aliasfile,"%*c") ;
    if (feof(aliasfile)) break ; 
    if (sw_nom==true) 
    {
     sw_nom=false ; num++ ; 
     aliasname.push_back(ligne) ; 
    }
    else
    {
      sw_nom=true ;       
      aliascommand.push_back(ligne) ;  
    } 
  }
  fclose(aliasfile) ; 
 }
 
 if (!strcmp(com, "set"))
 {
  for (i=0 ; i<num ; i++) {if (aliasname[i]==nom) {DISP_Warn("L'alias existe déjà. Impossible de l'ajouter.") ; return -2 ;}}
  aliasname.push_back(nom) ; 
  aliascommand.push_back("") ;
  for (i=1 ; i<argc-1 ; i++)
  {
    if (!strcmp(argv[i], "--alias")) {i++ ; continue ; }
    else if (!strcmp(argv[i], "--chainforce")) {argc-- ; }
    else if (!strcmp(argv[i], "--vtk")) {argc-- ; }
    else if(!strcmp(argv[i], "--wallchainforce")) {argc-- ; }
    
    aliascommand[num].append(argv[i]) ; 
    if (i!= argc-2) aliascommand[num].append(" ") ; 
  }
  num++ ; 
 }
 else if (!strcmp(com, "reset"))
 {
  for (i=0 ; i<num ; i++) {if (aliasname[i]==nom) break ; }
  if (i==num) {DISP_Warn("Impossible de changer un alias qui n'existe pas. Utiliser set.\n") ; return -3 ;}
  aliascommand[i]=("") ;
  for (i=1 ; i<argc-1 ; i++)
  {
    if (!strcmp(argv[i], "--alias")) {i++ ; continue ; }
    else if (!strcmp(argv[i], "--chainforce")) {argc-- ; }
    else if (!strcmp(argv[i], "--vtk")) {argc-- ; }
    else if(!strcmp(argv[i], "--wallchainforce")) {argc-- ; }
    
    aliascommand[i].append(argv[i]) ; 
    if (i!= argc-2) aliascommand[i].append(" ") ; 
  }
  i++ ; 
 }
 else if (!strcmp(com, "del"))
 {
   for (i=0 ; i<num-1 ; i++)
     if (aliasname[i]==nom) 
        {aliasname.erase(aliasname.begin()+i) ; aliascommand.erase(aliascommand.begin()+i) ; break ;}
   num-- ;   
 }
 else if (!strcmp(com, "call"))
 {
  for (i=0 ; i<num ; i++) {if (aliasname[i]==nom) break ; }
  if(i==num) {DISP_Warn("Impossible de charger l'alias qui n'existe pas") ; return -4 ; }
 
  // On commence par recopier argv !
  char ** argv2 ;
  argv2=(char**)malloc(sizeof(char*)*argc) ; 
  for (j=0 ; j<argc ; j++) argv2[j]=argv[j] ; 
  argv=argv2 ; 
 
  char argument[1000] ; 
  size_t found ; 
  argument[0]='a' ; argument[1]=0 ; 
  for (j=0 ; j<argc ; j++) {if (!strcmp("--alias", argv[i])) break ; }
  j+=2 ;
 
  while (strlen(argument)!=0)
  {
   sscanf(aliascommand[i].c_str(), "%s", argument) ; 
   found=aliascommand[i].find_first_of(" ",0);
   aliascommand[i].erase(aliascommand[i].begin(), aliascommand[i].begin()+found+1) ; 
   argc++ ;
   argv=(char**)realloc(argv, sizeof(char*)*(argc)) ; 
   for (k=argc-1 ; k>j ; k--) 
     argv[k]=argv[k-1] ; 
   argv[j]=(char*)malloc(sizeof(char)*strlen(argument)) ; 
   strcpy(argv[j], argument) ; 
  }
 }
 else if (!strcmp(com, "list"))
 {
  printf("Listing des alias définis\n") ; 
  for (i=0 ; i<num ; i++)
    printf("%s\n  => %s", aliasname[i].c_str(), aliascommand[i].c_str()) ; 
 }
 //else if (!strcmp(com, "info"))
 else {DISP_Warn("Argument d'alias inconnu") ; return -5 ; }
 
 // Writing du fichier d'alias 
 aliasfile=fopen(aliaspath, "w") ; 
 for (i=0 ; i<num ; i++)
   fprintf(aliasfile, "%s\n%s\n", aliasname[i].c_str(), aliascommand[i].c_str()) ; 
 fclose(aliasfile) ; 
 
 com[strlen(com)]=':' ;
 return argc ; 
}

//-----------------------------------------------------
int Actions::copyarg (string from, string to)
{
int id1, id2 ; unsigned int i ;
for (i=0 ; i<arguments.size() ; i++)
    {
    if(arguments[i].commande==from) id1=i ;
    if(arguments[i].commande==to) id2=i ;
    }

arguments[id2].set=arguments[id1].set ;
arguments[id2].params=arguments[id1].params ;

for (i=0 ; i<arguments[id2].p_names.size() ; i++)
	{if (arguments[id2].p_names[i]!=arguments[id1].p_names[i]) DISP_Warn("Les paramètres des arguments ne sont pas les mêmes, la copie est plus que discutable ...") ; }

return 0 ;
}

//-----------------------------------------------------
int Actions::check_dependences (void)
{
// TODO : il faudrait tester que l'on a pas de boucle dans les dépendances. Mais j'ai la flemme de le programmer maintenant et c'est pas si simple non plus ...
unsigned int i, j ;
int err=0 ;
for (i=0 ; i<arguments.size() ; i++)
 {
 if (arguments[i].set)
  { 
  for (j=0 ; j<arguments[i].dependances.size() ; j++)
      {
      if (arguments[arguments[i].dependances[j]-1].set==false)
         {cout << "ERR : l'argument "<< arguments[i].commande << " nécessite l'argument " << arguments[arguments[i].dependances[j]-1].commande << ".\n" ; err++ ; }
      }
  }
 }
return err ; 
}
//-----------------------------------------------------
void Actions::aide (void) //TODO
{
unsigned int i, j ;

for (i=0 ; i<arguments.size() ; i++)
{
	cout << "\n\t[" << arguments[i].id << "]" ;
	cout << "--" << arguments[i].commande << " " ;
	for (j=0 ; j<arguments[i].p_names.size() ; j++)
        {cout << arguments[i].p_names[j] << " " ; }
	cout << ":\n\t\t" << arguments[i].description ;
	if (arguments[i].dependances.size()>0)
	{cout << "\n\t\tDep : " ; for (j=0 ; j<arguments[i].dependances.size() ; j++) {cout <<"|"<<arguments[i].dependances[j];} cout <<"|" ;}
}
cout << "\n\n" ;
}

//----------------------------------------------------
void Actions::write_history (int todo, int argc, char * argv[])
{
time_t t;
static time_t told ;
int i, j ;
struct tm * timeinfo;
FILE *log ;
char *path=NULL;
char logpath[5000] ; 
size_t size=0;

time(&t);
timeinfo = localtime (&t);

char *temp  = getenv("HOME");
if (temp==NULL) {DISP_Warn("La variable environnement $HOME n'est pas accessible, impossible d'enregistrer un log de lancement de Post Processing\n") ; return ;}
logpath[0]=0 ; 
strcat(logpath, temp) ; 
strcat(logpath, "/bin/PostProcessing.log") ;
log=fopen(logpath, "a") ;
if (log==NULL) {printf("WARN : impossible d'enregistrer un log de lancement de PostProcessing") ; return ; }

if (todo==0)
  {
  path=getcwd(path,size);
  fprintf(log, "%s", ctime(&t)) ;
  fprintf(log, "\t%s compilé le %s\n", argv[0], __DATE__) ;
  fprintf(log, "\t Travaille dans %s\n\t", path) ;
  for (i=0 ; i<actions.arguments.size() ; i++)
  {
   if (actions.arguments[i].set)
   {
    fprintf(log, "--%s ",  actions.arguments[i].commande.c_str()) ; 
    for (j=0 ; j<actions.arguments[i].params.size() ; j++)
      fprintf(log, "%f ", actions.arguments[i].params[j]) ; 
   }
  }
  for (i=actions.numdump-1 ; i>=0 ; i--)
    fprintf(log, "%s ", actions.dumpnames[i].c_str()) ; 
  fprintf(log, "\n") ;
  told=t ;
  }
else if (todo==1)
  {
  fprintf(log, "Fin normale (durée %ld:%ld:%ld).\n\n", (t-told)/3600, ((t-told)%3600)/60, (((t-told)%3600)%60)) ;
  }
free(path) ;   
fclose(log) ;
}

//-------------------------------------------------------------
void Actions::regenerate_restart2vtk_script()
{
  char chemin[500] ; 
  FILE *Script ; 
  
  FILE *fp = popen("which gawk", "r");
  while(fgets(chemin, 500, fp)!=NULL){}
  fclose(fp);
  chemin[strlen(chemin)-1]=0 ; 
  
  Script=fopen("restart2vtk", "r") ; 
  if (Script != NULL) {DISP_Err("Le fichier restart2vtk existe deja dans le dossier courant. Supprimez-le pour le recréer.\n") ; return ;}
  
  Script=fopen("restart2vtk", "w") ; 
  fprintf(Script,"#!%s -f\n\nBEGIN{start=0 ;i=1 ;}\n\n/atoms/ {nbatm=$1}\n/Atoms/ {start=1}\n\n", chemin) ; 
  fprintf(Script,"{\nif ((NF==7 && start==1))\n{\nid[i]=$1 ; type[i]=$2; rad[i]=$3 ; rho[i]=$4 ; x[i]=$5 ; y[i]=$6 ; z[i]=$7 ;\ni=i+1 ;\n}\n}\n\n") ; 
  fprintf(Script,"END{print \"# vtk DataFile Version 2.0\" ; print \"Made by PostProcessing (CC-BY-NC)\" ;\n print \"ASCII\" ; \n \
print \"DATASET POLYDATA\" ; print \"POINTS\",nbatm,\"double\" ; \n for (i=1;i<=nbatm;i++) print x[i],y[i],z[i] ; \n \
print \"VERTICES\",nbatm,2*nbatm ; \n for (i=1;i<=nbatm;i++) print \"1\",i-1 ; \n print \"POINT_DATA\", nbatm ; \n print \"SCALARS rayon double 1\";\
\nprint \"LOOKUP_TABLE default\" ; \nfor (i=1;i<=nbatm;i++) print rad[i] ;\n print \"SCALARS id double 1\";\nprint \"LOOKUP_TABLE default\" ; \
\nfor (i=1;i<=nbatm;i++) print id[i] ;\n}") ; 
  fclose(Script) ;
  
  int res ; 
  res=system("chmod +x restart2vtk") ; 
  DISP_Info("Le script restart2vtk a été créé et est prêt à être utilisé\n") ; 
}

//-------------------------------------------------------------
void Actions::genere_bash_complete (void)
{
  FILE *out ; int i, j ; 
  
  out=fopen("bashcomplete.sh", "w") ; 
  
  // beginning
  fprintf(out, "_PPcomplete()\n{\n  local cur prev opts\n  COMPREPLY=()\n  cur=\"${COMP_WORDS[COMP_CWORD]}\"\n  prev=\"${COMP_WORDS[COMP_CWORD-1]}\"\n");
  fprintf(out, "  opts=\"") ;
  
  for (i=0 ; i<arguments.size() ; i++)
    fprintf(out, "--%s ", arguments[i].commande.c_str()) ; 
  fprintf(out, "\"\n") ; 
  
  fprintf(out, "\n\n  case \"${prev}\" in\n") ; 
  for (i=0 ; i<arguments.size() ; i++)
  {
   if (arguments[i].p_names.size() >0)
   {
     fprintf(out, "    --%s) \n      COMPREPLY=( [%s] ", arguments[i].commande.c_str(), arguments[i].commande.c_str()) ; 
     for(j=0 ; j<arguments[i].p_names.size() ; j++)
      fprintf(out, "%s ", arguments[i].p_names[j].c_str()) ; 
     fprintf(out, ")\n      return 0\n      ;;\n") ; 
   }
  }
  fprintf(out, "    *)\n      ;;\n  esac") ; 
  
  fprintf(out, "\n\n  if [[ ${cur} == -* ]] ; then\n  COMPREPLY=( $(compgen -W \"${opts}\" -- ${cur}) )\n  return 0\n  fi\n\n  COMPREPLY=( $(compgen -f ${cur}) ) \n}") ; 
  
  fclose(out) ; 
}





