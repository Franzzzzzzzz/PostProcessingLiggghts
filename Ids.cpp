#include "Headers/Ids.h"

void IDSCLASS::default_assign()
{
idname.resize(256, "") ;
vector <string> nullstr ;
idalias.resize(256, nullstr) ;


idname[1]="ID" ; idname[2]="TYPE" ;
idname[3]="POSX" ; idname[4]="POSY" ; idname[5]="POSZ" ;
idname[6]="VX" ; idname[7]="VY" ; idname[8]="VZ" ;
idname[9]="RAYON" ;
idname[10]="FX" ; idname[11]="FY" ; idname[12]="FZ" ;
idname[13]="MASSE" ;
idname[14]="OMEGAX" ; idname[15]="OMEGAY" ; idname[16]="OMEGAZ" ;
idname[17]="SIGMAXX" ; idname[18]="SIGMAXY" ; idname[19]="SIGMAXZ" ;
idname[20]="SIGMAYX" ; idname[21]="SIGMAYY" ; idname[22]="SIGMAYZ" ;
idname[23]="SIGMAZX" ; idname[24]="SIGMAZY" ; idname[25]="SIGMAZZ" ;
idname[26]="FORCEWALLX" ; idname[27]="FORCEWALLY" ; idname[28]="FORCEWALLZ" ;
idname[29]="ATMPRESSURE" ; idname[30]="IDMULTISPHERE" ;
idname[31]="KX" ; idname[32]="KY" ; idname[33]="KZ" ;
idname[34]="CENTROIDX" ; idname[35]="CENTROIDY" ; idname[36]="CENTROIDZ" ;
idname[37]="QUATERNION1" ; idname[38]="QUATERNION2" ; idname[39]="QUATERNION3" ; idname[40]="QUATERNION4" ;

//idname[63]="CELLDATAMIN" ;
idname[64]="CFID1" ; idname[65]="CFID2" ;
idname[66]="CFFORCEX" ; idname[67]="CFFORCEY" ; idname[68]="CFFORCEZ" ;
idname[69]="CFMAG" ;
idname[70]="CFR" ; idname[71]="CFTHETA" ; idname[72]="CFPHI" ;
idname[73]="CFX" ; idname[74]="CFY" ; idname[75]="CFZ" ;
idname[76]="CFPERIOD" ;
idname[77]="CFID1X" ; idname[78]="CFID1Y" ; idname[79]="CFID1Z" ;
idname[80]="CFID2X" ; idname[81]="CFID2Y" ; idname[82]="CFID2Z" ;
idname[83]="CFPOSX" ; idname[84]="CFPOSY" ; idname[85]="CFPOSZ" ;
//idname[127]="CELLDATAMAX" ;

idname[128]="PRESSURE" ; idname[129]="SHEARSTRESS" ;
idname[130]="FORCEX" ; idname[131]="FORCEY" ; idname[132]="FORCEZ" ;
idname[133]="NORMALEX" ; idname[134]="NORMALEY" ; idname[135]="NORMALEZ" ;
idname[136]="POINTX" ; idname[137]="POINTY" ; idname[138]="POINTZ" ;
idname[139]="POLY1" ; idname[140]="POLY2" ; idname[141]="POLY3" ; idname[142]="POLY4" ;
idname[143]="STRESSX" ; idname[144]="STRESSY" ; idname[145]="STRESSZ" ;
idname[146]="CENTREX" ; idname[147]="CENTREY" ; idname[148]="CENTREZ" ;

idname[192]="COARPHI" ;      idname[193]="COARATM" ;      idname[194]="COARCONTACTS" ; idname[195]="COARRAD" ;
idname[196]="COARSIGKXX" ;   idname[197]="COARSIGKXY" ;   idname[198]="COARSIGKXZ" ;
idname[199]="COARSIGKYX" ;   idname[200]="COARSIGKYY" ;   idname[201]="COARSIGKYZ" ;
idname[202]="COARSIGKZX" ;   idname[203]="COARSIGKZY" ;   idname[204]="COARSIGKZZ" ;
idname[205]="COARSIGTOTXX" ; idname[206]="COARSIGTOTXY" ; idname[207]="COARSIGTOTXZ" ;
idname[208]="COARSIGTOTYX" ; idname[209]="COARSIGTOTYY" ; idname[210]="COARSIGTOTYZ" ;
idname[211]="COARSIGTOTZX" ; idname[212]="COARSIGTOTZY" ; idname[213]="COARSIGTOTZZ" ;
idname[214]="COARSTDPRES"  ; idname[215]="COARSTDTAU" ;   idname[216]="COARSTDPHI" ;
idname[217]="COARORIENTXX" ; idname[218]="COARORIENTXY" ; idname[219]="COARORIENTXZ" ;
idname[220]="COARORIENTYX" ; idname[221]="COARORIENTYY" ; idname[222]="COARORIENTYZ" ;
idname[223]="COARORIENTZX" ; idname[224]="COARORIENTZY" ; idname[225]="COARORIENTZZ" ;
idname[226]="COARNBGP" ;

idname[255]="UNKNOWN" ;

idalias[1].push_back("id") ; idalias[2].push_back("type") ;
idalias[13].push_back("masse") ; idalias[13].push_back("mass") ;
idalias[9].push_back("radius") ; idalias[9].push_back("rayon") ; 
idalias[3].push_back("x") ; idalias[4].push_back("y") ; idalias[5].push_back("z") ;
idalias[3].push_back("posx") ; idalias[4].push_back("posy") ; idalias[5].push_back("posz") ;
idalias[3].push_back("xu") ; idalias[4].push_back("yu") ; idalias[5].push_back("zu") ;
idalias[6].push_back("vx") ; idalias[7].push_back("vy") ; idalias[8].push_back("vz") ;
idalias[10].push_back("fx") ; idalias[11].push_back("fy") ; idalias[12].push_back("fz") ;
idalias[14].push_back("omegax") ; idalias[15].push_back("omegay") ; idalias[16].push_back("omegaz") ;
idalias[17].push_back("sigmaxx") ; idalias[18].push_back("sigmaxy") ; idalias[19].push_back("sigmaxz") ;
idalias[20].push_back("sigmayx") ; idalias[21].push_back("sigmayy") ; idalias[22].push_back("sigmayz") ;
idalias[23].push_back("sigmazx") ; idalias[24].push_back("sigmazy") ; idalias[25].push_back("sigmazz") ;
idalias[29].push_back("atmpressure") ;
idalias[30].push_back("id_multisphere") ; idalias[30].push_back("idmultisphere") ; idalias[30].push_back("multisphere") ; idalias[30].push_back("group") ;
idalias[31].push_back("kx") ; idalias[32].push_back("ky") ; idalias[33].push_back("kz") ;
idalias[34].push_back("centroidx") ; idalias[35].push_back("centroidy") ; idalias[36].push_back("centroidz") ;
idalias[17].push_back("sigxx") ; idalias[18].push_back("sigxy") ; idalias[19].push_back("sigxz") ;
idalias[20].push_back("sigyx") ; idalias[21].push_back("sigyy") ; idalias[22].push_back("sigyz") ;
idalias[23].push_back("sigzx") ; idalias[24].push_back("sigzy") ; idalias[25].push_back("sigzz") ;
idalias[26].push_back("forcewallx") ; idalias[27].push_back("forcewally") ; idalias[28].push_back("forcewallz") ;
idalias[26].push_back("f_force_cyl[1]") ; idalias[27].push_back("f_force_cyl[2]") ; idalias[28].push_back("f_force_cyl[3]") ;
idalias[26].push_back("f_contactforces_cyl[2]") ; idalias[27].push_back("f_contactforces_cyl[3]") ; idalias[28].push_back("f_contactforces_cyl[4]") ;
idalias[1].push_back("f_contactforces_cyl[1]") ;
idalias[37].push_back("quat1") ; idalias[38].push_back("quat2") ; idalias[39].push_back("quat3") ; idalias[40].push_back("quat4") ;

idalias[64].push_back("cfid1") ; idalias[65].push_back("cfid2") ;
idalias[66].push_back("cfforcex") ; idalias[67].push_back("cfforcey") ; idalias[68].push_back("cfforcez") ;
idalias[66].push_back("cffx") ; idalias[67].push_back("cffy") ; idalias[68].push_back("cffz") ;
idalias[69].push_back("cfmag") ; idalias[70].push_back("cfr") ;  idalias[71].push_back("cftheta") ; idalias[72].push_back("cfphi") ;
idalias[73].push_back("cfx") ; idalias[74].push_back("cfy") ;  idalias[75].push_back("cfz") ;
idalias[76].push_back("cfperiod") ; idalias[76].push_back("periodicity") ;
idalias[83].push_back("cfposx") ; idalias[84].push_back("cfposy") ;  idalias[85].push_back("cfposz") ;
idalias[64].push_back("c_cout[1]") ; idalias[65].push_back("c_cout[2]") ;
idalias[66].push_back("c_cout[4]") ; idalias[67].push_back("c_cout[5]") ; idalias[68].push_back("c_cout[6]") ;
idalias[76].push_back("c_cout[3]") ;
idalias[64].push_back("c_contacts[7]") ; idalias[65].push_back("c_contacts[8]") ;
idalias[66].push_back("c_contacts[10]") ; idalias[67].push_back("c_contacts[11]") ; idalias[68].push_back("c_contacts[12]") ;
idalias[76].push_back("c_contacts[9]") ; idalias[76].push_back("periodicity") ;
idalias[77].push_back("c_contacts[1]") ; idalias[78].push_back("c_contacts[2]") ; idalias[79].push_back("c_contacts[3]") ;
idalias[80].push_back("c_contacts[4]") ; idalias[81].push_back("c_contacts[5]") ; idalias[82].push_back("c_contacts[6]") ;

idalias[64].push_back("c_allcontacts[7]") ; idalias[65].push_back("c_allcontacts[8]") ;
idalias[66].push_back("c_allcontacts[10]") ; idalias[67].push_back("c_allcontacts[11]") ; idalias[68].push_back("c_allcontacts[12]") ;
idalias[76].push_back("c_allcontacts[9]") ;
idalias[77].push_back("c_allcontacts[1]") ; idalias[78].push_back("c_allcontacts[2]") ; idalias[79].push_back("c_allcontacts[3]") ;
idalias[80].push_back("c_allcontacts[4]") ; idalias[81].push_back("c_allcontacts[5]") ; idalias[82].push_back("c_allcontacts[6]") ;


idalias[133].push_back("normalex") ; idalias[134].push_back("normaley") ; idalias[135].push_back("normalez") ;
idalias[146].push_back("centrex") ; idalias[147].push_back("centrey") ; idalias[148].push_back("centrez") ;

idalias[196].push_back("sigmakxx") ;   idalias[197].push_back("sigmakxy") ;   idalias[198].push_back("sigmakxz") ;
idalias[199].push_back("sigmakyx") ;   idalias[200].push_back("sigmakyy") ;   idalias[201].push_back("sigmakyz") ;
idalias[202].push_back("sigmakzx") ;   idalias[203].push_back("sigmakzy") ;   idalias[204].push_back("sigmakzz") ;
idalias[205].push_back("sigmatotxx") ; idalias[206].push_back("sigmatotxy") ; idalias[207].push_back("sigmatotxz") ;
idalias[208].push_back("sigmatotyx") ; idalias[209].push_back("sigmatotyy") ; idalias[210].push_back("sigmatotyz") ;
idalias[211].push_back("sigmatotzx") ; idalias[212].push_back("sigmatotzy") ; idalias[213].push_back("sigmatotzz") ;
idalias[195].push_back("mean_radius") ; idalias[193].push_back("nb_atomes") ; idalias[194].push_back("nb_contacts") ; idalias[192].push_back("phi") ;
idalias[214].push_back("stdpressure") ; idalias[215].push_back("stdtau") ; idalias[216].push_back("stdphi") ; //! STDPRESSURE et STDTAU, en fait c'est l'écart-type de l'estimation de la moyenne ...
idalias[217].push_back("orientxx") ; idalias[218].push_back("orientxy") ; idalias[219].push_back("orientxz") ;
idalias[220].push_back("orientyx") ; idalias[221].push_back("orientyy") ; idalias[222].push_back("orientyz") ;
idalias[223].push_back("orientzx") ; idalias[224].push_back("orientzy") ; idalias[225].push_back("orientzz") ;
idalias[226].push_back("nbgp") ;

idalias[255].push_back("null") ;

}

//------------------------------
vector <string> IDSCLASS::alias(int id)
{vector <string> retour ;

  if (id<0 || id>=idalias.size()) return retour ;
  else
  {
    retour=idalias[id] ;
    return retour ;
  }
}

//------------------------------
int IDSCLASS::new_id(string name, int type)
{
int min, max ; int i ;
string canonic_name ; char info[500] ;
// Est-ce que la chaîne existe déjà (en ID ou alias)
  if (find(name) != -1 ) {DISP_Warn(_("IDS : la chaîne ou l'alias existe déjà. Impossible de l'ajouter\n")) ; return -1 ; }

  canonic_name=name ;
  canonic_name.erase(remove_if( canonic_name.begin(), canonic_name.end(), ::ispunct ), canonic_name.end() );
  canonic_name.erase(remove_if( canonic_name.begin(), canonic_name.end(), ::isspace ), canonic_name.end() );
  std::transform(canonic_name.begin(), canonic_name.end(),canonic_name.begin(), ::toupper);

  switch (type)
  {
    case TL : min=1 ; max=63 ; break ;
    case TCF : min=64 ; max=127 ; canonic_name = "CF" + canonic_name ;  break ;
    case TF : min=128 ; max=191 ; break ;
    case TCOARSE : min=192 ; max=254 ; canonic_name = "COAR" + canonic_name ;  break ;
    case TNONE : min=1 ; max=254 ; DISP_Warn(_("IDS : le type TNONE doit être utilisé avec précaution")) ;
    default : DISP_Err(_("IDS : type inconnu")) ; break ;
  }

  for (i=min ; i<max ; i++)
  {
    if (idname[i]=="") // on a trouvé une place vide
    {
      idalias[i].push_back(name) ;
      idname[i]=canonic_name ;
      break ;
    }
  }

sprintf(info, _("L'alias %s a été ajouté avec succès sous le nom %s avec l'id %d.\n"), idalias[i][0].c_str(), idname[i].c_str(), i) ;
DISP_Info (info) ;
return i ;
}

//------------------------------
int IDSCLASS::display_all ()
{
  printf(_("List of ID names and aliases as defined in IDS\n----------------------------------------------\n")) ;
  for (int i=0 ; i<idname.size() ; i++)
  {
   printf("[%3d] %s\n", i, idname[i].c_str() ) ;
   for (int j=0 ; j<idalias[i].size() ; j++)
     printf("\t%s\n",idalias[i][j].c_str()) ;
  }
  return 0 ;
}
