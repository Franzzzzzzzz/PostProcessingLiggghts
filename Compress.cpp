#include "Headers/Compress.h"
#include <math.h>

int Compresser::compress(Dump &dump, string nom, FILE * out)
{
  char header[27] ; unsigned int i,j,k,l ;
  unsigned int longueur ;
  vector <unsigned char> fmts ;
  float F ; unsigned int UI ; unsigned char UC ; unsigned short int USI ; double D ; 
  vector <unsigned char> UCdatas ; 
  unsigned short int U16;
  
  //if (!check_all(dump))
  //   {cout << "ERREUR DE COMPRESSION : le dump n'a pas les caractéristiques nécessaires pour être compressé" ; exit ; }
  
  strcpy(header,"AVFFLIGGGHTSDUMPCOMPRES2.3") ; 
  
  // Ecriture du header d'identification et du nom du fichier
  fwrite(header, 1, 26, out) ; 
  U16= nom.length() ;  
  fwrite(&U16,2,1,out) ;
  fwrite(nom.c_str(), 1, U16, out) ; 
  
  // Ecriture des propriétés du dump
  fwrite(&dump.nbsteps, 4,1,out) ; //nb de pas de temps
  dump.check_timestep(0) ; 
  fwrite(&dump.steps[0].nb_idx, 1,1,out) ;
  
  champs.compute_global(dump) ; 
  dump.check_timestep(0) ; 
  fmts.resize(dump.steps[0].nb_idx) ; 

  for (i=0 ; i<dump.steps[0].nb_idx ; i++)
  {
  fwrite(&dump.steps[0].idx_col[i], 1,1,out) ;
  fmts[i]=format(i) ;
  fwrite(&fmts[i],1,1,out) ; 
  
  if ((fmts[i]&(~MASK_ALWAYS_THE_SAME))>=17) // Il s'agît d'un format dénombrable. Il faut écrire les valeurs.
     {
      unsigned char bits ; 
      bits=(fmts[i]&(~MASK_ALWAYS_THE_SAME))>>4 ;
      switch(bits)
        {
        case 2 : bits=2 ; break ; 
   	case 3 : bits=4 ; break ;
   	case 4 : bits=16 ;  break ; 
 	}
      for (j=0 ; j<champs.denom[i] ; j++)
        {ecrirefmt(champs.denombrables[i][j], (fmts[i]&(~MASK_ALWAYS_THE_SAME)), out) ;}
      for (j=champs.denom[i] ; j<bits ; j++)
        {ecrirefmt(0.0, (fmts[i]&(~MASK_ALWAYS_THE_SAME)), out) ;}
     }
  }
  
  // tout ce qu'il fallait en header a été écrit. On continue avec le corps réel du fichier
  actions.total=dump.nbsteps ; 
  cout << "\nCompresser::compress          " ;
  actions.disp_progress() ; 
  for (i=0 ; i<(unsigned int)dump.nbsteps ; i++)
    {
    actions.valeur=i ; 
    dump.check_timestep(i) ;  
    if (dump.steps[i].Type==TL || dump.steps[i].Type==TCF) 
       {
       fwrite(&dump.steps[i].timestep, 4, 1, out) ;
       fwrite(&dump.steps[i].nb_atomes, 4, 1, out) ; 
       }
    else
       {
       fwrite(&dump.steps[i].nb_pts, 4, 1, out) ;
       fwrite(&dump.steps[i].nb_triangles, 4, 1, out) ; 
       }
    if (dump.steps[i].Type==TL)
     {
     for (j=0 ; j<3 ; j++)
      {
      for (k=0 ; k<2 ; k++)
        {
	ecrirefmt(dump.steps[i].box[j][k], FLOAT, out) ; 
        }
      }
     }
    else
     {for (j=0 ; j<6 ; j++) {ecrirefmt(0.0, FLOAT, out) ;}}
    
    for (j=0; j<champs.idx_col.size() ; j++)
       {
       fmts[j]=fmts[j]&(~MASK_ALWAYS_THE_SAME) ; 
       if (champs.same[j]==false || i==0)
         {
         if (fmts[j]==UCHAR)
          {
	  for (l=0 ; l<dump.steps[i].datas[j].size() ; l++) {UC=dump.steps[i].datas[j][l] ; fwrite(&UC, 1, 1, out) ; } 
          }
         if (fmts[j]==USINT)	   
          {
	  for (l=0 ; l<dump.steps[i].datas[j].size() ; l++) {USI=dump.steps[i].datas[j][l] ; fwrite(&USI, 2, 1, out) ; }
          }
         else if (fmts[j]==UINT)
	  {
	   for (l=0 ; l<dump.steps[i].datas[j].size() ; l++) {UI=dump.steps[i].datas[j][l] ; fwrite(&UI, 4,1, out) ;}
	  }
         else if (fmts[j]==FLOAT)
	  {
	  for (l=0 ; l<dump.steps[i].datas[j].size() ; l++) {F=dump.steps[i].datas[j][l] ; fwrite(&F, 4, 1, out) ;}
	  }
         else if (fmts[j]==DOUBLE)
          {
          for (l=0 ; l<dump.steps[i].datas[j].size() ; l++) {D=dump.steps[i].datas[j][l] ; fwrite(&D, 8, 1, out) ;}
          }
         else if (fmts[j]>=33) // C'est un format dénombrable
          {longueur=create_datas_denom (UCdatas, fmts[j], j, dump.steps[i].datas[j]) ; 
	   for (l=0 ; l<longueur ; l++) {fwrite(&(UCdatas[l]),1,1,out) ;} 
	  }
	}
       }
    }
return 1 ;
}

//-----------------------------------------
Donnees Compresser::ecrirefmt (double valeur, int fmt, FILE * out)
{
Donnees tmp ; 
int format ;
format=fmt & 15 ;

if ((format==1) || (format==2))
  {tmp.longueur=1 ; tmp.C=valeur ; fwrite (&tmp.C, 1, 1, out) ; }
else if ((format ==3) || (format==4))
  {tmp.longueur=2 ; tmp.SI=valeur ; fwrite(&tmp.SI, 2, 1, out);}
else if ((format==5) || (format==6) || (format==9))
  {tmp.longueur=4 ;
   if (format==9)
      {tmp.F=valeur ; fwrite(&tmp.F, 4, 1, out);}
   else
      {tmp.I=valeur ; fwrite(&tmp.I, 4, 1, out);}
  }
else if ((format==7) || (format==8) || (format==10))
  {tmp.longueur=8 ; 
   if (format==10)
     {tmp.D=valeur ; fwrite(&tmp.D, 8, 1, out);}
   else
     {tmp.LI=valeur ; fwrite(&tmp.LI, 8, 1, out);} 
  }
else
{cout << _("Format inconnu dans la sous fonction convertDBL"); std::exit(EXIT_FAILURE) ;}

return tmp ; 
}

//-----------------------------------------
char Compresser::format (int i)
{  
unsigned char mask ; 
if (champs.same[i])
   mask=MASK_ALWAYS_THE_SAME ; 
else
   mask=0 ; 

  if (champs.idx_col[i]==IDS("POLY1") || champs.idx_col[i]==IDS("POLY2") ||champs.idx_col[i]==IDS("POLY3") ||champs.idx_col[i]==IDS("POLY4") )
    cout <<"" ;
  else if (champs.idx_col[i]==IDS("ID")) // C'est un unsigned integer
  {
      if (champs.denom[i]==1)
	return(UINT_CST|mask) ;
      else if (champs.denom[i]==2)
	return(UINT_DENOM_2|mask) ;
      else if (champs.denom[i]>2 && champs.denom[i]<=4)
	return(UINT_DENOM_4|mask) ; 
      else if (champs.denom[i]>4 && champs.denom[i]<=16)
	return(UINT_DENOM_16|mask) ; 
      else if (champs.max[i]<65536)
	return(USINT|mask) ;
      else
	return(UINT|mask) ;
  }
  else if (champs.idx_col[i]==IDS("TYPE")) // C'est un unsigned char
  {
      if (champs.denom[i]==1)
	return(UCHAR_CST|mask) ;
      else if (champs.denom[i]==2)
	return(UCHAR_DENOM_2|mask) ;
      else if (champs.denom[i]>2 && champs.denom[i]<=4)
	return(UCHAR_DENOM_4|mask) ; 
      else if (champs.denom[i]>4 && champs.denom[i]<=16)
	return(UCHAR_DENOM_16|mask) ; 
      else
	return(UCHAR) ;
  }
  else if (champs.idx_col[i]==IDS("CFID1") || champs.idx_col[i]==IDS("CFID2"))
  {   if (champs.denom[i]==1)
	return(UINT_CST|mask) ;
      else if (champs.denom[i]==2)
	return(UINT_DENOM_2|mask) ;
      else if (champs.denom[i]>2 && champs.denom[i]<=4)
	return(UINT_DENOM_4|mask) ; 
      else if (champs.denom[i]>4 && champs.denom[i]<=16)
	return(UINT_DENOM_16|mask) ; 
      else if (champs.max[i]<65536)
	return(USINT|mask) ;
      else
	return(UINT|mask) ;
  }
  else if (champs.idx_col[i]==IDS("CFPERIOD"))
  {
      if (champs.denom[i]==1)
	return(UINT_CST|mask) ;
      else if (champs.denom[i]==2)
	return(UINT_DENOM_2|mask) ;
      else
        cout << _("This is strange : in compression, CFPERIOD should be 0 or 1 ! ") ; 
  }
  else if (champs.idx_col[i]==IDS("UNKNOWN"))
  {
    cout << _("ERREUR de compression : le type UNKNOWN ne peut être utilisé dans la compression.\n") ;
    std::exit(EXIT_FAILURE) ;
  }

  else 
  {// C'est un float, évidemment ... ou pas
      if (champs.canbefloat[i]==true)
         {
         if (champs.denom[i]==1)
	    return(FLOAT_CST|mask) ;
         else if (champs.denom[i]==2)
	    return(FLOAT_DENOM_2|mask) ;
         else if (champs.denom[i]>2 && champs.denom[i]<=4)
            return(FLOAT_DENOM_4|mask) ; 
         else if (champs.denom[i]>4 && champs.denom[i]<=16)
	    return(FLOAT_DENOM_16|mask) ; 
         else
	    return(FLOAT|mask) ;
         }
      else
         {
         cout << _("INFO : la compression du champ de type ")<< champs.idx_col[i]<< _("sera malheureusement en double et non en float.\n") ; 
         if (champs.denom[i]==1)
	    return(DOUBLE_CST|mask) ;
         else if (champs.denom[i]==2)
	    return(DOUBLE_DENOM_2|mask) ;
         else if (champs.denom[i]>2 && champs.denom[i]<=4)
            return(DOUBLE_DENOM_4|mask) ; 
         else if (champs.denom[i]>4 && champs.denom[i]<=16)
	    return(DOUBLE_DENOM_16|mask) ; 
         else
	    return(DOUBLE|mask) ;
         }
  }  
return 1 ; 
}

//-------------------------------------------------------
int Compresser::create_datas_denom (vector <unsigned char> & datasout, unsigned char fmt, int idx, vector <double> &datasin)
{
 int longueur, j, reelle=0 ;
 unsigned int i ;
 int bits, decalage ; 
 unsigned char octet, octtmp=0, mask=0 ;
 
 // Initialisation et calcul de la longueur résiduelle
 datasout.clear() ; 
 bits=fmt>>4 ;
 switch(bits)
 {
   case 2 : bits=1 ; mask=1 ; break ; 
   case 3 : bits=2 ; mask=3 ; break ;
   case 4 : bits=4 ; mask=7 ; break ; 
 }
 longueur=ceil(datasin.size()*bits/8.) ;
 
 decalage=0 ; octet=0 ; 
 for (i=0 ; i<datasin.size() ; i++)
   {
   for (j=0 ; j<champs.denom[idx] ; j++)
     {
     if (champs.denombrables[idx][j]==datasin[i])
       {
       octtmp=j ; break ; 
       }
     }
   if (j==champs.denom[idx]) {cout <<_("Une erreur est survenue dans la compression d'un type dénombrable") ; }
   
   octet=octet | (octtmp&mask) ; 
   decalage+=bits ; 
   if (decalage==8)
     {
     datasout.push_back(octet) ; reelle++ ; 
     decalage=0 ; octet=0 ; 
     }
   else
     {octet=octet << bits ; }   
   }
if (decalage>0) {octet=octet << (8-decalage-1) ; reelle++ ; datasout.push_back(octet) ;  }
//printf("|%d ; %d|", longueur, reelle) ; 
return longueur ;  
}

//----------------------
unsigned char Compresser::binary_swap (unsigned char octet, unsigned int bits)
{
 unsigned char retour=0, mask=0 ;
 unsigned int i ;
 
 switch(bits) {
   case 1 : mask=1 ; break ; 
   case 2 : mask=3 ; break ; 
   case 4 : mask=7 ; break ; 
   } 

 for (i=0 ; i<8/bits ; i++)
     {
     retour=retour | (octet & mask) ; 
     if (i<(8/bits-1)) {retour=retour << bits ; octet = octet >> bits ; }
     }  
 return retour ; 
}


bool Compresser::check_all (LucDump &dump) 
{return false ;
}


//==================================
int Compresser::uncompress (FILE * in) 
{
  // Action de décompression des compressés ...
  //LucDump dmp ;
  //ofstream out ;
  /*char datas[500] ;
  unsigned char UC ; 
  unsigned short int USI ; 
  unsigned int UI ; 
  unsigned long int ULI ; 
  double D ; float F ; 
  int I ; long int LI ; 
  int i, j, k, l ;
  */
  cout << _("!!!!!!!!!!!!!!!!!!!!\nCette fonction n'existe qu'à titre de souvenir de la décompression dans la version précédente. Elle ne devrait dans l'idéal plus être utilisée. Il faut construire un objet de classe LcpDump et utiliser la fonction write_dump\n!!!!!!!!!!!!!!\n") ;
  return 1 ;
  
/*  fread(datas, 1, 26, in) ; datas[26]=0 ;  // Extraction de l'header AVFF
  if (strcmp(datas,"AVFFLIGGGHTSDUMPCOMPRESSED")) cout << "WARN : le fichier compressé ne semble pas avir l'header correct\n" ;

  fread(&USI, 2, 1, in) ; 
  fread(datas,1, USI, in) ; datas[USI]=0 ; 
  fread(&dmp.nbsteps, 4, 1, in) ; 
  dmp.steps.resize(dmp.nbsteps) ;

  // Création des colonnes, extraction des formats et dénombrables
  fread(&dmp.steps[0].nb_idx, 1, 1, in) ;
  formats.resize(dmp.steps[0].nb_idx) ; 
  champs.denom.resize(dmp.steps[0].nb_idx) ; 
  champs.denombrables.resize(dmp.steps[0].nb_idx) ;
  for (i=0 ; i<dmp.nbsteps ; i++)
    {
    dmp.steps[i].idx_col.resize(dmp.steps[0].nb_idx) ; 
    dmp.steps[i].nb_idx=dmp.steps[0].nb_idx ; 
    }
  for (i=0 ; i<dmp.steps[0].nb_idx; i++)
    {
    fread(&dmp.steps[0].idx_col[i], 1, 1, in) ; 
    champs.denombrables[i].resize(MAX_DENOM) ; 
    fread(&UC, 1,1, in) ; formats[i]=UC ;
    champs.denom[i]=0 ; 
    if (formats[i]>=17) // C'est un dénombrable !
       {
       unsigned char nombre, format ; 
       nombre=formats[i]>>4 ;
       switch (nombre)
         {
	 case 1 : nombre=1 ; break ;
	 case 2 : nombre=2 ; break ;
	 case 3 : nombre=4 ; break ;
	 case 4 : nombre=16 ; break ; 
	 }
       champs.denom[i]=nombre ; 
       format=formats[i]&15 ; 
       for (j=0 ; j<nombre ; j++)
         {
	 if ((format==1) || (format==2))
	    {fread(&UC, 1, 1, in); champs.denombrables[i][j]=UC ;}     // TODO Séparer signed et unsigned
	 else if ((format ==3) || (format==4))
	    {fread(&USI, 2, 1, in) ; champs.denombrables[i][j]=USI ; }
	 else if ((format==5) || (format==6))
	    {fread(&UI, 4, 1,in)  ; champs.denombrables[i][j]=UI ; }
	 else if (format==9)
            {fread(&F,4,1,in) ; champs.denombrables[i][j]=F ; }
         else if ((format==7) || (format==8))
	    {fread(&LI,8,1,in) ; champs.denombrables[i][j]=LI ; }  
	 else if (format==10)
	    {fread(&D,8,1,in) ; champs.denombrables[i][j]=D ; }  
         else
            {cout << "Format inconnu dans la sous fonction convertDBL" ; std::exit(EXIT_FAILURE) ;}
	 }
       }
    }   
  //champs.idx_col.resize(dmp.steps[0].nb_idx) ; 
  //champs.disp() ; 
  
  // Ecriture des caractéristiques des données à tous les pas de temps
  for (i=1 ; i<dmp.nbsteps ; i++)
    {
    dmp.steps[i].idx_col=dmp.steps[0].idx_col ;  
    }
    
  // Lecture complète du fichier, et écriture du fichier décompressé en direct (sinon on va avoir des problèmes de mémoire pour les très gros fichiers.
  out.open(strcat(datas,".uc"), ios::binary | ios::out) ; 
  actions.total=dmp.nbsteps ; 
  actions.disp_progress() ; 
  for (i=0 ; i<dmp.nbsteps ; i++)
    {
    actions.valeur=i ; 
    fread(&dmp.steps[i].timestep,4,1,in) ; 
    fread(&dmp.steps[i].nb_atomes,4,1,in) ;  
    fread(&F,4,1,in) ; dmp.steps[i].box[0][0]=F; fread(&F,4,1,in) ; dmp.steps[i].box[0][1]=F;
    fread(&F,4,1,in) ; dmp.steps[i].box[1][0]=F; fread(&F,4,1,in) ; dmp.steps[i].box[1][1]=F; 
    fread(&F,4,1,in) ; dmp.steps[i].box[2][0]=F; fread(&F,4,1,in) ; dmp.steps[i].box[2][1]=F; 
    dmp.steps[i].datas.resize(dmp.steps[i].nb_idx) ;
    for (j=0 ; j<dmp.steps[i].nb_idx ; j++)
      {
      dmp.steps[i].datas[j].resize(dmp.steps[i].nb_atomes) ; 
      if (formats[j]<17) // Ce n'est pas un dénombrable
	 {
	 unsigned char format ; 
	 format=(formats[j])&15 ; 
         for (k=0 ; k<dmp.steps[i].nb_atomes ; k++)
	   {
	   if ((format==1) || (format==2))
	    {fread(&UC, 1, 1, in); dmp.steps[i].datas[j][k]=UC ;}     // TODO Séparer signed et unsigned
	   else if ((format==3) || (format==4))
	    {fread(&USI, 2, 1, in) ; dmp.steps[i].datas[j][k]=USI ; }
	   else if ((format==5) || (format==6))
	    {fread(&UI, 4, 1,in)  ; dmp.steps[i].datas[j][k]=UI ; }
	   else if (format==9)
            {fread(&F,4,1,in) ; dmp.steps[i].datas[j][k]=F ; }
           else if ((format==7) || (format==8))
	    {fread(&LI,8,1,in) ; dmp.steps[i].datas[j][k]=LI ; }  
	   else if (format==10)
	    {fread(&D,8,1,in) ; dmp.steps[i].datas[j][k]=D ; }  
           else
            {cout << "Format inconnu dans la sous fonction convertDBL" ; std::exit(EXIT_FAILURE) ;}
	   }
	 }
      else if (formats[j]>=CHAR_CST && formats[j]<CHAR_DENOM_2)
         {for (k=0 ; k<dmp.steps[i].nb_atomes ; k++){dmp.steps[i].datas[j][k]=champs.denombrables[j][0] ;}} 
      else if (formats[j]>=CHAR_DENOM_2 && formats[j]<=DOUBLE_DENOM_16)
         {
	 unsigned char bits, mask, octet ; int longueur ; 
	 bits=(formats[j]>>4)&15 ;
	 switch(bits)
	     {
	     case 2 : bits=1 ; mask=1 ; break ; 
	     case 3 : bits=2 ; mask=3 ; break ;
	     case 4 : bits=4 ; mask=7 ; break ; 
	     }
         longueur=ceil(dmp.steps[i].nb_atomes*bits/8.) ;
	 for (k=0 ; k<longueur; k++)
	     {
	     fread(&octet, 1,1, in) ; octet=binary_swap(octet, bits) ; 
	     for (l=0 ; l<(8/bits) ; l++)
	       {dmp.steps[i].datas[j][k*(8/bits)]=champs.denombrables[j][octet&mask] ; octet=octet>>bits ;  }
	     }     
	 }
      else
         {
	 cout << "ERR : format inconnu lors de la décompression des steps." ;
	 }
      }
    
    // On a fini de lire toutes les données du timestep je pense ...
    // On écrit le fichier de sortie et on vide le step.
    dmp.steps[i].write_asDUMP(out) ; 
    dmp.free_timestep(i) ; 
    }
out.close() ;   */
}


