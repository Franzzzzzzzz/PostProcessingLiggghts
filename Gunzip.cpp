#include "Headers/Gunzip.h"

void gzifstream::open (const char* filename)
{
  int l ; 
  l=strlen(filename) ; 
  strcpy(name, filename) ; 
  if (filename[l-3]=='.' && filename[l-2]=='g' && filename[l-1]=='z') 
  {
    isgz=true ; 
    if (!actions["dumpall"].set) DISP_Warn("Il vaut mieux utiliser l'option dumpall vec les dump compressés gz, si la mémoire le permet.\n") ; 
    gzopen(filename) ;
  }
  else 
  {
    isgz=false ;
    in.open(filename) ; 
  }
}
//============================================
int gzifstream::gzopen(const char * filename)
{
  int ret ; 
  isgz=true ; 
  nomoredata=false ; 
  cur=0 ; qttavail=0 ; curtot=0 ;
  
  gzin=fopen(filename, "rb") ; 
  if (gzin==NULL) DISP_Err("Impossible d'ouvrir le dump\n") ; 
  
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  strm.avail_in = 0;
  //strm.next_in = Z_NULL;
  strm.next_in = inbuffer ; 
 
  /* allocate inflate state */
  ret = inflateInit2(&strm, 16+MAX_WBITS); 
  return ret;
}

//============================================
unsigned long int gzifstream::size()
{
 unsigned long int res ; unsigned long int  oldpos ; 
 if (isgz)
 {
   oldpos=(unsigned long int)(ftell(gzin)) ; 
   fseek(gzin, 0L, SEEK_END);
   res = (unsigned long int)(ftell(gzin));
   res=round(res/0.2) ; 
   fseek(gzin, oldpos, SEEK_SET);
 }
 else
 {
  oldpos=(unsigned long int)(in.tellg()) ; 
  in.seekg(0,ios::end) ; 
  res=(unsigned long int)(in.tellg()) ; 
  in.seekg(oldpos, ios::beg) ; 
 }
 return res ; 
}

//============================================
void gzifstream::gzgetline (unsigned char** s, streamsize n )
{
  long int i ; 
  for (i=cur ; i<qttavail && mybuffer[i]!='\n' && i-cur<n-1  ; i++) ;
  
  if (mybuffer[i]=='\n')
  {
    mybuffer[i]=0 ; 
    *s=mybuffer+cur ; 
    curtot+=(i-cur+1) ; 
    cur=i+1 ;
  }
  else if (i==qttavail) 
  { //system("sleep 0.2") ; 
    if (!nomoredata)
    {
     getmoredata() ; 
     gzgetline(s, n) ; 
    }
    else
    {
     mybuffer[i]=0 ;  // WARNING risque de pb si le fichier se termine sur la bordure de CLUSTER
     *s=mybuffer+cur ; 
     curtot+=(i-cur+1) ; 
     cur=i+1 ; 
    }
  }
  else if (i-cur==n-1)
    DISP_Warn("Ligne trop longue ... Bizarre \n") ; 
  else
    DISP_Err("Ne peut pas arriver ...\n") ; 
}
//---------------------------------------------
void gzifstream::getline (string & ligne) 
{
  if (isgz)
  {
    unsigned char *l ; 
    gzgetline(&l, 5000) ; 
    ligne.assign((char *)l) ; 
  }
  else std::getline(in, ligne) ;  
}
//---------------------------------------------
unsigned char gzifstream::getnextbyte()
{
 if (cur<qttavail) 
 {
   curtot++ ; 
   return (mybuffer[cur++]) ;
 }
 else
 {
  if (!nomoredata)
  {
    getmoredata() ; 
    return (getnextbyte()) ; 
  }
  else
    return (0) ; // Won't happen hopefully.
 }
}
//---------------------------------------------
void gzifstream::getmot(char * mot) 
{
  bool hasnonspace=false ; int i ; 
  i=0 ;
  do
  {
   mot[i]=(char)getnextbyte() ; 
   if (hasnonspace==false && !isspace(mot[i])) hasnonspace=true ;
   
   i++ ;    
  } while (hasnonspace==false || (hasnonspace==true && !isspace(mot[i-1])) ) ; 
  mot[i]=0 ; 
  
}

//---------------------------------------------
gzifstream& gzifstream::operator>> (int& val) 
{
  if (isgz)
  {
   char mot[500] ; 
   getmot(mot) ;
   val=atoi(mot) ; 
  }
  else 
    in>>val ;
  
  return (*this) ; 
}
gzifstream& gzifstream::operator>> (float& val)
{
  if (isgz) 
  {
   char mot[500] ; 
   getmot(mot) ; 
   val=(float)atof(mot) ; 
  }
  else 
    in>>val ; 
  
  return (*this) ; 
  
}
gzifstream& gzifstream::operator>> (double& val)
{
  if (isgz) 
  {
   char mot[500] ; 
   getmot(mot) ; 
   val=atof(mot) ; 
  }
  else 
    in>>val ;
  
  return (*this) ;  
}
  
//===============================================
gzifstream& gzifstream::seekg (unsigned long int pos) 
{
  if (isgz)
    gzmove(pos) ;
  else
    in.seekg(pos) ;    
  return (*this)  ; 
}
//------------------------------------------------
gzifstream& gzifstream::seekg (long int off, ios_base::seekdir way)
{
  if (isgz) 
  {
    if (way==ios::end) DISP_Err("Non implémenté: impossible de se déplacer de la fin en gz\n") ;
    else if (way==ios::cur) gzmove(curtot+off) ; 
    else if (way==ios::beg) gzmove(off) ;
    else DISP_Err("Should not happen (bad way)") ; 
  }
  else 
  {
    in.seekg(off, way) ; 
  }   
  return (*this)  ; 
}
//------------------------------------------------
void gzifstream::gzmove (unsigned long int pos)
{
  if (pos<curtot)
  {
    if (cur>curtot-pos) 
    {
      //printf("%ld %ld %ld\n", cur, curtot, pos) ; 
      while (curtot>pos)
      {
	cur-- ; curtot-- ; 
	if (mybuffer[cur]==0) mybuffer[cur]='\n' ; 
      }
    }
    else
    {
     // printf("%ld %ld %ld\n", cur, curtot, pos) ; 
    DISP_Warn("Déplacement en arrière, pas assez de données dans le buffer: relecture du fichier (il faudrait éviter ça ...)\n") ; 
    fclose(gzin) ; 
    gzopen(name) ; 
    }
  }
  
  pos-=curtot ; 
  
  if (cur==qttavail) 
  {
    if (!nomoredata) getmoredata() ;
    else {printf("%ld %ld %ld\n", cur, curtot, pos) ; DISP_Err("Mouvement au delà de la fin du fichier") ; exit(1) ; }
  }
    
  while (pos>0)
  {
    cur++ ; pos-=1 ; curtot++ ; 
    if (cur==qttavail) 
    {
      if (!nomoredata) getmoredata() ;
      else {DISP_Err("Mouvement au delà de la fin du fichier") ; exit(1) ; }
    }
  }
}

/*int gzifstream::getdatabackward()
{
  
}*/

//============================================
int gzifstream::getmoredata()
{
  unsigned long int i ; 

  if (nomoredata) {DISP_Err("Impossible de décompresser plus de données\n") ; return 1; }
  
  // On recopie ce qui n'a pas encore été lu dans le début du buffer.
  for (i=cur ; i<qttavail ; i++) mybuffer[i-cur]=mybuffer[i] ; 
  qttavail-=cur ;
  cur=0 ; 
  
  unzip(CLUSTER-qttavail) ; 
  return 0 ; 
}


//==============================================
int gzifstream::unzip(unsigned long int qttwanted)
{
 int ret=Z_OK;
 unsigned have;

 unsigned long int reste=qttwanted ; 
 
/*  printf("\nNEW{") ; 
 for (int j=0 ; j<100 ; j++) printf("%c", mybuffer[j]) ; 
 printf("}\n") ;*/
 
 strm.avail_out=qttwanted ; 
 //printf("[%d %d]", qttwanted, qttavail) ; fflush(stdout) ; 
   
 while (qttwanted>0 && ret!=Z_STREAM_END)
 {
   strm.next_out=mybuffer+qttavail ; 
   ret = inflate(&strm, Z_NO_FLUSH);
   qttavail=qttavail+qttwanted-strm.avail_out ; 
   qttwanted=strm.avail_out ;
   //printf("[%d %d]", qttwanted, qttavail) ; fflush(stdout) ; 
   assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
   switch (ret) 
   {
     case Z_NEED_DICT: ret = Z_DATA_ERROR;     /* and fall through */
     case Z_DATA_ERROR:
     case Z_MEM_ERROR: (void)inflateEnd(&strm); DISP_Err("gz Erreur occured\n") ; return ret;
   }
   
   if (strm.avail_out>0 && ret!=Z_STREAM_END) // pas assez de data dans l'input pour décompresser ce qu'on a demandé
   {
     strm.avail_in = fread(inbuffer, 1, CLUSTER, gzin);
     if (ferror(gzin)) {(void)inflateEnd(&strm); DISP_Err("gz: Error reading file") ; return Z_ERRNO;}
     strm.next_in = inbuffer ;
   }
 }
//  printf("\n{") ; 
// for (int j=0 ; j<100 ; j++) printf("%c", mybuffer[j]) ; 
// printf("}\n") ;
 if (ret==Z_STREAM_END) {(void)inflateEnd(&strm); nomoredata=true ; }
 return ret;
}

