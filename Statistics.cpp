#include "Headers/Statistics.h"
//=========================================

int Stats::compute_global (Dump & dump)
{
 int i, j ; 
 bool always_the_same = true ; 
 vector <vector <double> > temp_datas ; 
 int compteur=0 ;  
  // Le downsampling et l'extract ne jouent pas ici

 dump.check_timestep(0) ;
 //idx_col.resize(dump.steps[0].nb_idx, UNKNOWN) ;
 idx_col=dump.steps[0].idx_col ; 

 min.resize(dump.steps[0].nb_idx, DBL_MAX) ; 
 max.resize(dump.steps[0].nb_idx, DBL_MIN) ; 
 same.resize(dump.steps[0].nb_idx, true) ;

 canbefloat.resize(dump.steps[0].nb_idx, true) ;  
 denom.resize(dump.steps[0].nb_idx,0) ; 
 denombrables.resize(dump.steps[0].nb_idx) ; 
 for (i=0 ; i<dump.steps[0].nb_idx; i++)
 {denombrables[i].resize(MAX_DENOM, 0) ; }
 
 actions.total=dump.nbsteps ;
 cout << "\nStats::compute_global          " ;
 actions.disp_progress() ; 

 temp_datas.resize(dump.steps[0].nb_idx, vector <double>(1,0)) ; 
 for (i=0 ; i<dump.nbsteps ; i++)
   {
   actions.valeur=i ; 
   dump.check_timestep(i) ; 
   compute(dump.steps[i]) ; 
   if (always_the_same && i>0)
      {
      // Comparaison par colonne
      for (j=0 ; j<dump.steps[i].nb_idx ; j++)
          {if (temp_datas[j] != dump.steps[i].datas[j]) same[j]=false ; }
      // Si toutes les colonnes sont fausses on ne s'embête plus ...
      compteur=0 ; 
      for (j=0 ; j<dump.steps[i].nb_idx ; j++)
          {if (same[j]==false) compteur++ ; }
      if (compteur==dump.steps[i].nb_idx) always_the_same=false ; 
      } 
   if (always_the_same)
      {
      // recopie des données
      for (j=0 ; j<dump.steps[i].nb_idx ; j++)
          {temp_datas[j]= dump.steps[i].datas[j];}
      }
  }
return 1 ;
}

//-------------------------------------
int Stats::compute_step (Dump & dump, int ts)
{
int i ; 
dump.check_timestep(ts) ;
idx_col=dump.steps[ts].idx_col ; 
min.resize(dump.steps[ts].nb_idx, DBL_MAX) ; 
max.resize(dump.steps[ts].nb_idx, DBL_MIN) ; 
denom.resize(dump.steps[ts].nb_idx,0) ; 
denombrables.resize(dump.steps[ts].nb_idx) ; 
for (i=0 ; i<dump.steps[ts].nb_idx; i++)
{denombrables[i].resize(MAX_DENOM, 0) ; }

dump.check_timestep(ts) ; 
compute(dump.steps[ts]) ; 
return 1 ;
}

//========================================
int Stats::compute(Step &step)
{
 unsigned int i,j,k ;
 double dbl ; float flt ; 
 
 for (k=0 ; k<denom.size() ; k++)
   {
   //cout << index << "   " ; fflush(stdout) ; 
   //idx=step.find_idx(k) ; 
 
   for (i=0 ; i<step.datas[k].size() ; i++)
     {
     max[k]=step.datas[k][i]>max[k]?step.datas[k][i]:max[k] ; 
     min[k]=step.datas[k][i]<min[k]?step.datas[k][i]:min[k] ;
     if (canbefloat[k]==true)
        {
        flt=(float) (step.datas[k][i]) ; dbl=(double) flt ; 
        if ((dbl-step.datas[k][i])>=0.000001 || (dbl-step.datas[k][i])<=-0.000001) canbefloat[k]=false ; 
        }
     //cout << denom[k] <<" "; 
     if (denom[k]<=MAX_DENOM)
       {
       for (j=0 ; j<denom[k] ; j++)
         {
         if (denombrables[k][j]==step.datas[k][i]) 
	    {j=MAX_DENOM+10 ; }
         }
         
       if (j<MAX_DENOM+10)
         {
         denom[k]++ ;
         if (denom[k]<=MAX_DENOM)
            denombrables[k][denom[k]-1]=step.datas[k][i] ; 
         }
     
       }   
     }
   }
return 1 ;
}


//=======================================
bool Stats::check_index(Dump &dump)
{
 int i,j ; 
 
 for (i=1 ; i<dump.nbsteps ; i++)
   {
   for (j=0 ; j<dump.steps[i].nb_idx ; j++)
     {
     if (dump.steps[i].idx_col[j]!=dump.steps[0].idx_col[j])
       return false ; 
     }
   }
  
  return true ;   
}

//=================================================
double Stats::minmaxdst(Dump &dump, int ts )
{
int i , j ;
double r ;  
double rmin=10000000, rmax=0 ; 
int idx[3] ; 

dump.check_timestep(ts) ;
idx[0]=dump.steps[ts].find_idx(IDS("POSX")) ;
idx[1]=dump.steps[ts].find_idx(IDS("POSY")) ; 
idx[2]=dump.steps[ts].find_idx(IDS("POSZ")) ;   
printf("%lf %d %d\n", dump.steps[ts].datas[5][0], idx[1], idx[2]) ; 
for (i=0 ; i<dump.steps[ts].nb_atomes ; i++)
  {
  printf("%d ", i) ; fflush(stdout) ; 
  for (j=0 ; j<i ; j++)
      {
      //printf("[%d]", j) ; fflush(stdout) ; 
      r=sqrt((dump.steps[ts].datas[idx[0]][i]-dump.steps[ts].datas[idx[0]][j])*(dump.steps[ts].datas[idx[0]][i]-dump.steps[ts].datas[idx[0]][j])+
             (dump.steps[ts].datas[idx[1]][i]-dump.steps[ts].datas[idx[1]][j])*(dump.steps[ts].datas[idx[1]][i]-dump.steps[ts].datas[idx[1]][j])+
	     (dump.steps[ts].datas[idx[2]][i]-dump.steps[ts].datas[idx[2]][j])*(dump.steps[ts].datas[idx[2]][i]-dump.steps[ts].datas[idx[2]][j])) ;
      //printf("[%lf]", r) ; fflush(stdout) ; 
      if (r<rmin) rmin=r ; 
      if (r>rmax) rmax=r ;
      }
  }

printf(_("Distance min entre atomes : %f\nDistance max entre atomes : %f\n"), rmin, rmax) ;  
return 0 ; 
}




//======================================
int Stats::disp ()
{
 unsigned int i ;
 cout << _("\n---------------------------\n Statistiques des données : \n---------------------------\n IDX  ") ;
 
 for (i=0 ; i<idx_col.size(); i++)
   {printf("%d ", idx_col[i]) ; }
 cout<<"\nDENOM " ; 
 for (i=0 ; i<idx_col.size(); i++)
   {printf("%d ",denom[i]) ; }
 cout<<"\nMIN   " ; 
 for (i=0 ; i<idx_col.size(); i++)
   {printf("%f ", min[i]) ; }
 cout<<"\nMAX   " ; 
 for (i=0 ; i<idx_col.size(); i++)
   {printf("%f ", max[i]) ; }
/* cout << "\n" ; 
 for (i=0 ; i<idx_col.size() ; i++)
   { 
   for (j=0 ; j<denom[i]-1 ; j++)
     {cout << denombrables[i][j] << "|" ; }
   cout << "\n" ;   
   }*/
 cout <<"\n\n" ; 
return 1 ;
}



