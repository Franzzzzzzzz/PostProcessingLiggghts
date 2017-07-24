#include <cstdlib>
#include <cstdio>
#include <map>

using namespace std ; 

int main ()
{
 map <int, int> hello ; 


 for (int i=10 ; i>0 ; i--)
   hello[i]+=i ; 
 
 for (int i=1 ; i<5 ; i++)
   printf("%d %d \n",i, hello[i]) ;  
 
}
