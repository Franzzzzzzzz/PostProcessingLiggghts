#define MASTERH
#include <math.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <pthread.h>
#include <unistd.h>
#include <stdarg.h>
#include <algorithm>
#include <assert.h>
#include <zlib.h>
#include <ctime>
//#include <limits>

#ifdef USETIFF
#include <tiffio.h>
#endif
//#include "../newmat10/newmat.h"

#ifdef MATLAB 
#include <mat.h>
#include <matrix.h>
#endif

#ifdef BOOST
#define MCMAXITER 100000
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#endif

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Unsupported/KroneckerProduct"

using namespace std ; 
using namespace Eigen ;

class Step ; 
class Dump ; 

// Fonctions sympatiques pour de la couleur ...
#define DISP_Info(param) printf("\033[32m%s\033[0m",param)
#define DISP_Warn(param) printf("\033[33m%s\033[0m",param)
#define DISP_Err(param) printf("\033[31m%s\033[0m",param)

//#define Debug(a) printf(a);fflush(stdout)
#include "Gunzip.h"
#include "Ids.h"
#include "Calculs.h"
#include "Multisphere.h"
#include "Step.h"
#include "Dump.h" // Inclus aussi Compress.h et Stat.h
#include "FDump.h"
#include "LDump.h"
#include "LcfDump.h"
#include "Surface.h"
#include "Filter.h"
#include "Coarse.h"
#include "Writing.h"
#include "Main.h"

//#include "Ssh.h"
