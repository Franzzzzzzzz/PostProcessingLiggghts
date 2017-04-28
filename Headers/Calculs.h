#ifndef MASTERH
#include "MasterHeader.h"

#else 
class Calcul ; 

//=================================================================
class Vector ;
class Matrix3x3 {
public :
  Matrix3x3() ;
  Matrix3x3(const char * str) ; 
  // Operateurs
    //Recherche
    double & operator() (int i, int j) {return dat[(i-1)*3+(j-1)] ;}
    double & operator[] (int i) {return dat[i] ;}
    //affectation
    //void operator= (Matrix3x3 & m) ; 
    void operator= (double x) {for (i=0 ; i<9 ; i++) dat[i]=x ;}
    // Addition
    Matrix3x3  operator+ (double x ) ;
    Matrix3x3  operator+ (Matrix3x3 m) ;
    // Soustraction et opposé
    Matrix3x3  operator- () ;
    Matrix3x3  operator- (double x) {return ((*this)+(-x));};    
    Matrix3x3  operator- (Matrix3x3 & m) {Matrix3x3 tmp=-m ; return ((*this)+(tmp)) ;}
    // Multiplication
    Matrix3x3  operator* (double x) ;
    Matrix3x3  operator* (Matrix3x3 & m) ;
    //Vector operator* (Vector & v) ;
    // Division
    Matrix3x3  operator/ (double x) ;    
  //Fonctions
  Matrix3x3  t (void) ; // transpose
  void disp (void) ; 
  
  bool tmp ;  
private :
  double dat[9] ; 
  int i ; 
} ;
//------------------------------------------------------------
class Vector {
public :
  // Constructeurs
  Vector() {dat[0]=dat[1]=dat[2]=0 ; iscolumn=true ; } ;
  Vector(double a)   {dat[0]=dat[1]=dat[2]=a ; iscolumn=true ; } ;
  Vector(const Vector & v)   {dat[0]=v.dat[0] ; dat[1]=v.dat[1] ; dat[2]=v.dat[2] ; iscolumn=v.iscolumn ; } ;
  Vector(double x, double y, double z) {dat[0]=x ; dat[1]=y ; dat[2]=z ; iscolumn=true ; } ;

  //Recherche
  double & operator() (int i) {return dat[i-1] ;}
  double & operator[] (int i) {return dat[i] ;}

  //affectations
  Vector & operator= (const Vector & m) {dat[0]=m.dat[0] ; dat[1]=m.dat[1] ; dat[2]=m.dat[2] ; iscolumn=m.iscolumn ; return *this ; };
  Vector & operator= (double x) {dat[0]=dat[1]=dat[2]=x ; iscolumn=true ; return *this ; };

  // Arithmétique + affectation (l'arithmétique de base est en dehors de la classe)
  // NB : la multiplication & la division sont par une constante
  // 	  Le modulo est remplacé par la multiplication terme à terme (.*)
  Vector & operator+= (const Vector & m) {dat[0]+=m.dat[0] ; dat[1]+=m.dat[1] ; dat[2]+=m.dat[2] ; return *this ; }
  Vector & operator-= (const Vector & m) {dat[0]-=m.dat[0] ; dat[1]-=m.dat[1] ; dat[2]-=m.dat[2] ; return *this ; }
  Vector & operator*= (const Vector & m) {dat[0]*=m.dat[0] ; dat[1]*=m.dat[1] ; dat[2]*=m.dat[2] ; return *this ; }
  Vector & operator/= (const Vector & m) {dat[0]/=m.dat[0] ; dat[1]/=m.dat[1] ; dat[2]/=m.dat[2] ; return *this ; }
  Vector & operator%= (const Vector & m) {dat[0]*=m.dat[0] ; dat[1]*=m.dat[1] ; dat[2]*=m.dat[2] ; return *this ; }
  // + Autres fonctions arithmétiques ...
  Vector div (void) { Vector c(*this) ; c.dat[0]=1/dat[0] ; c.dat[1]=1/dat[1] ; c.dat[2]=1/dat[2] ; return c ;} 			// Inverse
  Vector div (const Vector & m) {Vector c(*this) ; if (m.iscol()!=iscolumn && !m.isconst()) {printf("VECERR:WAYdiv") ;return 0 ;} ; c/=m ; return c ;}	// Division terme à terme ./

  // Comparaison & Test
  bool operator< (Vector const & v) {if (dat[0]<v.dat[0] && dat[1]<v.dat[1] && dat[2]<v.dat[2]) return true ; return false ;}
  bool operator> (Vector const & v) {if (dat[0]>v.dat[0] && dat[1]>v.dat[1] && dat[2]>v.dat[2]) return true ; return false ;}
  bool operator<= (Vector const & v) {if (dat[0]<=v.dat[0] && dat[1]<=v.dat[1] && dat[2]<=v.dat[2]) return true ; return false ;}
  bool operator>= (Vector const & v) {if (dat[0]>=v.dat[0] && dat[1]>=v.dat[1] && dat[2]>=v.dat[2]) return true ; return false ;}
  bool operator== (Vector const & v) {if (dat[0]==v.dat[0] && dat[1]==v.dat[1] && dat[2]==v.dat[2] && iscolumn==v.iscolumn) return true ; return false ;}
  bool operator!= (Vector const & v) {if (!(*this==v)) return true ; return false ;}
  bool isconst (void) const {if (dat[0]==dat[1] && dat[1]==dat[2]) return true ; return false ;}
  bool iscol (void) const { return iscolumn ; }


  // Orientation
  Vector  t (void) {Vector res(*this) ; res.autochange() ; return res ; } // transpose
  void autochange (void) {iscolumn=!iscolumn;}
  void set (double a, double b, double c) {dat[0]=a ; dat[1]=b ; dat[2]=c ;} 
  bool isnan (void) {if (isnan(dat[0]) || isnan(dat[1]) || isnan(dat[2])) return true ; return false ;}
  bool isnan (double & a) {return (std::isnan(a)) ;}
  
  // Fonctions vectorielles
  double dot(Vector v) {return (v[0]*dat[0]+v[1]*dat[1]+v[2]*dat[2]) ;}
  Vector cross(Vector w) {Vector res ; res.dat[0]=dat[1]*w.dat[2]-dat[2]*w.dat[1] ;
						     res.dat[1]=dat[2]*w.dat[0]-dat[0]*w.dat[2] ;
  	  	  	  	  	  	     res.dat[2]=dat[0]*w.dat[1]-dat[1]*w.dat[0] ; 
						     res.iscolumn=iscolumn ; return res ; }
  double norm () { return (sqrt(dat[0]*dat[0]+dat[1]*dat[1]+dat[2]*dat[2])) ; }
  void normalise () { double d=norm() ; dat[0]/=d ; dat[1]/=d ; dat[2]/=d ; }
  Matrix3x3 multiply (Vector & v) ;

  void disp (void) ;

  double dat[3] ;  
  bool iscolumn ;
};

inline Vector operator+ (Vector const& a, Vector const& b) {Vector c(a) ; if (a.iscol()!=b.iscol()) {printf("VECERR:WAY+") ; c.dat[0]=c.dat[1]=c.dat[2]=NAN ; return c ;} c+=b ; return c ; }
inline Vector operator- (Vector const& a) {Vector c(a) ; c-=a ; c-=a ; return c ; }
inline Vector operator- (Vector const& a, Vector const& b) {Vector c(a) ; if (a.iscol()!=b.iscol()) {printf("VECERR:WAY-") ; c.dat[0]=c.dat[1]=c.dat[2]=NAN ; return c ;} c-=b ; return c ; }
inline Vector operator* (Vector const& a, Vector const& b) {Vector c(a) ; if (!(a.isconst()||b.isconst())) {/*printf("VECERR:WAY*") ;*/ c.dat[0]=c.dat[1]=c.dat[2]=NAN ; return c ;} c*=b ; return c ; }
inline Vector operator* (Matrix3x3 & m, Vector & v) ;
inline Vector operator* (Vector & v, Matrix3x3 & m) ;
inline Vector operator/ (Vector const& a, Vector const& b) {Vector c(a) ; if (!(b.isconst())) {printf("VECERR:WAY/") ; c.dat[0]=c.dat[1]=c.dat[2]=NAN ; return c ;} c/=b ; return c ; }
inline Vector operator% (Vector const& a, Vector const& b) {Vector c(a) ; if (a.iscol()!=b.iscol()) {printf("VECERR:WAY%%") ; c.dat[0]=c.dat[1]=c.dat[2]=NAN ; return c ;} c%=b ; return c ; }
inline Vector cross(Vector & a, Vector & b) {return (a.cross(b)) ; }
inline double triple_product (Vector a, Vector b, Vector c) {return (a.dot(b.cross(c))) ; }
inline double mixed_product (Vector a, Vector b, Vector c) {return (triple_product(a,b,c)) ; } // alias for triple_product 

//Surcharge des operateurs sur réels
inline Matrix3x3  operator+ (double x, Matrix3x3 & m) {return (m+x) ; }
inline Matrix3x3  operator- (double x, Matrix3x3 & m) {return ((-m)+x) ; }
inline Matrix3x3  operator* (double x, Matrix3x3 & m) {return (m*x) ; }
//inline void operator/ (double x, Matrix3x3 m) {std::cerr << "MATRIXERR : double/matrix" ; }
/*inline Vector  operator+ (double x, Vector & m) {return (m+x) ; }
inline Vector  operator- (double x, Vector & m) {return ((-m)+x) ; }
inline Vector  operator* (double x, Vector & m) {return (m*x) ; }*/
//inline void operator/ (double x, Vector m) {cerr << "MATRIXERR : double/vector\n" ; }
//============================================================================

typedef struct {
  Vector axe ; 
  Vector centre ; 
  double D, L ; 
} Cylindre ; 

class Icosahedre {
public:
  Icosahedre() ; 
  int write_vtk (string chemin) ;
  int belonging_tri (Vector V) ;
  int deform () ; 
  int subdivide() ;
  int subdivide(int n) {if (n>1) subdivide(n-1) ; subdivide() ; }
  double *data ; 
  int npts, nedges, nfaces ; 
  double *solidangle ; 
  
private :
  double **pts ;
  int **edges; 
  int **faces ;
} ; 

class Geometrie
{
public :
//static Matrix gravicentre (Matrix pt1, Matrix pt2, Matrix pt3) ;
static Vector gravicentre (Vector pt1, Vector pt2, Vector pt3) ;
// With newmat
//static Matrix rotation (Matrix vecteur, double angle, Matrix axe) ; 
//static Matrix rotation (Matrix vecteur, Matrix rot) ; 
//static Matrix rotation_axeZ (Matrix vecteur, double angle) ;
//static Matrix get_rot_matrix (double angle, Matrix axe) ; 
//static Matrix cart2sph (Matrix cart) ; 
static Vector cart2sph (Vector cart) ; 
static Vector sph2cart (Vector sph) ; 

// Without newmat
static Matrix3x3 rotation (Matrix3x3 mat, Matrix3x3 rotation) ;
static	  Vector rotation (Vector vec, Matrix3x3 rotation) ;
static Matrix3x3 rotation (Matrix3x3 mat, double angle, Vector axe) ;
static    Vector rotation (Vector vecteur, double angle, Vector axe) ;
static Matrix3x3 rotation_axeZ (Matrix3x3 vecteur, double angle) ; 
static    Vector rotation_axeZ (Vector vec, double angle) ;
static Matrix3x3 get_rot_matrix (double angle, Vector axe) ; 

static double distance_cylindre (Cylindre & C, Vector x) ; 
static Vector contact_cylindre (Cylindre & C, Vector x) ;
static double tri_surface (Vector a, Vector b, Vector c) ;
static double solid_angle (Vector a, Vector b, Vector c) ; 
} ; 

/*class Convert
{
public :
static Matrix toVect (double a, double b, double c) {Matrix t(3,1) ; t<<a<<b<<c ; return t ; }
static Matrix toVect (double a, double b) {Matrix t(2,1) ; t<<a<<b ; return t ; }
} ;*/ 

class Calcul
{
public :
static Vector cross_product (Vector a, Vector b) ; 
static double norm (Vector x) {return sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3));}
static double norm2D (Vector x) {return sqrt(x(1)*x(1)+x(2)*x(2));}
static double arctan (double x, double y) { return (x==0?(y>0?(M_PI/2):(-M_PI/2)):(x>0?atan(y/x):atan(y/x)+M_PI)) ; }
static double sgn (double x) { return x<0?-1:(x==0?0:1) ;}
static double angle_0_2pi (double angle) ;

static double true_random () ; 
static double fborne (double v, double a, double b) {double tmp; if (b<a) {tmp=a ; a=b ; b=tmp ;} if (v<=a) return a ; else if (v>=b) return b ; else return v ; }

static double blockstd (double * x, int N) ;
static double mean(double *x, int N) ; 
static double var(double *x, int N) ; 

static int eigen (Matrix3d A, Vector3d & eigval, Matrix3d & eigvec) ;


};

class Fonction
{
public:
static void gaussienne1D (double & res,const double & x, const double & sigma) ;
static void creneau1D (double & res,const double & x, const double & sigma) ;

static void gaussienne3D (double & res,const double & x, const double & w, double radius) ; 
static void gaussienne3D_multivariee (double & res,  Vector & v, Vector & w, double radius) ; 
static void creneau3D (double & res,  Vector & x,  Vector & w, double radius) ; 
static void creneau3D_fin (double & res, Vector x, Vector w, double radius) ;
static void intersect_sphere (double & res, Vector &x, Vector &w, double radius) ;
static void intersect_cercle (double & res, Vector &x, Vector &w, double radius) ;
static void intersect_creneau2D (double & res, Vector &x, Vector &w, double radius) ;
static void intersect_tranche(double &res, Vector &x, Vector &w, double radius, char dir) ; 
static void int_creneau3D (double & res,  Vector & xa,  Vector & xb,   Vector & w) ;
static void int_creneau2D (double & res, Vector &xa, Vector &xb, Vector &w) ;
static void int_gaussienne3D(double & res, Vector & xa,  Vector & xb,  Vector & w) ; 
static void int_intersect_sphere (double & res, Vector &xa, Vector &xb, Vector & w) ;
static double volume_intersect_sphere_cylinder (Vector xsphere, double Rsphere, Vector xcyl, Vector dircyl, double Rcyl, double Lcyl) ;
static double volume_intersect_pave_cylindre(Vector pos, Vector sigmafenetre, double R) ;
static double aireleft_intersect_cercle (double R, double r, double x) ;
// Overload
static void gaussienne3D (double & res,  Vector & x,  Vector & w, double radius) ;
static void intersect_tranche_x(double &res, Vector &x, Vector &w, double radius) {Fonction::intersect_tranche(res, x, w, radius, 'x') ; }
static void intersect_tranche_y(double &res, Vector &x, Vector &w, double radius) {Fonction::intersect_tranche(res, x, w, radius, 'y') ; }
static void intersect_tranche_z(double &res, Vector &x, Vector &w, double radius) {Fonction::intersect_tranche(res, x, w, radius, 'z') ; }
static int  intersect3D_RayTriangle(Vector ray, Vector t1, Vector t2, Vector t3) ; 
} ; 

#endif






