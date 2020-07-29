#include "Headers/Calculs.h"

/*Matrix Geometrie::cart2sph (Matrix cart)
{
Matrix sph(3,1) ;
sph(1,1)=sqrt(cart(1,1)*cart(1,1)+cart(2,1)*cart(2,1)+cart(3,1)*cart(3,1)) ;
sph(2,1)=acos(cart(3,1)/sph(1,1)) ;
if (cart(2,1)>=0)
   sph(3,1)=acos(cart(1,1)/sqrt(cart(1,1)*cart(1,1)+cart(2,1)*cart(2,1))) ;
else
   sph(3,1)=2*M_PI-acos(cart(1,1)/sqrt(cart(1,1)*cart(1,1)+cart(2,1)*cart(2,1))) ;
return sph ;
}*/
Vector Geometrie::cart2sph (Vector cart)
{
Vector sph ;
sph(1)=sqrt(cart(1)*cart(1)+cart(2)*cart(2)+cart(3)*cart(3)) ;
sph(2)=acos(cart(3)/sph(1)) ;
if (cart(2)>=0)
   sph(3)=acos(cart(1)/sqrt(cart(1)*cart(1)+cart(2)*cart(2))) ;
else
   sph(3)=2*M_PI-acos(cart(1)/sqrt(cart(1)*cart(1)+cart(2)*cart(2))) ;
return sph ;
}

Vector Geometrie::sph2cart (Vector sph)
{
Vector cart ;

cart(3)=sph(1)*cos(sph(2)) ;
cart(1)=sph(1)*sin(sph(2))*cos(sph(3)) ;
cart(2)=sph(1)*sin(sph(2))*sin(sph(3)) ;

return cart ;
}


Vector Geometrie::gravicentre (Vector pt1, Vector pt2, Vector pt3)
{
Vector *res= new Vector ;
(*res)=((pt1+pt2+pt3)/3.) ;
return (*res) ;
}
//------------------------------------------------------------
// Rotation avec Vecteur/Matrice et MatriceRotation => Resutat Vecteur/Matrice
/*Matrix Geometrie::rotation (Matrix vecteur, Matrix rotation)
{ return rotation*vecteur ;}*/
Matrix3x3 Geometrie::rotation (Matrix3x3 mat, Matrix3x3 rotation)
{ return rotation*mat ;}
Vector Geometrie::rotation (Vector vecteur, Matrix3x3 rotation)
{ return rotation*vecteur ; }
//--------------------------------------------------------------------
// Rotation avec Vecteur/Matrice et angle + axe
/*Matrix Geometrie::rotation (Matrix vecteur, double angle, Matrix axe)
{
Matrix rot ;
rot=Geometrie::get_rot_matrix(angle, axe) ;
return Geometrie::rotation(vecteur, rot) ;
}*/
Vector Geometrie::rotation (Vector vecteur, double angle, Vector axe)
{
Matrix3x3 rot ;
rot=Geometrie::get_rot_matrix(angle, axe) ;
return Geometrie::rotation(vecteur, rot) ;
}
Matrix3x3 Geometrie::rotation (Matrix3x3 vecteur, double angle, Vector axe)
{
Matrix3x3 rot ;
rot=Geometrie::get_rot_matrix(angle, axe) ;
return Geometrie::rotation(vecteur, rot) ;
}
//-------------------------------------------------------------
// Rotations d'axe Z
/*Matrix Geometrie::rotation_axeZ (Matrix vecteur, double angle)
{
Matrix res ;
Matrix axe(3,1) ;
axe << 0.0 << 0.0 << 1.0 ;
res=Geometrie::rotation(vecteur, angle, axe) ;
return res;
}*/
Matrix3x3 Geometrie::rotation_axeZ (Matrix3x3 mat, double angle)
{
Matrix3x3 res ;
Vector axe(3) ;
axe(1)=axe(2)=0 ; axe(3)=1 ;
res=Geometrie::rotation(mat, angle, axe) ;
return res;
}
Vector Geometrie::rotation_axeZ (Vector vec, double angle)
{
Vector res ;
Vector axe(3) ;
axe(1)=axe(2)=0 ; axe(3)=1 ;
res=Geometrie::rotation(vec, angle, axe) ;
return res;
}
//---------------------------------------------------------------
double Geometrie::distance_cylindre (Cylindre & C, Vector x)
{
 Vector x2, n_axe, contact ;
 double dst_axe, norm_axe ; int sgn ;

 x2=x-C.centre ;
 dst_axe=x2.dot(C.axe) ;
 if (dst_axe<0) {dst_axe=-dst_axe ; C.axe=-C.axe ; }
 n_axe=x2-(dst_axe*C.axe) ;
 norm_axe=n_axe.norm() ;


 if (norm_axe!=0)
    n_axe=n_axe/norm_axe ;
 else
    n_axe=0 ;

 if (dst_axe<C.L/2 && norm_axe<C.D/2) // inside
 {
   sgn=-1 ;
   if (dst_axe>(C.L/2-C.D/2) && (C.D/2-norm_axe)>(C.L/2-dst_axe))
     contact=(C.L/2)*C.axe + norm_axe*n_axe ;
   else
     contact=(C.D/2)*n_axe+dst_axe*C.axe ;
 }
 else
 {
   sgn=1 ;
   contact=n_axe*min(C.D/2, norm_axe)+C.axe*min(C.L/2, dst_axe) ;
 }
 x2=x2-contact ;
 return (x2.norm()*sgn) ;
}

//---------------------------------------------------------------
Vector Geometrie::contact_cylindre (Cylindre & C, Vector x)
{
 Vector x2, n_axe, contact ;
 double dst_axe, norm_axe ; int sgn ;

 x2=x-C.centre ;
 dst_axe=x2.dot(C.axe) ;
 if (dst_axe<0) {dst_axe=-dst_axe ; C.axe=-C.axe ; }
 n_axe=x2-(dst_axe*C.axe) ;
 norm_axe=n_axe.norm() ;


 if (norm_axe!=0)
    n_axe=n_axe/norm_axe ;
 else
    n_axe=0 ;

 if (dst_axe<C.L/2 && norm_axe<C.D/2) // inside
 {
   sgn=-1 ;
   if (dst_axe>(C.L/2-C.D/2) && (C.D/2-norm_axe)>(C.L/2-dst_axe))
     contact=(C.L/2)*C.axe + norm_axe*n_axe ;
   else
     contact=(C.D/2)*n_axe+dst_axe*C.axe ;
 }
 else
 {
   sgn=1 ;
   contact=n_axe*min(C.D/2, norm_axe)+C.axe*min(C.L/2, dst_axe) ;
 }

 contact=C.centre+contact ;
 return (contact) ;
}

//-------------------------------------------------------
/*Matrix Geometrie::get_rot_matrix (double angle, Matrix axe)
{
Matrix rotation(3,3), utu(3,3), antisymu(3,3) ;
IdentityMatrix Id(3) ;

utu=axe*(axe.t()) ;
antisymu << 0.0       << -axe(3,1)  << axe(2,1)
         << axe(3,1)  << 0.0        << -axe(1,1)
         << -axe(2,1) << axe(1,1)   << 0.0 ;
rotation=utu+(Id-utu)*cos(angle)+antisymu*sin(angle);
return rotation ;
}*/
//-------------------------------------------------------
Matrix3x3 Geometrie::get_rot_matrix (double angle, Vector axe)
{
Matrix3x3 rotation, utu, antisymu, cosine;
Matrix3x3 Id("Id") ;
Vector tmp ;
tmp=axe.t() ;
utu=axe.multiply(tmp) ;
antisymu(1,1)=antisymu(2,2)=antisymu(3,3)=0.0 ;
antisymu(1,2)=-axe(3) ;
antisymu(1,3)=axe(2) ;
antisymu(2,1)=axe(3) ;
antisymu(3,1)=-axe(2) ;
antisymu(2,3)=-axe(1) ;
antisymu(3,2)=axe(1) ;
rotation=utu+(Id-utu)*cos(angle)+antisymu*sin(angle) ;
return rotation ;
}
//---------------------------------------------------------
double Geometrie::tri_surface (Vector a, Vector b, Vector c)
{
 Vector u, v ;
 u=b-a ;
 v=c-a ;
 u=u.cross(v) ;
 return ((1/2.0)*u.norm()) ;
}

double Geometrie::solid_angle (Vector a, Vector b, Vector c)
{
 double d, n, r ;
 n=abs(triple_product(a,b,c)) ;
 d=a.norm()*b.norm()*c.norm()+ c.norm()*(a.dot(b))+ b.norm()*(a.dot(c)) + a.norm()*(b.dot(c)) ;
 r=atan(n/d) ; if (r<0) r+=M_PI ;
 return (2*r) ;
}

//==================================================
/*Matrix Calcul::cross_product (Matrix a, Matrix b)
{
Matrix res(3,1) ;
res(1,1)=a(2,1)*b(3,1) - a(3,1)*b(2,1) ;
res(2,1)=a(3,1)*b(1,1) - a(1,1)*b(3,1) ;
res(3,1)=a(1,1)*b(2,1) - a(2,1)*b(1,1) ;
return res ;
}*/
Vector Calcul::cross_product (Vector a, Vector b)
{
Vector res(3) ;
res(1)=a(2)*b(3) - a(3)*b(2) ;
res(2)=a(3)*b(1) - a(1)*b(3) ;
res(3)=a(1)*b(2) - a(2)*b(1) ;
return res ;
}

double Calcul::angle_0_2pi(double angle)
{
while (angle>2*M_PI) {angle-=2*M_PI ; }
while (angle<0) {angle+=2*M_PI ; }
return angle ;
}

double Calcul::true_random (void)
{
 static FILE * rf = NULL ;  int r ; size_t res ;
 if (rf==NULL) rf=fopen("/dev/urandom", "rb") ;
 res=fread(&r, 4,1, rf) ; r=r& (~(1<<31)) ;
 return ((double)r/(double)(~(1<<31))) ;
}

double Calcul::blockstd (double * x, int N)
{
// Based on H. Flyvbjerg and H. G. Petersen, "Error estimates on averages of
// correlated data", J. Chem. Phys. 91, 461 (1989); doi: 10.1063/1.457480

int i=0,j ;
double *c0, *sm, max ;

c0=(double*)malloc(sizeof(double)*ceil(log2(N))) ;
sm=(double*)malloc(sizeof(double)*ceil(log2(N))) ;

while (N>4)
{
 c0[i]=var(x, N) ;
 sm[i]=c0[i]/(double)(N-1) ;

 for (j=0 ; j<N/2 ; j++) x[j]=(x[2*j]+x[2*j+1])/2 ;
 N=N/2 ;
 i=i+1 ;
}

max=0 ;
for (j=0 ; j<round(0.8*log2(N)) ; j++)
  if (fabs(sm[j])>max) max=sm[j] ;

free(c0) ; free(sm) ;
return (sqrt(max)) ;
}
//-----------------------------
double Calcul::mean (double *x, int N)
{
 double r=0 ;
 for (int i=0 ; i<N ; i++) r+=x[i] ;
 return (r/(double(N))) ;
}
//------------------------------
double Calcul::var (double *x, int N)
{
 double std=0 ;
 double moy ;
 moy=mean(x, N) ;
 for (int i=0 ; i<N ; i++) std+=((x[i]-moy)*(x[i]-moy)) ;
 return (std/(double)(N)) ;
}

int Calcul::eigen (Matrix3d A, Vector3d & eigval, Matrix3d & eigvec)
{
  Matrix3d modif ;
  EigenSolver<Matrix3d> es;
  es.compute(A, true);
  eigval=es.eigenvalues().real() ;
  eigvec=es.eigenvectors().real() ; // Columns are the eigenvectors

  if (eigval(0)<=eigval(1) && eigval(1)<=eigval(2))
  {modif<< 0,0,1, 0,1,0, 1,0,0 ;  }
  else if (eigval(0)<=eigval(2) && eigval(2)<=eigval(1))
  { modif<< 0,0,1, 1,0,0, 0,1,0 ; }
  else if (eigval(1)<=eigval(0) && eigval(0)<=eigval(2))
  {  modif<< 0,1,0, 0,0,1, 1,0,0 ; }
  else if (eigval(1)<=eigval(2) && eigval(2)<=eigval(0))
   { modif<< 1,0,0, 0,0,1, 0,1,0 ;  }
  else if (eigval(2)<=eigval(0) && eigval(0)<=eigval(1))
   { modif<< 0,1,0, 1,0,0, 0,0,1 ; }
  else if (eigval(2)<=eigval(1) && eigval(1)<=eigval(0))
   { modif<< 1,0,0, 0,1,0, 0,0,1 ; }
  else modif<< 1,0,0, 0,1,0, 0,0,1 ;

  eigvec=eigvec*modif ;
  eigval=(eigval.transpose()*modif).transpose() ;
return 0 ;
}


//============================================================


//=================================================
void Fonction::gaussienne1D (double & res, const double & x, const double & sigma)
{
res=1/(sigma*sqrt(2*M_PI))*exp(-(x*x)/(2*sigma*sigma)) ;
}
//..........
void Fonction::creneau1D (double & res, const double & x, const double & w)
{
if (x<=-w/2 || x>w/2) {res=0 ; return ;}
else {res=1/w ; return ;}
}
//..........
void Fonction::gaussienne3D (double & res, const double & x, const double & w, double radius)
{
res=(1/(pow(M_PI,3./2.)*w*w*w)*exp(-(x/w)*(x/w))) ;
}
//..........
void Fonction::gaussienne3D_multivariee (double & res, Vector &v, Vector &w, double radius )
{
double prefactor ; prefactor=1/(pow(M_PI,3./2.)*(w(1))*w(2)*w(3)) ;
res=(prefactor*exp(-(v(1)/w(1))*(v(1)/w(1))-(v(2)/w(2))*(v(2)/w(2))-(v(3)/w(3))*(v(3)/w(3)))) ;
}
//..........
void Fonction::gaussienne3D (double & res, Vector &x, Vector &w, double radius) {double norme=Calcul::norm(x) ; Fonction::gaussienne3D(res, norme, w(1), radius) ; }
//..........
void Fonction::creneau3D (double & res, Vector &x, Vector &w, double radius)
{
double V=w(1)*w(2)*w(3) ;
if (x(1)>-w(1)/2. && x(1)<w(1)/2. && x(2)>-w(2)/2. && x(2)<w(2)/2. && x(3)>-w(3)/2. && x(3)<w(3)/2.)
   {res=1/V ; return ;}
else {res=0 ; return ;}
}
//..........
void Fonction::creneau3D_fin (double &res, Vector x, Vector w, double radius)
{
double V=w(1)*w(2)*w(3) ;
double d=radius*2 ;
x(1)=fabs(x(1)) ; x(2)=fabs(x(2)) ; x(3)=fabs(x(3)) ;
w=w/2. ;

if (x(1)<(w(1)-d) && x(2)<(w(2)-d) && x(3)<(w(3)-d)) // (1 cas) : complètement à l'intérieur
   {res=1/V ; return ;}
else if (x(1)>(w(1)+d) || x(2)>(w(2)+d) || x(3)>(w(3)+d)) //(19 cas) : complètement à l'extérieur
   {res=0 ; return ;}
else if ( (x(1)<(w(1)+d) && x(1)>w(1)-d) ^ (x(2)<(w(2)+d) && x(2)>w(2)-d) ^ (x(3)<(w(3)+d) && x(3)>w(3)-d) ) { return ;}  //(3 cas) : un seul plan de coupe
else if ( x(1)<(w(1)-d) ^ x(2)<(w(2)-d) ^ x(3)<w(3)-d ) {return ;} // (3 cas) : 2 plans de coupe
else {return ;}		// (1cas) : 3 plans de coupe

}
//..........
void Fonction::intersect_sphere (double & res, Vector &x, Vector &w, double radius)
{
double V=4/3.0*M_PI*w(1)*w(1)*w(1) ;
double xnorm=Calcul::norm(x) ;
double h1, h2, V1, V2, r1, r2, Vp ;

if (actions["delcyl"].set) V=res ;

if (xnorm<w(1)-radius) {res=1/V ; return ; }
else if (xnorm>w(1)+radius) {res=0 ; return ; }
else
	{
	r1 = w(1) ; r2 = radius ;
	Vp=4/3.0*M_PI*r2*r2*r2 ;
	h1 = (r1*r1-r2*r2+xnorm*xnorm)/(2*xnorm) ;
	h2 = (r2*r2-r1*r1+xnorm*xnorm)/(2*xnorm) ;
	V1 = M_PI*(r1*r1*r1*2/3.0 - r1*r1*h1 + h1*h1*h1/3.0) ;
	V2 = M_PI*(r2*r2*r2*2/3.0 - r2*r2*h2 + h2*h2*h2/3.0) ;
	res=(V1+V2) / (Vp*V) ; return ;
	}
}
//..........
double Fonction::aireleft_intersect_cercle (double R, double r, double x)
{
// 1er cercle de rayon R à l'origine, second de rayon r centré en (x,0,0), x>0.
 double Ainter ;
 double h1, h2, theta1, theta2, V1, V2 ;
 if (x>r+R) Ainter=0 ;
 else
 {
  if (R<r && x<r-R ) Ainter=M_PI*R*R ;
  else if (R>r && x<R-r) Ainter=M_PI*r*r ;
  else
  {
    h1 = (R*R-r*r+x*x)/(2*x) ;
    h2 = (r*r-R*R+x*x)/(2*x) ;
    theta1 = acos (h1/R) ; theta2 = acos (h2/r) ;
    V1 = R*R*theta1-sin(theta1)*R*h1 ;
    V2 = r*r*theta2-sin(theta2)*r*h2 ;
    Ainter=V1+V2 ;
  }
 }
 return (M_PI*R*R-Ainter) ;
}
double Fonction::volume_intersect_sphere_cylinder (Vector xsphere, double Rsphere, Vector xcyl, Vector dircyl, double Rcyl, double Lcyl)
{
  int i, maxstep=100000 ;
  double b, zcyl, dz, zplane ;
  double theta, Rcercle, volume ;
  // Axis change
  xcyl=xcyl-xsphere ; // Changement d'origine
  theta=acos(xcyl.dot(dircyl)/xcyl.norm()) ;
  b=sin(theta)*xcyl.norm() ; // Calcul du parametre d'impact
  zcyl=sqrt(xcyl.norm()*xcyl.norm()-b*b) ;

 // Intégration numérique
  dz=2*Rsphere/maxstep ; zplane=-Rsphere+dz/2 ;
  volume=0 ;
  for (i=0 ; i<maxstep ; i++)
  {
    Rcercle=sqrt(Rsphere*Rsphere-zplane*zplane) ;
    if (zplane>zcyl-Lcyl/2 && zplane<zcyl+Lcyl/2)
      volume+=aireleft_intersect_cercle(Rcercle,Rcyl,b)*dz ;
    else
      volume+=M_PI*Rcercle*Rcercle*dz ;

    zplane+=dz ;
  }
  return volume ;
}
//.......................................
double Fonction::volume_intersect_pave_cylindre(Vector pos, Vector sigmafenetre, double R)
{
  //if (actions["is2D"].set) DISP_Err("2D not IMPLEMENTED for volume_intersect_pave_cylindre") ;
  Vector p1, p2, p3, p4 ;
  double Volume=sigmafenetre(1)*sigmafenetre(2)*sigmafenetre(3) ;

  p1(1)=pos(1)-sigmafenetre(1)/2. ; p1(2)=pos(2)-sigmafenetre(2)/2. ; p1(3)=0 ;
  p2(1)=pos(1)-sigmafenetre(1)/2. ; p2(2)=pos(2)+sigmafenetre(2)/2. ; p2(3)=0 ;
  p3(1)=pos(1)+sigmafenetre(1)/2. ; p3(2)=pos(2)-sigmafenetre(2)/2. ; p3(3)=0 ;
  p4(1)=pos(1)+sigmafenetre(1)/2. ; p4(2)=pos(2)+sigmafenetre(2)/2. ; p4(3)=0 ;

  //  printf("%g %g %g %g %g\n", p1.norm(), p2.norm(), p3.norm(), p4.norm(), R ) , fflush(stdout) ;
  if (p1.norm()<=R && p2.norm()<=R && p3.norm()<=R && p4.norm()<=R) return Volume ;
  if (p1.norm()>=R && p2.norm()>=R && p3.norm()>=R && p4.norm()>=R) return 0 ;

  // Sinon, MonteCarlo methode
#ifndef BOOST
  DISP_Err("Impossible d'effectuer le calcul si le programme n'a pas été compilé avec BOOST") ;
  return Volume ;
#else
  static boost::random::mt19937 gen(std::time(0));
  boost::random::uniform_real_distribution<> distx(p1(1),p3(1)), disty(p1(2),p2(2))  ;
  int cpt=0, i ; Vector u ;
  u(3)=0 ;
  for (i=0 ; i<MCMAXITER ; i++)
  {
    u(1)=distx(gen) ; u(2)=disty(gen) ;
    if (u.norm()<R) cpt++ ;
  }
  return (cpt/((double)MCMAXITER)*Volume) ;
#endif
}

//...........................................................................................
void Fonction::intersect_cercle (double & res, Vector &x, Vector &w, double radius)
{
double V=M_PI*w(1)*w(1) ;
double xnorm=Calcul::norm2D(x) ;
double h1, h2, V1, V2, r1, r2, Vp, theta1, theta2 ;
static bool warn_once=true ;

//if (radius>w(1)) {DISP_Warn("Attention : rayon trop grand pour le coarse !") ; printf("%g %g ", radius, w(1)) ; }
if (xnorm<fabs(w(1)-radius))
{
  if (radius<=w(1)) {res=1/V ; return ; }
  else
  {
    if (warn_once) DISP_Warn("Attention : rayon trop grand pour le coarse (on devrait s'en sortir, mais la mesure de phi risque d'être bof ...)!\n") ;
    warn_once=false ;
    res=1/(M_PI*radius*radius) ; return ;
  }
}
else if (xnorm>w(1)+radius) {res=0 ; return ; }
else
	{
	r1 = w(1) ; r2 = radius ;
	Vp=M_PI*r2*r2 ;
	h1 = (r1*r1-r2*r2+xnorm*xnorm)/(2*xnorm) ;
	h2 = (r2*r2-r1*r1+xnorm*xnorm)/(2*xnorm) ;
	theta1 = acos (h1/r1) ; theta2 = acos (h2/r2) ;
	V1 = r1*r1*theta1-sin(theta1)*r1*h1 ;
	V2 = r2*r2*theta2-sin(theta2)*r2*h2 ;
	res=(V1+V2) / (Vp*V) ; return ;
	}
}
//..........
void Fonction::intersect_creneau2D (double &res, Vector &x, Vector &w, double radius)
{
double V=w(1)*w(2) ;
double Vp=M_PI*radius*radius ;
double A, X, Y ;

if (2*radius > w(1) || 2*radius > w(2)) DISP_Warn ("Attention : rayon trop grand pour le coarse rect2D !") ;

// On se place dans le cadran x+y+
X=fabs(x(1)) ; Y=fabs(x(2)) ;

if (X>w(1)/2+radius || Y>w(2)/2+radius) {res=0 ; return ;}
else if (X<w(1)/2-radius && Y<w(2)/2-radius) {res=1/V ; return ; }
else
{
 double l, l1, l2, angle ;
 A=2*Vp ;
 if (Y<w(2)/2-radius) // l'arrête verticale coupe seule
 {
   l=w(1)/2-X ;
   A=acos(Calcul::fborne(-l/radius,-1,1))*radius*radius+l*sqrt(fmax(0,radius*radius-l*l)) ;
 }
 else if (X<w(1)/2-radius) // l'arrête horizontale coupe seule
 {
   l=w(2)/2-Y ;
   A=acos(Calcul::fborne(-l/radius,-1,1))*radius*radius+l*sqrt(fmax(0,radius*radius-l*l)) ;
 }
 else // On est près du sommet
 {
   l1=w(1)/2-X ; l2=w(2)/2-Y ;
   if (l1*l1+l2*l2 > radius*radius) // Le sommet est en dehors du cercle
   {
     if ( l1<0 && l2<0) { res=0 ; return ; }
     else if (l1<0 && l2>0) // En fait, il n'y a d'un segment qui coupe
       A=acos(fmin(1,-l1/radius))*radius*radius+l1*sqrt(fmax(0,radius*radius-l1*l1)) ;
     else if (l2<0 && l1>0)
       A=acos(fmin(1,-l2/radius))*radius*radius+l2*sqrt(fmax(0,radius*radius-l2*l2)) ;
     else // 2 segments de coupe mais sommet en dehors du cercle quand même
       A=radius*radius*(M_PI - acos(fmin(1,l1/radius)) - acos(fmin(1,l2/radius))) + l1*sqrt(fmax(0,radius*radius-l1*l1)) + l2*sqrt(fmax(0,radius*radius-l2*l2)) ;
  }
   else
   {
     angle=M_PI/2.-asin(-l2/radius)-asin(-l1/radius) ;
     A=angle/2.*radius*radius - cos(angle/2.)*sin(angle/2.)*radius*radius + 0.5*(sqrt(fmax(0,radius*radius-l1*l1))+l2)*(sqrt(fmax(0,radius*radius-l2*l2))+l1) ;
   }
 }

res=A/(Vp*V) ; return ;
}
}
//..........
void Fonction::intersect_tranche(double &res, Vector &x, Vector &w, double radius, char dir)
{
  double V=w(1)*w(2)*w(3), Vin ; double h ;
  int orient ;
  double Vp=M_PI*radius*radius*radius*4/3. ;
  switch (dir)
  {
    case 'x': orient=1 ; break ;
    case 'y': orient=2 ; break ;
    case 'z': orient=3 ; break ;
    default : orient=1 ; DISP_Err("Orientation tranche incorrecte") ; break ;
  }
  if (2*radius > w(orient)) DISP_Warn("Attention : rayon trop grand pour le coarse tranche !") ;
  if (x(orient)>w(orient)/2.+radius || x(orient)<-w(orient)/2.-radius) {res=0 ; return ; }
  else if (x(orient)<w(orient)/2.-radius && x(orient)>-w(orient)/2.+radius) {res=1/V ; return ; }
  else
  {
    h=radius-(fabs(x(orient))-w(orient)/2.) ;
    Vin=M_PI*h*h/3*(3*radius-h) ;
    res=Vin/Vp/V ;
    return ;
  }
}
//..........
void Fonction::int_intersect_sphere (double & res , Vector &xa, Vector &xb, Vector &w)
{
cout << "ERR : int_intersect_sphere : NOT IMPLEMENTED YET" ;
}

void Fonction::int_gaussienne3D(double & res, Vector &xa, Vector &xb, Vector &w)
{
double normea=Calcul::norm(xa), normeb=Calcul::norm(xb), ra, rb ;
//Fonction::gaussienne3D(ra,normea, w(1), double radius) ; Fonction::gaussienne3D(rb, normeb, w(1), double radius) ;
//res=0.5*(ra+rb) ;
DISP_Err("NOT IMPLEMENTED");
}
//..........
void Fonction::int_creneau3D (double & res, Vector &xa, Vector &xb, Vector &w)
{
Vector alpha ;
Vector a,b,d, moins1, plus1, w2, delta ;
double V ;

// Calcul du volume
moins1=-1 ;
plus1=1 ;
V=w(1)*w(2)*w(3) ;

if (actions["deltank"].set) V=res ;
// Les deux sont dedans ? les deux sont dehors ?
w2=w/2. ;
a=xa.div(w2) ;
b=xb.div(w2) ;
if (a>moins1 && a<plus1 && b>moins1 && b<plus1) {res=1/V ; return ;}
// C'est à cheval donc ...
else if (a>moins1 && a<plus1) // C'est xa qui est dans le cube, xb est en dehors
 {
  // Recherche du alpha qui coupe sur les bornes à +w/2
  delta=xb-xa ;
  w2=w/2.-xa ;
  alpha=w2.div(delta) ;
  if (alpha(1)>=0 && alpha(1)<=1) {res=alpha(1)/V ; return ; }
  if (alpha(2)>=0 && alpha(2)<=1) {res=alpha(2)/V ; return ; }
  if (alpha(3)>=0 && alpha(3)<=1) {res=alpha(3)/V ; return ; }
  w2=-w/2.-xa ;
  alpha=w2.div(delta) ;
  if (alpha(1)>=0 && alpha(1)<=1) {res=alpha(1)/V ; return ; }
  if (alpha(2)>=0 && alpha(2)<=1) {res=alpha(2)/V ; return ; }
  if (alpha(3)>=0 && alpha(3)<=1) {res=alpha(3)/V ; return ; }
  cout << "INTCRENEAU3DERRA" << fflush(stdout) ;
 }
else if (b>moins1 && b<plus1) // C'est xa qui est en dehors du cube, xb est dedans
 {
  delta=xa-xb ;
  // Recherche du alpha qui coupe sur les bornes à +w/2
  w2=w/2.-xb ;
  alpha=w2.div(delta) ;
  if (alpha(1)>=0 && alpha(1)<=1) {res=alpha(1)/V ; return ; }
  if (alpha(2)>=0 && alpha(2)<=1) {res=alpha(2)/V ; return ; }
  if (alpha(3)>=0 && alpha(3)<=1) {res=alpha(3)/V ; return ; }
  w2=-w/2.-xb ;
  alpha=w2.div(delta) ;
  if (alpha(1)>=0 && alpha(1)<=1) {res=alpha(1)/V ; return ; }
  if (alpha(2)>=0 && alpha(2)<=1) {res=alpha(2)/V ; return ; }
  if (alpha(3)>=0 && alpha(3)<=1) {res=alpha(3)/V ; return ; }
  cout << "INTCRENEAU3DERRB" << fflush(stdout) ;
 }
else
 {
   res=0 ; return ;
 }
/*
alpha=d(1,1)/b(1,1) ;
if (alpha>0 && alpha<1)
   {res=1/V*(b(1,1)<0?(1-alpha):alpha) ; return ;}
alpha=d(2,1)/b(2,1) ;
if (alpha>0 && alpha<1)
   {res=1/V*(b(2,1)<0?(1-alpha):alpha) ;return ;}
alpha=d(3,1)/b(3,1) ;
if (alpha>0 && alpha<1)
   {res=1/V*(b(3,1)<0?(1-alpha):alpha) ; return ;}

d=-w/2.-xa ;
alpha=d(1,1)/b(1,1) ;
if (alpha>0 && alpha<1)
   {res=1/V*(b(1,1)>0?(1-alpha):alpha) ; return ;}
alpha=d(2,1)/b(2,1) ;
if (alpha>0 && alpha<1)
   {res=1/V*(b(2,1)>0?(1-alpha):alpha) ;return ;}
alpha=d(3,1)/b(3,1) ;
if (alpha>0 && alpha<1)
   {res=1/V*(b(3,1)>0?(1-alpha):alpha) ;return ;}
b=(xb+xa)/2. ;
if (b(1,1)<w(1,1)/2. && b(1,1)>-w(1,1)/2.)
   {res=1/V ; return ;}
else
   {res=0 ;return ;}*/
}

void Fonction::int_creneau2D (double & res, Vector &xa, Vector &xb, Vector &w)
{
Vector alpha ;
Vector a,b,d, moins1, plus1, w2, delta ;
double V ;

// Calcul du volume
moins1=-1 ;
plus1=1 ;
V=w(1)*w(2);
// Les deux sont dedans ? les deux sont dehors ?
w2=w/2. ;
a=xa.div(w2) ;
b=xb.div(w2) ;
if (a>moins1 && a<plus1 && b>moins1 && b<plus1) {res=1/V ; return ;}
// C'est à cheval donc ...
else if (a>moins1 && a<plus1) // C'est xa qui est dans le cube, xb est en dehors
 {
  // Recherche du alpha qui coupe sur les bornes à +w/2
  delta=xb-xa ;
  w2=w/2.-xa ;
  alpha=w2.div(delta) ;
  if (alpha(1)>=0 && alpha(1)<=1) {res=alpha(1)/V ; return ; }
  if (alpha(2)>=0 && alpha(2)<=1) {res=alpha(2)/V ; return ; }
  w2=-w/2.-xa ;
  alpha=w2.div(delta) ;
  if (alpha(1)>=0 && alpha(1)<=1) {res=alpha(1)/V ; return ; }
  if (alpha(2)>=0 && alpha(2)<=1) {res=alpha(2)/V ; return ; }
  cout << "INTCRENEAU2DERRA" << fflush(stdout) ;
 }
else if (b>moins1 && b<plus1) // C'est xa qui est en dehors du cube, xb est dedans
 {
  delta=xa-xb ;
  // Recherche du alpha qui coupe sur les bornes à +w/2
  w2=w/2.-xb ;
  alpha=w2.div(delta) ;
  if (alpha(1)>=0 && alpha(1)<=1) {res=alpha(1)/V ; return ; }
  if (alpha(2)>=0 && alpha(2)<=1) {res=alpha(2)/V ; return ; }
  w2=-w/2.-xb ;
  alpha=w2.div(delta) ;
  if (alpha(1)>=0 && alpha(1)<=1) {res=alpha(1)/V ; return ; }
  if (alpha(2)>=0 && alpha(2)<=1) {res=alpha(2)/V ; return ; }
  cout << "INTCRENEAU2DERRB" << fflush(stdout) ;
 }
else
 {
   res=0 ; return ;
 }
}

int Fonction::intersect3D_RayTriangle(Vector ray, Vector t1, Vector t2, Vector t3)
{
//    Input:  a ray R, and a triangle T
//    Output: *I = intersection point (when it exists)
//    Return: -1 = triangle is degenerate (a segment or point)
//             0 =  disjoint (no intersect)
//             1 =  intersect in unique point I1
//             2 =  are in the same plane

    Vector    u, v, n;              // triangle vectors
    Vector    w0, w;           // ray vectors
    float     r, a, b;              // params to calc ray-plane intersect

    Vector I ;

    // get triangle edge vectors and plane normal
    u = t2-t1;
    v = t3-t1;
    n = u.cross(v);              // cross product

    w0 = - t1; // w0=origin-t1 pour être complet, avec des vecteur ray dont l'origine n'est pas en 0
    a = -n.dot(w0);
    b = n.dot(ray);
    if (fabs(b) < 0.0000001) {     // ray is  parallel to triangle plane
        if (a == 0)                 // ray lies in triangle plane
            return 2;
        else return 0;              // ray disjoint from plane
    }

    // get intersect point of ray with triangle plane
    r = a / b;
    if (r < 0.0)                    // ray goes away from triangle
        return 0;                   // => no intersect
    // for a segment, also test if (r > 1.0) => no intersect

    //*I = R.P0 + r * dir;            // intersect point of ray and plane
    I=r*ray ;

    // is I inside T?
    float    uu, uv, vv, wu, wv, D;
    uu = u.dot(u);
    uv = u.dot(v);
    vv = v.dot(v);
    w = I - t1;
    wu = w.dot(u);
    wv = w.dot(v);
    D = uv * uv - uu * vv;

    // get and test parametric coords
    float s, t;
    s = (uv * wv - vv * wu) / D;
    if (s < 0.0 || s > 1.0)         // I is outside T
        return 0;
    t = (uv * wu - uu * wv) / D;
    if (t < 0.0 || (s + t) > 1.0)  // I is outside T
        return 0;

    return 1;                       // I is in T
}



//========================================================================================
//=======================Matrix et Vector=================================================
//-------------------------
Matrix3x3 Id("Id") ;
//--------------
Matrix3x3::Matrix3x3()
{for (i=0 ; i<9 ; i++) {dat[i]=0;} tmp=false; }
Matrix3x3::Matrix3x3(const char * str)
{for (i=0 ; i<9 ; i++) {dat[i]=0;} if (!strcmp(str, "Id")) {dat[0]=dat[4]=dat[8]=1 ;} tmp=false ;}
//--------------
void Matrix3x3::disp(void)
{for (i=0 ; i<9 ; i++) {cout << dat[i] ; cout << (i%3==2?"\n":" ") ; } cout << "\n" ; }
//--------------
/*void Matrix3x3::operator= (Matrix3x3 & m)
{for (i=0 ; i<9 ; i++) {dat[i]=m[i];} }*/
//-------------
Matrix3x3 Matrix3x3::t (void)
{
Matrix3x3 res;
res[0]=dat[0] ; res[1]=dat[3] ; res[2]=dat[6] ;
res[3]=dat[1] ; res[4]=dat[4] ; res[5]=dat[7] ;
res[6]=dat[2] ; res[7]=dat[5] ; res[8]=dat[8] ;
res.tmp=true ;
return res ;
}
//--------------
Matrix3x3 Matrix3x3::operator+ (Matrix3x3 m)
{
 Matrix3x3 res ;
 for (i=0 ; i<9 ; i++)
 {res[i]=dat[i]+m[i] ; }
 res.tmp=true ;
 return res ;
}
//--------------
Matrix3x3 Matrix3x3::operator+ (double x)
{
 Matrix3x3 res ;
 for (i=0 ; i<9 ; i++)
 {res[i]=dat[i]+x ; }
 res.tmp=true ;
 return res ;
}
//---------------
Matrix3x3 Matrix3x3::operator- ()
{
 Matrix3x3 res ;
 for (i=0 ; i<9 ; i++)
 {res[i]=-dat[i] ; }
 res.tmp=true ;
 return res ;
}
//--------------
Matrix3x3 Matrix3x3::operator* (double x)
{
 Matrix3x3 res ;
 for (i=0 ; i<9 ; i++)
 {res[i]=dat[i]*x ; }
 res.tmp=true ;
 return res ;
}
//--------------
Matrix3x3 Matrix3x3::operator* (Matrix3x3 & m)
{
 Matrix3x3 res ;
 res[0]=dat[0]*m[0]+dat[1]*m[3]+dat[2]*m[6] ;
 res[1]=dat[0]*m[1]+dat[1]*m[4]+dat[2]*m[7] ;
 res[2]=dat[0]*m[2]+dat[1]*m[5]+dat[2]*m[8] ;
 res[3]=dat[3]*m[0]+dat[4]*m[3]+dat[5]*m[6] ;
 res[4]=dat[3]*m[1]+dat[4]*m[4]+dat[5]*m[7] ;
 res[5]=dat[3]*m[2]+dat[4]*m[5]+dat[5]*m[8] ;
 res[6]=dat[6]*m[0]+dat[7]*m[3]+dat[8]*m[6] ;
 res[7]=dat[6]*m[1]+dat[7]*m[4]+dat[8]*m[7] ;
 res[8]=dat[6]*m[2]+dat[7]*m[5]+dat[8]*m[8] ;
 res.tmp=true ;
 return res ;
}
//--------------
Matrix3x3 Matrix3x3::operator/ (double x)
{
 Matrix3x3 res ;
 for (i=0 ; i<9 ; i++)
 {res[i]=dat[i]/x ; }
 res.tmp=true ;
 return res ;
}
//--------------
/*Vector Matrix3x3::operator* (Vector & v)
{
 if (v.iscol()==false) {cerr << "MATRIXERR : multiplication impossible matrix*line" ;}
 else
 {
  Vector res ;
  res[0]=dat[0]*v[0]+dat[1]*v[1]+dat[2]*v[2] ;
  res[1]=dat[3]*v[0]+dat[4]*v[1]+dat[5]*v[2] ;
  res[2]=dat[6]*v[0]+dat[7]*v[1]+dat[8]*v[2] ;
  return res ;
 }
 return v ;
}*/
//--------------
void Vector::disp(void)
{if (iscol()) printf("%f\n%f\n%f\n\n",dat[0],dat[1],dat[2]) ;
 else printf("%f %f %f\n\n",dat[0],dat[1], dat[2]); }
//--------------
Vector operator* (Vector & v, Matrix3x3 & m)
{
 if (v.iscol()) {cerr << "MATRIXERR : multiplication impossible matrix*ligne" ;}
 else
 {
  Vector res(v);
  res[0]=v[0]*m[0]+v[1]*m[3]+v[2]*m[6] ;
  res[1]=v[0]*m[1]+v[1]*m[4]+v[2]*m[7] ;
  res[2]=v[0]*m[2]+v[1]*m[5]+v[2]*m[8] ;
  return res ;
 }
 return 0 ;
}
//--------------
Vector operator* (Matrix3x3 & m, Vector & v)
{
 if (!v.iscol()) {cerr << "MATRIXERR : multiplication impossible matrix*ligne" ;}
 else
 {
  Vector res(v);
  res[0]=v[0]*m[0]+v[1]*m[1]+v[2]*m[2] ;
  res[1]=v[0]*m[3]+v[1]*m[4]+v[2]*m[5] ;
  res[2]=v[0]*m[6]+v[1]*m[7]+v[2]*m[8] ;
  return res ;
 }
 return 0 ;
}
//-------------
Matrix3x3 Vector::multiply (Vector & v)
{
 if (iscol()==false && v.iscol()==true) {cerr << "MATRIXERR : utiliser la fonction .dot() pour faire line*column\n" ;}
 else
 {
  Matrix3x3 res ;
  res[0]=dat[0]*v[0] ;res[1]=dat[0]*v[1] ;res[2]=dat[0]*v[2] ;
  res[3]=dat[1]*v[0] ;res[4]=dat[1]*v[1] ;res[5]=dat[1]*v[2] ;
  res[6]=dat[2]*v[0] ;res[7]=dat[2]*v[1] ;res[8]=dat[2]*v[2] ;
  return res ;
 }
 return 0 ;
}
//-----------
/*Vector Vector::operator% (Vector & v)
{
if (v.iscol()!=iscol()) {cerr <<"MATRIXERR:l'operateur % demande des vecteurs de même orientation (equivalent matlab .*)" ; }
else
 {
 Vector res;
 if (v.iscol()==false) res.autochange() ;
 res[0]=dat[0]*v[0] ;res[1]=dat[1]*v[1] ;res[2]=dat[2]*v[2] ;
 res.tmp=true ;
 return res ;
 }
}*/
//-----------
/*Vector Vector::operator/ (Vector & v)
{
if (v.iscol()!=iscol()) {cerr <<"MATRIXERR:l'operateur / demande des vecteurs de même orientation (equivalent matlab ./)" ; Vector v_null ; return v_null ; }
else
 {
 Vector res;
 if (v.iscol()==false) res.autochange() ;
 res[0]=dat[0]/v[0] ;res[1]=dat[1]/v[1] ;res[2]=dat[2]/v[2] ;
 res.tmp=true ;
 return res ;
 }
}
*/


//==================
Icosahedre::Icosahedre ()
{
  int i ;
  double Phi=(1+sqrt(5))/2. ;
  npts=12 ; nedges=30 ; nfaces=20 ;

  pts=(double **)malloc(sizeof(double*)*npts) ;
  for (i=0 ; i<npts ; i++) pts[i]=(double*)malloc(sizeof(double)*3) ;
  edges=(int **)malloc(sizeof(int*)*nedges) ;
  for (i=0 ; i<nedges ; i++) edges[i]=(int*)malloc(sizeof(int)*2) ;
  faces=(int **)malloc(sizeof(int*)*nfaces) ;
  for (i=0 ; i<nfaces ; i++) faces[i]=(int*)malloc(sizeof(int)*3) ;

  data=(double*) malloc(sizeof(double)*nfaces) ;
  solidangle=(double*) malloc(sizeof(double)*nfaces) ;

  pts[0][0]=0 ; pts[0][1]=1 ;  pts[0][2]=Phi ;
  pts[1][0]=0 ; pts[1][1]=1 ;  pts[1][2]=-Phi ;
  pts[2][0]=0 ; pts[2][1]=-1 ; pts[2][2]=Phi ;
  pts[3][0]=0 ; pts[3][1]=-1 ; pts[3][2]=-Phi ;
  pts[4][0]=Phi ; pts[4][1]=0 ;  pts[4][2]=1 ;
  pts[5][0]=-Phi ; pts[5][1]=0 ;  pts[5][2]=1 ;
  pts[6][0]=Phi ; pts[6][1]=0 ; pts[6][2]=-1 ;
  pts[7][0]=-Phi ; pts[7][1]=0 ; pts[7][2]=-1 ;
  pts[8][0]=1 ; pts[8][1]=Phi ;  pts[8][2]=0 ;
  pts[9][0]=1 ; pts[9][1]=-Phi ;  pts[9][2]=0 ;
  pts[10][0]=-1 ; pts[10][1]=Phi ; pts[10][2]=0 ;
  pts[11][0]=-1 ; pts[11][1]=-Phi ; pts[11][2]=0 ;

  edges[0][0]=0 ; edges[0][1]=2 ;
  edges[1][0]=0 ; edges[1][1]=4 ;
  edges[2][0]=0 ; edges[2][1]=8 ;
  edges[3][0]=0 ; edges[3][1]=10 ;
  edges[4][0]=0 ; edges[4][1]=5 ;
  edges[5][0]=2 ; edges[5][1]=4 ;
  edges[6][0]=4 ; edges[6][1]=8 ;
  edges[7][0]=8 ; edges[7][1]=10 ;
  edges[8][0]=10; edges[8][1]=5 ;
  edges[9][0]=5 ; edges[9][1]=2 ;

  faces[0][0]=0 ; faces[0][1]=2 ; faces[0][2]=4 ;
  faces[1][0]=0 ; faces[1][1]=4 ; faces[1][2]=8 ;
  faces[2][0]=0 ; faces[2][1]=8 ; faces[2][2]=10 ;
  faces[3][0]=0 ; faces[3][1]=10 ; faces[3][2]=5 ;
  faces[4][0]=0 ; faces[4][1]=5 ; faces[4][2]=2 ;

  edges[10][0]=2 ;  edges[10][1]=9 ;
  edges[11][0]=9 ;  edges[11][1]=4 ;
  edges[12][0]=4 ;  edges[12][1]=6 ;
  edges[13][0]=6 ;  edges[13][1]=8 ;
  edges[14][0]=8 ;  edges[14][1]=1 ;
  edges[15][0]=1 ;  edges[15][1]=10 ;
  edges[16][0]=10 ; edges[16][1]=7 ;
  edges[17][0]=7 ;  edges[17][1]=5 ;
  edges[18][0]=5 ;  edges[18][1]=11 ;
  edges[19][0]=11 ; edges[19][1]=2 ;

  faces[5][0]=2 ;  faces[5][1]=9 ;  faces[5][2]=4 ;
  faces[6][0]=9 ;  faces[6][1]=4 ;  faces[6][2]=6 ;
  faces[7][0]=4 ;  faces[7][1]=6 ;  faces[7][2]=8 ;
  faces[8][0]=6 ;  faces[8][1]=8 ;  faces[8][2]=1 ;
  faces[9][0]=8 ;  faces[9][1]=1 ;  faces[9][2]=10 ;
  faces[10][0]=1 ; faces[10][1]=10; faces[10][2]=7 ;
  faces[11][0]=10; faces[11][1]=7 ; faces[11][2]=5 ;
  faces[12][0]=7 ; faces[12][1]=5 ; faces[12][2]=11 ;
  faces[13][0]=5 ; faces[13][1]=11; faces[13][2]=2 ;
  faces[14][0]=11; faces[14][1]=2 ; faces[14][2]=9 ;

  edges[20][0]=9 ;  edges[20][1]=6 ;
  edges[21][0]=6 ;  edges[21][1]=1 ;
  edges[22][0]=1 ;  edges[22][1]=7 ;
  edges[23][0]=7 ;  edges[23][1]=11 ;
  edges[24][0]=11 ; edges[24][1]=9 ;

  edges[25][0]=3 ; edges[25][1]=9 ;
  edges[26][0]=3 ; edges[26][1]=6 ;
  edges[27][0]=3 ; edges[27][1]=1 ;
  edges[28][0]=3 ; edges[28][1]=7 ;
  edges[29][0]=3 ; edges[29][1]=11;

  faces[15][0]=3; faces[15][1]=9; faces[15][2]=6 ;
  faces[16][0]=3; faces[16][1]=6 ; faces[16][2]=1 ;
  faces[17][0]=3; faces[17][1]=1 ; faces[17][2]=7 ;
  faces[18][0]=3; faces[18][1]=7; faces[18][2]=11 ;
  faces[19][0]=3; faces[19][1]=11 ; faces[19][2]=9 ;

  double norme ;
  for (i=0 ; i<npts ; i++)
  {
    norme=sqrt(pts[i][0]*pts[i][0]+pts[i][1]*pts[i][1]+pts[i][2]*pts[i][2]) ;
    pts[i][0]/=norme ; pts[i][1]/=norme ; pts[i][2]/=norme ;
  }

  Vector t1, t2, t3 ;
  for (i=0 ; i<20 ; i++)
  {
    data[i]=0 ;
    t1.set(pts[faces[i][0]][0], pts[faces[i][0]][1], pts[faces[i][0]][2]) ;
    t2.set(pts[faces[i][1]][0], pts[faces[i][1]][1], pts[faces[i][1]][2]) ;
    t3.set(pts[faces[i][2]][0], pts[faces[i][2]][1], pts[faces[i][2]][2]) ;
    solidangle[i]=Geometrie::solid_angle(t1,t2,t3) ;
  }
}



int Icosahedre::write_vtk (string chemin)
{
  FILE *out ; int i, j ;

  out=fopen(chemin.c_str(), "w") ;
  fprintf(out,"# vtk DataFile Version 2.0\nMade by PostProcessing (CC-BY-NC)\nASCII \n") ;
  fprintf(out,"DATASET POLYDATA\n") ;
  fprintf(out,"POINTS %d double\n", npts) ;

  for (i=0 ; i<npts ; i++)
   fprintf(out, "%g %g %g\n", pts[i][0], pts[i][1], pts[i][2]) ;
  fprintf(out, "POLYGONS %d %d\n", nfaces, nfaces*3+nfaces) ;
  for (i=0 ; i<nfaces ; i++)
    fprintf(out, "3 %d %d %d\n", faces[i][0], faces[i][1], faces[i][2]) ;

  fprintf(out, "\nCELL_DATA %d\n", nfaces) ;
  fprintf(out, "SCALARS cell_scalars double 1\n") ;
  fprintf(out, "LOOKUP_TABLE default\n") ;
  for (i=0 ; i<nfaces ; i++)
    fprintf(out, "%g\n", data[i]) ;

  fclose(out) ;
  return 0 ; 
}

int Icosahedre::subdivide()
{
 Vector t1, t2, t3, newpt ; int i,j ;
 int p[3], e[2] ; int lastpt=npts, lastedge=nedges, lastface=nfaces ;
 int *milieux ;

 milieux=(int*)malloc(sizeof(int)*nedges) ;
 for (i=0 ; i<nedges ; i++) milieux[i]=-1 ;

 pts=(double**)realloc(pts, sizeof(double*)*(npts+nedges)) ;
 for (i=npts ; i<npts+nedges ; i++) pts[i]=(double*)malloc(sizeof(double)*3) ;
 edges=(int**)realloc(edges, sizeof(double*)*(2*nedges+3*nfaces)) ;
 for (i=nedges ; i<2*nedges+3*nfaces ; i++) edges[i]=(int*)malloc(sizeof(int)*2) ;
 faces=(int**)realloc(faces, sizeof(double*)*(4*nfaces)) ;
 for (i=nfaces ; i<4*nfaces ; i++) faces[i]=(int*)malloc(sizeof(int)*3) ;

 for (i=0 ; i<nfaces ; i++)
 {
   t1.set(pts[faces[i][0]][0], pts[faces[i][0]][1], pts[faces[i][0]][2]) ;
   t2.set(pts[faces[i][1]][0], pts[faces[i][1]][1], pts[faces[i][1]][2]) ;
   t3.set(pts[faces[i][2]][0], pts[faces[i][2]][1], pts[faces[i][2]][2]) ;

   for (j=0 ; j<nedges ; j++)
   {
     if ((faces[i][0]==edges[j][0] && faces[i][1]==edges[j][1]) || (faces[i][0]==edges[j][1] && faces[i][1]==edges[j][0]))
       e[0]=j ;
     if ((faces[i][1]==edges[j][0] && faces[i][2]==edges[j][1]) || (faces[i][1]==edges[j][1] && faces[i][2]==edges[j][0]))
       e[1]=j ;
     if ((faces[i][2]==edges[j][0] && faces[i][0]==edges[j][1]) || (faces[i][2]==edges[j][1] && faces[i][0]==edges[j][0]))
       e[2]=j ;
   }

   if (milieux[e[0]]==-1)
   {
     newpt=t1+(t2-t1)/2. ; newpt=newpt/newpt.norm() ;
     pts[lastpt][0]=newpt(1) ; pts[lastpt][1]=newpt(2) ; pts[lastpt][2]=newpt(3) ;
     milieux[e[0]]=lastpt ;
     lastpt++ ;
   }
   if (milieux[e[1]]==-1)
   {
     newpt=t2+(t3-t2)/2. ; newpt=newpt/newpt.norm() ;
     pts[lastpt][0]=newpt(1) ; pts[lastpt][1]=newpt(2) ; pts[lastpt][2]=newpt(3) ;
     milieux[e[1]]=lastpt ;
     lastpt++ ;
   }
   if (milieux[e[2]]==-1)
   {
     newpt=t3+(t1-t3)/2. ; newpt=newpt/newpt.norm() ;
     pts[lastpt][0]=newpt(1) ; pts[lastpt][1]=newpt(2) ; pts[lastpt][2]=newpt(3) ;
     milieux[e[2]]=lastpt ;
     lastpt++ ;
   }

   //faces[i][1]=milieux[e[0]] ; faces[i][2]=milieux[e[3]] ;
   faces[lastface][0]=faces[i][0] ;   faces[lastface][1]=milieux[e[0]] ; faces[lastface][2]=milieux[e[2]] ; lastface++ ;
   faces[lastface][0]=milieux[e[0]] ; faces[lastface][1]=faces[i][1] ;   faces[lastface][2]=milieux[e[1]] ; lastface++ ;
   faces[lastface][0]=milieux[e[2]] ; faces[lastface][1]=milieux[e[1]] ; faces[lastface][2]=faces[i][2] ; lastface++ ;
   faces[i][0]=milieux[e[0]] ; faces[i][1]=milieux[e[1]] ; faces[i][2]=milieux[e[2]] ;

   edges[lastedge][0]=milieux[e[0]] ; edges[lastedge][1]=milieux[e[1]] ; lastedge++ ;
   edges[lastedge][0]=milieux[e[1]] ; edges[lastedge][1]=milieux[e[2]] ; lastedge++ ;
   edges[lastedge][0]=milieux[e[2]] ; edges[lastedge][1]=milieux[e[0]] ; lastedge++ ;

 }

 for (i=0 ; i<nedges ; i++)
 {
  edges[lastedge][0]=edges[i][1] ;
  edges[i][1]=milieux[i] ;
  edges[lastedge][1]=milieux[i] ;
  lastedge++ ;
 }

 npts=npts+nedges ;
 nedges=2*nedges+3*nfaces ;
 nfaces*=4 ;
 printf("Icosahedre subdivision : %d %d %d\n", lastpt , lastface, lastedge) ;

 data=(double*)realloc(data, sizeof(double)*nfaces) ;
 solidangle=(double*)realloc(solidangle, sizeof(double)*nfaces) ;
 for (i=nfaces/4 ; i<nfaces ; i++) data[i]=0 ;
 for (i=0 ; i<nfaces ; i++)
 {
  t1.set(pts[faces[i][0]][0], pts[faces[i][0]][1], pts[faces[i][0]][2]) ;
  t2.set(pts[faces[i][1]][0], pts[faces[i][1]][1], pts[faces[i][1]][2]) ;
  t3.set(pts[faces[i][2]][0], pts[faces[i][2]][1], pts[faces[i][2]][2]) ;
  solidangle[i]=Geometrie::solid_angle(t1,t2,t3) ;
 }

 free(milieux) ;

 // TEST
 /*double sum=0 ;
 for (i=0 ; i<nfaces ; i++)
 {
   t1.set(pts[faces[i][0]][0], pts[faces[i][0]][1], pts[faces[i][0]][2]) ;
   t2.set(pts[faces[i][1]][0], pts[faces[i][1]][1], pts[faces[i][1]][2]) ;
   t3.set(pts[faces[i][2]][0], pts[faces[i][2]][1], pts[faces[i][2]][2]) ;
   printf("[%g %g]", Geometrie::tri_surface(t1,t2,t3), solidangle[i]) ;
   sum+=solidangle[i];
 }
 printf("{{{%g}}}", sum) ; */
 // END TEST

 return 0 ;
}

int Icosahedre::belonging_tri (Vector V)
{
  int i ;
  Vector t1, t2, t3 ;

  for (i=0 ; i<nfaces ; i++)
  {
   t1.set(pts[faces[i][0]][0], pts[faces[i][0]][1], pts[faces[i][0]][2]) ;
   t2.set(pts[faces[i][1]][0], pts[faces[i][1]][1], pts[faces[i][1]][2]) ;
   t3.set(pts[faces[i][2]][0], pts[faces[i][2]][1], pts[faces[i][2]][2]) ;
   if (Fonction::intersect3D_RayTriangle(V, t1, t2, t3)==1)
     return i ;
  }

  return -1 ;
}

int Icosahedre::deform (void)
{
int i,j,k ; double max ; double factor ; int n ;

/*max=0 ;
for (i=0 ; i<nfaces ; i++)
  if (data[i]>max)
    max=data[i] ;
*/
max=1 ;

for (i=0 ; i<npts ; i++)
{
  factor=0 ; n=0 ;
  for (j=0 ; j<nfaces ; j++)
  {for (k=0 ; k<3 ; k++)
    {
      if (faces[j][k]==i)
      {factor+=data[j] ; n++ ; }
    }
  }
  if (n!=5 && n!=6) {DISP_Warn("WARN: nombre de triangle auxquels appartient le point incorrect (5 attendu) dans Icosahedre:deform.") ; printf("%d\n", n) ; }
  else
  {
    factor=factor/n ;
    pts[i][0]*=(factor/max) ;
    pts[i][1]*=(factor/max) ;
    pts[i][2]*=(factor/max) ;
  }
}

return 1 ;
}
