#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#define NMAX 5000
#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0

#define GET_PSUM for (j=1;j<=ndim;j++) { for (i=1,sum=0.0;i<=mpts;i++)	\
            sum += p[i][j]; psum[j]=sum;}

void amoeba(p,y,ndim,ftol,funk,nfunk)
    float **p,y[],ftol,(*funk)();
int ndim,*nfunk;
{
    int i,j,ilo,ihi,inhi,mpts=ndim+1;
    float ytry,ysave,sum,rtol,amotry(),*psum,*vector();
    void nrerror(),free_vector();

    psum=vector(1,ndim);
    *nfunk=0;
    GET_PSUM
	for (;;) {
	    ilo=1;
	    ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
	    for (i=1;i<=mpts;i++) {
		if (y[i] < y[ilo]) ilo=i;
		if (y[i] > y[ihi]) {
		    inhi=ihi;
		    ihi=i;
		} else if (y[i] > y[inhi])
		    if (i != ihi) inhi=i;
	    }
	    rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
	    if (rtol < ftol) break;
	    if (*nfunk >= NMAX) nrerror("Too many iterations in AMOEBA");
	    ytry=amotry(p,y,psum,ndim,funk,ihi,nfunk,-ALPHA);
	    if (ytry <= y[ilo])
		ytry=amotry(p,y,psum,ndim,funk,ihi,nfunk,GAMMA);
	    else if (ytry >= y[inhi]) {
		ysave=y[ihi];
		ytry=amotry(p,y,psum,ndim,funk,ihi,nfunk,BETA);
		if (ytry >= ysave) {
		    for (i=1;i<=mpts;i++) {
			if (i != ilo) {
			    for (j=1;j<=ndim;j++) {
				psum[j]=0.5*(p[i][j]+p[ilo][j]);
				p[i][j]=psum[j];
			    }
			    y[i]=(*funk)(psum);
			}
		    }
		    *nfunk += ndim;
		    GET_PSUM
			}
	    }
	}
    free_vector(psum,1,ndim);
}

float amotry(p,y,psum,ndim,funk,ihi,nfunk,fac)
    float **p,*y,*psum,(*funk)(),fac;
int ndim,ihi,*nfunk;
{
    int j;
    float fac1,fac2,ytry,*ptry,*vector();
    void nrerror(),free_vector();

    ptry=vector(1,ndim);
    fac1=(1.0-fac)/ndim;
    fac2=fac1-fac;
    for (j=1;j<=ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
    ytry=(*funk)(ptry);
    ++(*nfunk);
    if (ytry < y[ihi]) {
	y[ihi]=ytry;
	for (j=1;j<=ndim;j++) {
	    psum[j] += ptry[j]-p[ihi][j];
	    p[ihi][j]=ptry[j];
	}
    }
    free_vector(ptry,1,ndim);
    return ytry;
}

#undef ALPHA
#undef BETA
#undef GAMMA
#undef NMAX

void spline(x,y,n,yp1,ypn,y2)
    float x[],y[],yp1,ypn,y2[];
int n;
{
    int i,k;
    float p,qn,sig,un,*u,*vector();
    void free_vector();

    u=vector(1,n-1);
    if (yp1 > 0.99e30)
	y2[1]=u[1]=0.0;
    else {
	y2[1] = -0.5;
	u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
    }
    for (i=2;i<=n-1;i++) {
	sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
	p=sig*y2[i-1]+2.0;
	y2[i]=(sig-1.0)/p;
	u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
	u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    if (ypn > 0.99e30)
	qn=un=0.0;
    else {
	qn=0.5;
	un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
    }
    y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
    for (k=n-1;k>=1;k--)
	y2[k]=y2[k]*y2[k+1]+u[k];
    free_vector(u,1,n-1);
}

void splint(xa,ya,y2a,n,x,y)
    float xa[],ya[],y2a[],x,*y;
int n;
{
    int klo,khi,k;
    float h,b,a;
    void nrerror();

    klo=1;
    khi=n;
    while (khi-klo > 1) {
	k=(khi+klo) >> 1;
	if (xa[k] > x) khi=k;
	else klo=k;
    }
    h=xa[khi]-xa[klo];
    if (h == 0.0) nrerror("Bad XA input to routine SPLINT");
    a=(xa[khi]-x)/h;
    b=(x-xa[klo])/h;
    *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

#define FUNC(x) ((*func)(x))

float trapzd(func,a,b,n)
    float a,b;
    float (*func)();  /* ANSI: float (*func)(float); */
int n;
{
    float x,tnm,sum,del;
    static float s;
    static int it;
    int j;

    if (n == 1) {
	it=1;
	return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
    } else {
	tnm=it;
	del=(b-a)/tnm;
	x=a+0.5*del;
	for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
	it *= 2;
	s=0.5*(s+(b-a)*sum/tnm);
	return s;
    }
}

#include <math.h>

#define EPS 0.5e-5
#define JMAX 20
#define JMAXP JMAX+1
#define K 5

float qromb(func,a,b)
    float a,b;
    float (*func)();
{
    float ss,dss,trapzd();
    float s[JMAXP+1],h[JMAXP+1];
    int j;
    void polint(),nrerror();

    h[1]=1.0;
    for (j=1;j<=JMAX;j++) {
	s[j]=trapzd(func,a,b,j);
	if (j >= K) {
	    polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
	    if (fabs(dss) < EPS*fabs(ss)) return ss;
	}
	s[j+1]=s[j];
	h[j+1]=0.25*h[j];
    }
    nrerror("Too many steps in routine QROMB");
    return 0.0;
}

#undef EPS
#undef JMAX
#undef JMAXP
#undef K

#include <math.h>

void polint(xa,ya,n,x,y,dy)
    float xa[],ya[],x,*y,*dy;
int n;
{
    int i,m,ns=1;
    float den,dif,dift,ho,hp,w;
    float *c,*d,*vector();
    void nrerror(),free_vector();

    dif=fabs(x-xa[1]);
    c=vector(1,n);
    d=vector(1,n);
    for (i=1;i<=n;i++) {
	if ( (dift=fabs(x-xa[i])) < dif) {
	    ns=i;
	    dif=dift;
	}
	c[i]=ya[i];
	d[i]=ya[i];
    }
    *y=ya[ns--];
    for (m=1;m<n;m++) {
	for (i=1;i<=n-m;i++) {
	    ho=xa[i]-x;
	    hp=xa[i+m]-x;
	    w=c[i+1]-d[i];
	    if ( (den=ho-hp) == 0.0) nrerror("Error in routine POLINT");
	    den=w/den;
	    d[i]=hp*den;
	    c[i]=ho*den;
	}
	*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
    }
    free_vector(d,1,n);
    free_vector(c,1,n);
}

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(idum)
    int *idum;
{
    static int inext,inextp;
    static long ma[56];
    static int iff=0;
    long mj,mk;
    int i,ii,k;

    if (*idum < 0 || iff == 0) {
	iff=1;
	mj=MSEED-(*idum < 0 ? -*idum : *idum);
	mj %= MBIG;
	ma[55]=mj;
	mk=1;
	for (i=1;i<=54;i++) {
	    ii=(21*i) % 55;
	    ma[ii]=mk;
	    mk=mj-mk;
	    if (mk < MZ) mk += MBIG;
	    mj=ma[ii];
	}
	for (k=1;k<=4;k++)
	    for (i=1;i<=55;i++) {
		ma[i] -= ma[1+(i+30) % 55];
		if (ma[i] < MZ) ma[i] += MBIG;
	    }
	inext=0;
	inextp=31;
	*idum=1;
    }
    if (++inext == 56) inext=1;
    if (++inextp == 56) inextp=1;
    mj=ma[inext]-ma[inextp];
    if (mj < MZ) mj += MBIG;
    ma[inext]=mj;
    return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

#include <math.h>

static double at,bt,ct;
#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ?			\
		     (ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))

static double maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?	\
		  (maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void dsvdcmp(a,m,n,w,v)
    double **a,*w,**v;
    int m,n;
{
    int flag,i,its,j,jj,k,l,nm;
    double c,f,h,s,x,y,z;
    double anorm=0.0,g=0.0,scale=0.0;
    double *rv1,*dvector();
    void nrerror(),free_dvector();

    if (m < n) nrerror("SVDCMP: You must augment A with extra zero rows");
    rv1=dvector(1,n);
    for (i=1;i<=n;i++) {
	l=i+1;
	rv1[i]=scale*g;
	g=s=scale=0.0;
	if (i <= m) {
	    for (k=i;k<=m;k++) scale += fabs(a[k][i]);
	    if (scale) {
		for (k=i;k<=m;k++) {
		    a[k][i] /= scale;
		    s += a[k][i]*a[k][i];
		}
		f=a[i][i];
		g = -SIGN(sqrt(s),f);
		h=f*g-s;
		a[i][i]=f-g;
		if (i != n) {
		    for (j=l;j<=n;j++) {
			for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
			f=s/h;
			for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
		    }
		}
		for (k=i;k<=m;k++) a[k][i] *= scale;
	    }
	}
	w[i]=scale*g;
	g=s=scale=0.0;
	if (i <= m && i != n) {
	    for (k=l;k<=n;k++) scale += fabs(a[i][k]);
	    if (scale) {
		for (k=l;k<=n;k++) {
		    a[i][k] /= scale;
		    s += a[i][k]*a[i][k];
		}
		f=a[i][l];
		g = -SIGN(sqrt(s),f);
		h=f*g-s;
		a[i][l]=f-g;
		for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
		if (i != m) {
		    for (j=l;j<=m;j++) {
			for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
			for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
		    }
		}
		for (k=l;k<=n;k++) a[i][k] *= scale;
	    }
	}
	anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
    }
    for (i=n;i>=1;i--) {
	if (i < n) {
	    if (g) {
		for (j=l;j<=n;j++)
		    v[j][i]=(a[i][j]/a[i][l])/g;
		for (j=l;j<=n;j++) {
		    for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
		    for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
		}
	    }
	    for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
	}
	v[i][i]=1.0;
	g=rv1[i];
	l=i;
    }
    for (i=n;i>=1;i--) {
	l=i+1;
	g=w[i];
	if (i < n)
	    for (j=l;j<=n;j++) a[i][j]=0.0;
	if (g) {
	    g=1.0/g;
	    if (i != n) {
		for (j=l;j<=n;j++) {
		    for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
		    f=(s/a[i][i])*g;
		    for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
		}
	    }
	    for (j=i;j<=m;j++) a[j][i] *= g;
	} else {
	    for (j=i;j<=m;j++) a[j][i]=0.0;
	}
	++a[i][i];
    }
    for (k=n;k>=1;k--) {
	for (its=1;its<=30;its++) {
	    flag=1;
	    for (l=k;l>=1;l--) {
		nm=l-1;
		if (fabs(rv1[l])+anorm == anorm) {
		    flag=0;
		    break;
		}
		if (fabs(w[nm])+anorm == anorm) break;
	    }
	    if (flag) {
		c=0.0;
		s=1.0;
		for (i=l;i<=k;i++) {
		    f=s*rv1[i];
		    if (fabs(f)+anorm != anorm) {
			g=w[i];
			h=PYTHAG(f,g);
			w[i]=h;
			h=1.0/h;
			c=g*h;
			s=(-f*h);
			for (j=1;j<=m;j++) {
			    y=a[j][nm];
			    z=a[j][i];
			    a[j][nm]=y*c+z*s;
			    a[j][i]=z*c-y*s;
			}
		    }
		}
	    }
	    z=w[k];
	    if (l == k) {
		if (z < 0.0) {
		    w[k] = -z;
		    for (j=1;j<=n;j++) v[j][k]=(-v[j][k]);
		}
		break;
	    }
	    if (its == 60) nrerror("No convergence in 60 SVDCMP iterations");
	    x=w[l];
	    nm=k-1;
	    y=w[nm];
	    g=rv1[nm];
	    h=rv1[k];
	    f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
	    g=PYTHAG(f,1.0);
	    f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
	    c=s=1.0;
	    for (j=l;j<=nm;j++) {
		i=j+1;
		g=rv1[i];
		y=w[i];
		h=s*g;
		g=c*g;
		z=PYTHAG(f,h);
		rv1[j]=z;
		c=f/z;
		s=h/z;
		f=x*c+g*s;
		g=g*c-x*s;
		h=y*s;
		y=y*c;
		for (jj=1;jj<=n;jj++) {
		    x=v[jj][j];
		    z=v[jj][i];
		    v[jj][j]=x*c+z*s;
		    v[jj][i]=z*c-x*s;
		}
		z=PYTHAG(f,h);
		w[j]=z;
		if (z) {
		    z=1.0/z;
		    c=f*z;
		    s=h*z;
		}
		f=(c*g)+(s*y);
		x=(c*y)-(s*g);
		for (jj=1;jj<=m;jj++) {
		    y=a[jj][j];
		    z=a[jj][i];
		    a[jj][j]=y*c+z*s;
		    a[jj][i]=z*c-y*s;
		}
	    }
	    rv1[l]=0.0;
	    rv1[k]=f;
	    w[k]=x;
	}
    }
    free_dvector(rv1,1,n);
}

#undef SIGN
#undef MAX
#undef PYTHAG
/* Local Variables:  */
/* mode: c */
/* c-file-style:"stroustrup" */
/* End: */
/* vim: set sw=4 ts=4 et tw=80 smartindent :*/
