// Symmetric simplicial 4-polytopes with symmetric vertex maps (topological degree, n = 4).
// Points: m orbits of Z_n acting on R^4 by R_A = rot(2 pi p1/n) + rot(2 pi q1/n):  A[j*n+i] = R_A^i x_j.
// Polytope = convex hull (brute force; rejected unless simplicial with all V = m n points as vertices).
// Vertex map: C[j*n+i] = R_C^i c_j with R_C = rot(2 pi p2/n) + rot(2 pi q2/n), all pairs (p2,q2) tried.
// degree = sum_t s_t [u in cone(C_t)], s_t = (outward orientation of t) * sign det(C_t); evaluated at NU generic u (must agree).
// usage: sympoly n m p1 q1 seed nA nSamp nClimb sigmaR [startfile]
//   nA = number of point configurations sampled; per new combinatorial type and per (p2,q2): nSamp random maps + nClimb climb moves.
//   startfile (optional): m lines of 4 numbers = base points x_j; then configurations are perturbations of it.
// output: one line per new type "T <hash> F E best p2 q2", and "BEST ..." blocks with base points and base vectors.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#define MAXV 64
#define MAXF 4096
#define NU 3
#define MAXM 8
static int n, m, V, F, fac[MAXF][4], ori[MAXF];
static double A[MAXV][4], C[MAXV][4], U[NU][4];
static uint64_t rs;
static inline double rnd(void){ rs ^= rs<<13; rs ^= rs>>7; rs ^= rs<<17; return (rs>>11)*(1.0/9007199254740992.0); }
static double gauss(void){ double a=rnd(), b=rnd(); if(a<1e-300)a=1e-300; return sqrt(-2*log(a))*cos(6.283185307179586*b); }
static inline double det3(const double*a,const double*b,const double*c,int i,int j,int k){
  return a[i]*(b[j]*c[k]-b[k]*c[j]) - a[j]*(b[i]*c[k]-b[k]*c[i]) + a[k]*(b[i]*c[j]-b[j]*c[i]); }
static inline double det4(const double*a,const double*b,const double*c,const double*d){
  return a[0]*det3(b,c,d,1,2,3) - a[1]*det3(b,c,d,0,2,3) + a[2]*det3(b,c,d,0,1,3) - a[3]*det3(b,c,d,0,1,2); }
static void orbit(double X[][4], const double base[][4], int p, int q){
  for(int j=0;j<m;j++) for(int i=0;i<n;i++){ double t1=6.283185307179586*p*i/n, t2=6.283185307179586*q*i/n; double *x=X[j*n+i]; const double*b=base[j];
    x[0]=cos(t1)*b[0]-sin(t1)*b[1]; x[1]=sin(t1)*b[0]+cos(t1)*b[1]; x[2]=cos(t2)*b[2]-sin(t2)*b[3]; x[3]=sin(t2)*b[2]+cos(t2)*b[3]; } }
static unsigned char seen[1<<21];   // bitset over packed 4-subsets (6 bits per vertex)
static int packed[MAXF];
static inline int pack(int a,int b,int c,int d){ return a|(b<<6)|(c<<12)|(d<<18); }
static void sort4(int*x){ for(int i=0;i<4;i++) for(int j=i+1;j<4;j++) if(x[j]<x[i]){ int t=x[i]; x[i]=x[j]; x[j]=t; } }
// convex hull; returns 1 if simplicial, all points vertices; fills fac, ori, F
static int hull(void){
  double cen[4]={0,0,0,0}; for(int v=0;v<V;v++) for(int k=0;k<4;k++) cen[k]+=A[v][k]/V;
  for(int t=0;t<F;t++) seen[packed[t]>>3]&=~(1<<(packed[t]&7));
  F=0; int isv[MAXV]={0};
  for(int j=0;j<m;j++){ int a=j*n;
    for(int b=0;b<V;b++){ if(b==a) continue; for(int c=b+1;c<V;c++){ if(c==a) continue; for(int d=c+1;d<V;d++){ if(d==a) continue;
      { int q[4]={a,b,c,d}; sort4(q); int pk=pack(q[0],q[1],q[2],q[3]); if(seen[pk>>3]&(1<<(pk&7))) continue; }
      double e1[4],e2[4],e3[4]; for(int k=0;k<4;k++){ e1[k]=A[b][k]-A[a][k]; e2[k]=A[c][k]-A[a][k]; e3[k]=A[d][k]-A[a][k]; }
      double N[4]={ det3(e1,e2,e3,1,2,3), -det3(e1,e2,e3,0,2,3), det3(e1,e2,e3,0,1,3), -det3(e1,e2,e3,0,1,2) };
      double nn=sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]+N[3]*N[3]); if(nn<1e-9) continue;
      int pos=0,neg=0,bad=0;
      for(int v=0;v<V&&!(pos&&neg);v++){ if(v==a||v==b||v==c||v==d) continue; double s=0; for(int k=0;k<4;k++) s+=N[k]*(A[v][k]-A[a][k]); s/=nn;
        if(fabs(s)<1e-7){ bad=1; } if(s>0) pos++; else neg++; }
      if(pos&&neg) continue;
      if(bad) return 0;
      // facet: add its whole orbit
      for(int i=0;i<n;i++){ int q[4]={a,b,c,d}; for(int k=0;k<4;k++) q[k]=(q[k]/n)*n+(q[k]%n+i)%n; sort4(q); int pk=pack(q[0],q[1],q[2],q[3]);
        if(seen[pk>>3]&(1<<(pk&7))) continue; seen[pk>>3]|=1<<(pk&7); if(F>=MAXF) return 0; packed[F]=pk; for(int k=0;k<4;k++){ fac[F][k]=q[k]; isv[q[k]]=1; }
        double f1[4],f2[4],f3[4],f4[4]; for(int k=0;k<4;k++){ f1[k]=A[q[1]][k]-A[q[0]][k]; f2[k]=A[q[2]][k]-A[q[0]][k]; f3[k]=A[q[3]][k]-A[q[0]][k]; f4[k]=cen[k]-A[q[0]][k]; }
        double D=det4(f1,f2,f3,f4); if(fabs(D)<1e-12) return 0; ori[F]=D>0?1:-1; F++; }
    }}}
  }
  for(int v=0;v<V;v++) if(!isv[v]) return 0;
  return 1;
}
static int cmpi(const void*a,const void*b){ return *(const int*)a-*(const int*)b; }
static uint64_t typehash(int *E){
  static int tmp[MAXF]; memcpy(tmp,packed,sizeof(int)*F); qsort(tmp,F,sizeof(int),cmpi);
  uint64_t h=1469598103934665603ULL; for(int t=0;t<F;t++){ h^=(uint64_t)tmp[t]; h*=1099511628211ULL; }
  static unsigned char ad[MAXV][MAXV]; memset(ad,0,sizeof(ad)); int e=0;
  for(int t=0;t<F;t++) for(int i=0;i<4;i++) for(int j=i+1;j<4;j++) if(!ad[fac[t][i]][fac[t][j]]){ ad[fac[t][i]][fac[t][j]]=1; e++; }
  *E=e; return h;
}
// degree of current C; returns 0 on degeneracy, else 1 and *deg
static int degree(int *deg){
  int tot[NU]={0};
  for(int t=0;t<F;t++){ const double *a=C[fac[t][0]], *b=C[fac[t][1]], *c=C[fac[t][2]], *d=C[fac[t][3]];
    double D=det4(a,b,c,d); if(fabs(D)<1e-10) return 0; int sd=D>0?1:-1;
    for(int j=0;j<NU;j++){ const double*u=U[j];
      double d0=det4(u,b,c,d); if(fabs(d0)<1e-12) return 0; if(d0*sd<0) continue;
      double d1=det4(a,u,c,d); if(fabs(d1)<1e-12) return 0; if(d1*sd<0) continue;
      double d2=det4(a,b,u,d); if(fabs(d2)<1e-12) return 0; if(d2*sd<0) continue;
      double d3=det4(a,b,c,u); if(fabs(d3)<1e-12) return 0; if(d3*sd<0) continue;
      tot[j]+=ori[t]*sd; } }
  for(int j=1;j<NU;j++) if(tot[j]!=tot[0]) return 0;
  *deg=tot[0]; return 1;
}
#define HS (1<<20)
static uint64_t hkey[HS]; static int hval[HS];
int main(int argc,char**argv){
  n=atoi(argv[1]); m=atoi(argv[2]); int p1=atoi(argv[3]), q1=atoi(argv[4]); rs=strtoull(argv[5],0,10)*0x9E3779B97F4A7C15ULL+88172645463325252ULL;
  long nA=atol(argv[6]); int nSamp=atoi(argv[7]), nClimb=atoi(argv[8]); double sigR=atof(argv[9]);
  V=n*m; if(V>MAXV||m>MAXM) return 1;
  double start[MAXM][4]; int have=0;
  if(argc>10){ FILE*g=fopen(argv[10],"r"); if(g){ have=1; for(int j=0;j<m;j++) for(int k=0;k<4;k++) if(fscanf(g,"%lf",&start[j][k])!=1) have=0; fclose(g);} }
  for(int i=0;i<20;i++) rnd();
  for(int j=0;j<NU;j++){ double s=0; for(int k=0;k<4;k++){ U[j][k]=gauss(); s+=U[j][k]*U[j][k]; } for(int k=0;k<4;k++) U[j][k]/=sqrt(s); }
  int gbest=0; long ntypes=0, nok=0;
  for(long it=0; it<nA; it++){
    double xb[MAXM][4];
    for(int j=0;j<m;j++){
      if(have){ double sc = (it==0)?0: (rnd()<0.5?0.02:0.15); for(int k=0;k<4;k++) xb[j][k]=start[j][k]+sc*gauss(); }
      else { double s=0; for(int k=0;k<4;k++){ xb[j][k]=gauss(); s+=xb[j][k]*xb[j][k]; } double r=(sigR>0&&rnd()<0.5)?exp(sigR*gauss()):1.0; for(int k=0;k<4;k++) xb[j][k]*=r/sqrt(s); } }
    orbit(A,xb,p1,q1);
    if(!hull()) continue;
    nok++;
    int E; uint64_t h=typehash(&E); int slot=h&(HS-1); int found=0;
    while(hkey[slot]){ if(hkey[slot]==h){ found=1; break; } slot=(slot+1)&(HS-1); }
    if(found) continue;
    hkey[slot]=h; ntypes++;
    int tbest=0, tp=0, tq=0; double tc[MAXM][4];
    for(int p2=0;p2<=n/2;p2++) for(int q2=p2;q2<=n/2;q2++){
      if(p2==0) continue;                     // a trivial summand: the image lies in a half-space unless the fixed part vanishes; skip
      double cb[MAXM][4], cur[MAXM][4]; int curd=0, sgn=1;
      for(int s=0;s<nSamp;s++){ for(int j=0;j<m;j++) for(int k=0;k<4;k++) cb[j][k]=gauss(); orbit(C,cb,p2,q2); int d; if(!degree(&d)) continue;
        if(abs(d)>curd||s==0){ curd=abs(d); sgn=d<0?-1:1; memcpy(cur,cb,sizeof(cb)); } }
      for(int s=0;s<nClimb;s++){ memcpy(cb,cur,sizeof(cb)); int j=(int)(rnd()*m); double r2=rnd(); double sc=r2<0.34?0.03:(r2<0.67?0.15:0.6);
        for(int k=0;k<4;k++) cb[j][k]+=sc*gauss(); orbit(C,cb,p2,q2); int d; if(!degree(&d)) continue; d*=sgn;
        if(d>=curd){ curd=d; memcpy(cur,cb,sizeof(cb)); } }
      if(curd>tbest){ tbest=curd; tp=p2; tq=q2; memcpy(tc,cur,sizeof(cur)); }
    }
    hval[slot]=tbest;
    printf("T %016llx F %d E %d best %d p2 %d q2 %d X",(unsigned long long)h,F,E,tbest,tp,tq);
    for(int j=0;j<m;j++) for(int k=0;k<4;k++) printf(" %.12g",xb[j][k]); printf("\n");
    if(tbest>gbest){ gbest=tbest;
      printf("BEST n %d m %d p1 %d q1 %d F %d E %d deg %d p2 %d q2 %d\n",n,m,p1,q1,F,E,tbest,tp,tq);
      for(int j=0;j<m;j++) printf("X %.17g %.17g %.17g %.17g\n",xb[j][0],xb[j][1],xb[j][2],xb[j][3]);
      for(int j=0;j<m;j++) printf("Y %.17g %.17g %.17g %.17g\n",tc[j][0],tc[j][1],tc[j][2],tc[j][3]);
      fflush(stdout); }
  }
  printf("DONE n %d m %d p1 %d q1 %d configs %ld simplicial %ld types %ld best %d\n",n,m,p1,q1,nA,nok,ntypes,gbest);
  return 0;
}
