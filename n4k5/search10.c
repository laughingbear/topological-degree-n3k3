// max |deg| of vertex maps (K,C): K from pool (oriented 3-spheres on <=9 labels), C in R^4.
// deg = sum over facets f with u in cone(C_f) of eps_f * sign det(C_f).   usage: search pool V seed steps [restart]
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define MAXS 6000
#define NV 10
#define NSUB 210
static int ns, nf[MAXS], fi[MAXS][36], fe[MAXS][36];
static int sub[NSUB][4], nsub;
static double C[NV][4], u[4];
static int T[NSUB], V;
static unsigned long long rs;
static double rnd(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return (rs>>11)*(1.0/9007199254740992.0); }
static double gauss(void){ double a=rnd()+1e-300,b=rnd(); return sqrt(-2*log(a))*cos(6.283185307179586*b); }
static double det4(const double *a,const double *b,const double *c,const double *d){
  double m[4][4]; for(int i=0;i<4;i++){m[i][0]=a[i];m[i][1]=b[i];m[i][2]=c[i];m[i][3]=d[i];}
  double s0=m[0][0]*m[1][1]-m[1][0]*m[0][1], s1=m[0][0]*m[1][2]-m[1][0]*m[0][2], s2=m[0][0]*m[1][3]-m[1][0]*m[0][3];
  double s3=m[0][1]*m[1][2]-m[1][1]*m[0][2], s4=m[0][1]*m[1][3]-m[1][1]*m[0][3], s5=m[0][2]*m[1][3]-m[1][2]*m[0][3];
  double c5=m[2][2]*m[3][3]-m[3][2]*m[2][3], c4=m[2][1]*m[3][3]-m[3][1]*m[2][3], c3=m[2][1]*m[3][2]-m[3][1]*m[2][2];
  double c2=m[2][0]*m[3][3]-m[3][0]*m[2][3], c1=m[2][0]*m[3][2]-m[3][0]*m[2][2], c0=m[2][0]*m[3][1]-m[3][0]*m[2][1];
  return s0*c5-s1*c4+s2*c3+s3*c2-s4*c1+s5*c0; }
static void tab(int only){
  for(int s=0;s<nsub;s++){ int *q=sub[s];
    if(q[3]>=V){T[s]=0;continue;}
    if(only>=0 && q[0]!=only&&q[1]!=only&&q[2]!=only&&q[3]!=only) continue;
    double d=det4(C[q[0]],C[q[1]],C[q[2]],C[q[3]]);
    double d0=det4(u,C[q[1]],C[q[2]],C[q[3]]), d1=det4(C[q[0]],u,C[q[2]],C[q[3]]);
    double d2=det4(C[q[0]],C[q[1]],u,C[q[3]]), d3=det4(C[q[0]],C[q[1]],C[q[2]],u);
    int sg=d>0?1:-1; T[s]=(d0*sg>0&&d1*sg>0&&d2*sg>0&&d3*sg>0)?sg:0; } }
static int best_s;
static int eval(int *per){ int b=-1; for(int s=0;s<ns;s++){ int t=0; for(int f=0;f<nf[s];f++) t+=fe[s][f]*T[fi[s][f]];
    if(t<0)t=-t; if(per && t>per[s]) per[s]=t; if(t>b){b=t;best_s=s;} } return b; }
int main(int argc,char**argv){
  FILE*fp=fopen(argv[1],"r"); V=atoi(argv[2]); rs=88172645463325252ULL^(strtoull(argv[3],0,10)*2654435761ULL); long steps=atol(argv[4]);
  int restart=argc>5?atoi(argv[5]):2000;
  if(fscanf(fp,"%d",&ns)!=1) return 1;
  for(int s=0;s<ns;s++){ if(fscanf(fp,"%d",&nf[s])!=1) return 1; for(int f=0;f<nf[s];f++) if(fscanf(fp,"%d %d",&fi[s][f],&fe[s][f])!=2) return 1; }
  for(int a=0;a<NV;a++)for(int b=a+1;b<NV;b++)for(int c=b+1;c<NV;c++)for(int d=c+1;d<NV;d++){sub[nsub][0]=a;sub[nsub][1]=b;sub[nsub][2]=c;sub[nsub][3]=d;nsub++;}
  static int per[MAXS]; int gbest=-1; double Cb[NV][4], ub[4]; int sb=-1;
  long hist[32]={0}; int cur=-1, stale=0; int Ts[NSUB]; double save[4];
  for(long it=0;it<steps;it++){
    if(cur<0||stale>restart){ for(int i=0;i<NV;i++)for(int j=0;j<4;j++)C[i][j]=gauss(); for(int j=0;j<4;j++)u[j]=gauss(); tab(-1); cur=eval(per); stale=0; hist[cur]++; }
    int v=(int)(rnd()*V); memcpy(save,C[v],sizeof save); memcpy(Ts,T,sizeof T);
    double r=rnd();
    if(r<0.3){ for(int j=0;j<4;j++)C[v][j]=gauss(); }
    else { double sc=r<0.6?0.5:(r<0.85?0.15:0.04); for(int j=0;j<4;j++)C[v][j]+=sc*gauss(); }
    tab(v); int e=eval(per);
    if(e>=cur){ if(e>cur)stale=0; else stale++; cur=e; } else { memcpy(C[v],save,sizeof save); memcpy(T,Ts,sizeof T); stale++; }
    if(cur>gbest){ gbest=cur; eval(0); sb=best_s; memcpy(Cb,C,sizeof C); memcpy(ub,u,sizeof u);
      printf("it %ld new best %d sphere %d\n",it,gbest,sb); fflush(stdout);} }
  printf("BEST %d sphere %d\nC\n",gbest,sb);
  for(int i=0;i<V;i++) printf("%.17g %.17g %.17g %.17g\n",Cb[i][0],Cb[i][1],Cb[i][2],Cb[i][3]);
  printf("u %.17g %.17g %.17g %.17g\n",ub[0],ub[1],ub[2],ub[3]);
  int h[32]={0}; for(int s=0;s<ns;s++) h[per[s]]++; printf("per-sphere max hist:"); for(int i=0;i<32;i++) if(h[i]) printf(" %d:%d",i,h[i]); printf("\n");
  printf("PER"); for(int s=0;s<ns;s++) printf(" %d",per[s]); printf("\n");
  return 0; }
