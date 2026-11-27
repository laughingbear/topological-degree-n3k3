// Forward RUP (DRUP) proof checker, server mode.  Same checking core as rupcheck.c.
// usage: rupserve base.cnf     (base.cnf: clauses shared by all jobs, here the Grassmann-Pluecker clauses)
// stdin, per job:   J <id>        then one clause per line "l1 l2 ... 0"   (the job's own clauses)
//                   P             then the proof: one lemma per line "l1 ... 0" ; lines starting with d are ignored
//                   E
// stdout, per job:  "<id> VERIFIED" iff the empty clause is derived by RUP steps from base + job clauses, else "<id> FAIL ...".
// Soundness: the database only ever contains base clauses, job clauses, and lemmas that passed the RUP test against it.
// Deletions are ignored (keeping clauses cannot make a wrong lemma pass).  State is fully reset between jobs.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define NV 60000
#define LINE (1<<20)
static int *lits; static long nlits, caplits; static long *cstart; static int *clen; static long ncl, capcl;
static long **watch; static int *nw, *capw; static signed char *val; static int *trail; static int tsz;
static inline int code(int l){ return l>0 ? 2*l : 2*(-l)+1; }
static inline int litval(int l){ int v=val[abs(l)]; return l>0? v : -v; }
static void addwatch(int l,long c){ int k=code(l); if(nw[k]==capw[k]){ capw[k]=capw[k]?2*capw[k]:4; watch[k]=realloc(watch[k],capw[k]*sizeof(long)); } watch[k][nw[k]++]=c; }
static void assign(int l){ val[abs(l)]= l>0?1:-1; trail[tsz++]=l; }
static int propagate(int qhead){
  while(qhead<tsz){
    int l=trail[qhead++]; int fl=-l; int k=code(fl); long *ws=watch[k]; int n=nw[k]; int j=0;
    for(int i=0;i<n;i++){
      long c=ws[i]; int *cl=&lits[cstart[c]]; int len=clen[c];
      if(cl[0]==fl){ cl[0]=cl[1]; cl[1]=fl; }
      if(litval(cl[0])==1){ ws[j++]=c; continue; }
      int found=0;
      for(int t=2;t<len;t++) if(litval(cl[t])!=-1){ cl[1]=cl[t]; cl[t]=fl; addwatch(cl[1],c); found=1; break; }
      if(found){ ws=watch[k]; continue; }
      ws[j++]=c;
      if(litval(cl[0])==-1){ for(i++;i<n;i++) ws[j++]=ws[i]; nw[k]=j; return 1; }
      assign(cl[0]);
    }
    nw[k]=j;
  }
  return 0;
}
static int level0_conflict=0, qhead0=0;
static int tmp[LINE/2];
static void addclause(int *buf,int n){
  int m=0;
  for(int i=0;i<n;i++){ int v=litval(buf[i]); if(v==1) return; if(v==0) tmp[m++]=buf[i]; }
  if(m==0){ level0_conflict=1; return; }
  if(m==1){ assign(tmp[0]); if(propagate(qhead0)) level0_conflict=1; qhead0=tsz; return; }
  if(nlits+m>caplits){ caplits=2*(nlits+m)+1024; lits=realloc(lits,caplits*sizeof(int)); }
  if(ncl==capcl){ capcl=2*capcl+1024; cstart=realloc(cstart,capcl*sizeof(long)); clen=realloc(clen,capcl*sizeof(int)); }
  memcpy(&lits[nlits],tmp,m*sizeof(int)); cstart[ncl]=nlits; clen[ncl]=m; nlits+=m;
  addwatch(tmp[0],ncl); addwatch(tmp[1],ncl); ncl++;
}
static int buf[LINE/2];
static int parse(char *s){ int n=0; char *e; for(;;){ long x=strtol(s,&e,10); if(e==s) return -1; s=e; if(x==0) return n; if(labs(x)>=NV){ fprintf(stderr,"variable out of range\n"); exit(2);} buf[n++]=(int)x; } }
static long ncl0, nlits0;
static void reset(void){
  ncl=ncl0; nlits=nlits0; for(int k=0;k<2*NV+2;k++) nw[k]=0;
  for(long c=0;c<ncl0;c++){ addwatch(lits[cstart[c]],c); addwatch(lits[cstart[c]+1],c); }
  for(int i=0;i<tsz;i++) val[abs(trail[i])]=0; tsz=0; qhead0=0; level0_conflict=0;
}
int main(int argc,char**argv){
  static char line[LINE];
  val=calloc(NV+1,1); trail=malloc((NV+1)*sizeof(int)); watch=calloc(2*NV+2,sizeof(long*)); nw=calloc(2*NV+2,sizeof(int)); capw=calloc(2*NV+2,sizeof(int));
  FILE*f=fopen(argv[1],"r"); if(!f){ fprintf(stderr,"no base file\n"); return 2; }
  while(fgets(line,LINE,f)){ if(line[0]=='c'||line[0]=='p') continue; int n=parse(line); if(n<2){ fprintf(stderr,"base clause too short\n"); return 2; } addclause(buf,n); }
  fclose(f);
  if(tsz!=0||level0_conflict){ fprintf(stderr,"base formula propagates\n"); return 2; }
  ncl0=ncl; nlits0=nlits; fprintf(stderr,"base: %ld clauses\n",ncl0);
  char id[64]="?"; int mode=0, verdict=0; long nlem=0;       // mode 0 idle, 1 clauses, 2 proof ; verdict 1 verified, -1 failed
  while(fgets(line,LINE,stdin)){
    if(line[0]=='J'){ reset(); sscanf(line+1,"%63s",id); mode=1; verdict=0; nlem=0; continue; }
    if(line[0]=='P'){ mode=2; if(level0_conflict) verdict=1; continue; }
    if(line[0]=='E'){ if(verdict==1) printf("%s VERIFIED %ld\n",id,nlem); else printf("%s FAIL %s\n",id,verdict==0?"no empty clause":"lemma not RUP"); fflush(stdout); mode=0; continue; }
    if(line[0]=='d'||line[0]=='c'||line[0]=='\n') continue;
    int n=parse(line); if(n<0){ verdict=-1; continue; }
    if(mode==1){ addclause(buf,n); continue; }
    if(mode!=2||verdict!=0) continue;
    nlem++;
    int sat=0; for(int i=0;i<n;i++) if(litval(buf[i])==1){ sat=1; break; }
    if(!sat){
      int save=tsz; for(int i=0;i<n;i++){ if(litval(buf[i])==0) assign(-buf[i]); }
      int conflict=propagate(save);
      for(int i=tsz-1;i>=save;i--) val[abs(trail[i])]=0; tsz=save;
      if(!conflict){ verdict=-1; continue; }
    }
    if(n==0){ verdict=1; continue; }
    addclause(buf,n); if(level0_conflict) verdict=1;
  }
  return 0;
}
