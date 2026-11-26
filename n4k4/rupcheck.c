// Minimal forward RUP (DRUP) proof checker. Written for this project.
// Usage: rupcheck formula.cnf proof.drup
// Soundness: every clause in the database is a consequence of the formula (original clauses, or
// lemmas verified RUP against the database at the time they were added). Deletions are ignored
// (keeping clauses can only make checks easier, never unsound). Accepts iff the empty clause is derived.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int nvars;
static int *lits; static long nlits, caplits;        // clause literal storage
static long *cstart; static int *clen; static long ncl, capcl;
static long **watch; static int *nw, *capw;           // watch lists indexed by literal code
static signed char *val;                              // val[var] in {-1,0,1}
static int *trail; static int tsz;

static inline int code(int l){ return l>0 ? 2*l : 2*(-l)+1; }
static inline int litval(int l){ int v=val[abs(l)]; return l>0? v : -v; }
static void addwatch(int l, long c){ int k=code(l); if(nw[k]==capw[k]){ capw[k]=capw[k]?2*capw[k]:4; watch[k]=realloc(watch[k],capw[k]*sizeof(long)); } watch[k][nw[k]++]=c; }
static void assign(int l){ val[abs(l)] = l>0?1:-1; trail[tsz++]=l; }

// propagate from trail position qhead; return 1 on conflict
static int propagate(int qhead){
  while(qhead<tsz){
    int l=trail[qhead++]; int fl=-l; int k=code(fl);   // clauses watching fl (now false)
    long *ws=watch[k]; int n=nw[k]; int j=0;
    for(int i=0;i<n;i++){
      long c=ws[i]; int *cl=&lits[cstart[c]]; int len=clen[c];
      if(cl[0]==fl){ cl[0]=cl[1]; cl[1]=fl; }          // ensure cl[1] is the false watch
      if(litval(cl[0])==1){ ws[j++]=c; continue; }
      int found=0;
      for(int t=2;t<len;t++) if(litval(cl[t])!=-1){ cl[1]=cl[t]; cl[t]=fl; addwatch(cl[1],c); found=1; break; }
      if(found) continue;
      ws[j++]=c;
      if(litval(cl[0])==-1){ for(i++;i<n;i++) ws[j++]=ws[i]; nw[k]=j; return 1; }
      assign(cl[0]);
    }
    nw[k]=j;
  }
  return 0;
}
static int level0_conflict=0, qhead0=0;
// add clause to DB (literals in buf). Handles units at level 0.
static void addclause(int *buf, int n){
  // remove level-0 false literals; skip if satisfied at level 0
  int m=0; int tmp[4096];
  for(int i=0;i<n;i++){ int v=litval(buf[i]); if(v==1) return; if(v==0) tmp[m++]=buf[i]; }
  if(m==0){ level0_conflict=1; return; }
  if(m==1){ assign(tmp[0]); if(propagate(qhead0)) level0_conflict=1; qhead0=tsz; return; }
  if(nlits+m>caplits){ caplits=2*(nlits+m)+1024; lits=realloc(lits,caplits*sizeof(int)); }
  if(ncl==capcl){ capcl=2*capcl+1024; cstart=realloc(cstart,capcl*sizeof(long)); clen=realloc(clen,capcl*sizeof(int)); }
  memcpy(&lits[nlits],tmp,m*sizeof(int)); cstart[ncl]=nlits; clen[ncl]=m; nlits+=m;
  addwatch(tmp[0],ncl); addwatch(tmp[1],ncl); ncl++;
}
static int readclause(FILE*f, int *buf, int *isdel){
  int c, n=0, x; *isdel=0;
  for(;;){ c=fgetc(f); if(c==EOF) return -1; if(c=='c'){ while((c=fgetc(f))!='\n'&&c!=EOF); continue; }
    if(c=='p'){ while((c=fgetc(f))!='\n'&&c!=EOF); continue; } if(c=='d'){ *isdel=1; break; }
    if(c=='-'||(c>='0'&&c<='9')){ ungetc(c,f); break; } }
  while(fscanf(f,"%d",&x)==1){ if(x==0) return n; buf[n++]=x; if(abs(x)>nvars){ fprintf(stderr,"var out of range\n"); exit(2);} }
  return -1;
}
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); char line[256]; long ncls0=0;
  while(fgets(line,sizeof line,f)) if(line[0]=='p'){ sscanf(line,"p cnf %d %ld",&nvars,&ncls0); break; }
  val=calloc(nvars+1,1); trail=malloc((nvars+1)*sizeof(int));
  watch=calloc(2*nvars+2,sizeof(long*)); nw=calloc(2*nvars+2,sizeof(int)); capw=calloc(2*nvars+2,sizeof(int));
  int buf[4096], isdel, n; long nread=0;
  while((n=readclause(f,buf,&isdel))>=0){ addclause(buf,n); nread++; if(level0_conflict) break; }
  fclose(f);
  fprintf(stderr,"formula: %d vars, %ld clauses read (header %ld)\n",nvars,nread,ncls0);
  if(level0_conflict){ printf("s VERIFIED (formula trivially UNSAT by unit propagation)\n"); return 0; }
  FILE*g=fopen(argv[2],"r"); long nlem=0, ndel=0;
  while((n=readclause(g,buf,&isdel))>=0){
    if(isdel){ ndel++; continue; }
    nlem++;
    // RUP check: assume negation of lemma
    int sat=0; for(int i=0;i<n;i++) if(litval(buf[i])==1){ sat=1; break; }
    if(!sat){
      int save=tsz, conflict=0;
      for(int i=0;i<n && !conflict;i++){ int v=litval(buf[i]); if(v==0) assign(-buf[i]); }
      conflict=propagate(save);
      for(int i=tsz-1;i>=save;i--) val[abs(trail[i])]=0; tsz=save;
      if(!conflict){ printf("s NOT VERIFIED: lemma %ld is not RUP\n",nlem); return 1; }
    }
    if(n==0){ printf("s VERIFIED: empty clause derived at lemma %ld (deletions ignored: %ld)\n",nlem,ndel); return 0; }
    addclause(buf,n);
    if(level0_conflict){ printf("s VERIFIED: level-0 conflict after lemma %ld (deletions ignored: %ld)\n",nlem,ndel); return 0; }
  }
  printf("s NOT VERIFIED: proof ended without empty clause\n"); return 1;
}
