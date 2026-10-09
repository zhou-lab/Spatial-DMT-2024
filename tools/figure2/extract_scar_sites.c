/* Extract covered CpG ordinal/M/U triples for selected physical YAME records.
 * Binary output: uint64 record, covered count, then covered * (ordinal,M,U).
 * Ordinal follows format-3 zero runs; no genomic coordinate remapping.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
static uint64_t le(unsigned char *p,int n){uint64_t v=0;for(int j=0;j<n;j++)v|=(uint64_t)p[j]<<(8*j);return v;}
static void die(char*s){fprintf(stderr,"%s\n",s);exit(1);}
int main(int argc,char**argv){
 if(argc!=3)die("cg path and record list required");
 int keep[20000]={0},r;FILE *ids=fopen(argv[2],"r");if(!ids)die("record list missing");
 while(fscanf(ids,"%d",&r)==1){if(r<0||r>=20000)die("record outside limit");keep[r]=1;}fclose(ids);
 gzFile f=gzopen(argv[1],"rb");if(!f)die("cg missing");unsigned char h[17];uint64_t rec=0;int got;
 while((got=gzread(f,h,17))!=0){
  if(got!=17||le(h,8)!=266563789635ULL||h[8]!='3')die("invalid header");
  uint64_t n=le(h+9,8);if(n>1073741824ULL)die("record too large");
  unsigned char*p=malloc(n?n:1);if(!p||gzread(f,p,n)!=(int)n)die("truncated record");
  if(rec<20000&&keep[rec]){
   uint64_t *out=malloc((n+1)*3*sizeof(uint64_t)),k=0,pos=0;
   if(!out)die("allocation failed");
   for(uint64_t i=0;i<n;){int t=p[i]&3,w=t==1?1:t==3?8:2;if(i+w>n)die("truncated value");
    uint64_t v=le(p+i,w),m=0,u=0;i+=w;
    if(!t){pos+=v>>2;continue;}
    if(t==1){m=v>>5;u=(v>>2)&7;}else if(t==2){m=v>>9;u=(v>>2)&127;}else{m=v>>33;u=(v>>2)&2147483647ULL;}
    if(m+u){out[3*k]=pos;out[3*k+1]=m;out[3*k+2]=u;k++;}pos++;
   }
   fwrite(&rec,8,1,stdout);fwrite(&k,8,1,stdout);fwrite(out,24,k,stdout);free(out);
  }free(p);rec++;
 }if(!gzeof(f))die("stream error");gzclose(f);return 0;
}
