#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <malloc.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include "date2mjd_ld.c"
#include "ascii_header.c"
#include <fftw.h>

#define DADAHDR_SIZE 4096

void usage(char *prg_name)
{
  fprintf(stdout,
           "%s [options]\n"
           " --xo  First pol0 file base\n"
	   " --xe  Second pol0 file base\n"
	   " --yo  First pol1 file base\n"
	   " --ye  Second pol1 file base\n"
	   " -T    Starting UT in Yr-Mon-Dat-hr:min:sec (e.g., 2017-01-01-01:03:05)\n" 
           " -S    Dada header file\n"
	   " -f    Central frequency (MHz)\n"
	   " -b    Bandwidth (MHz)\n"
	   " -s    Band sense (-1 for lower, 1 for upper, by default 1)\n"
	   " -N    Source name\n"
	   " -i    Start index (by default 0)\n"
	   " -R    RA (by default 00:00:00.00)\n"
	   " -D    Dec (by default -00:00:00.00)\n"
           " -O    Route for output\n"
	   " -h    Available options\n"
	   "\n"
	   " -c    Enable cutting option\n"
	   " -l    Lower end to keep\n"
	   " -u    Upper end to keep\n",
	  prg_name);
  exit(0);
}

main(int argc, char *argv[])
{
  FILE *bb[4],*dadahdr,*odada;
  int arg,len,ibg,i,j,f,ied,ndim,fct,nblk,k,npol,ctblk,bs,ct,optct,flowidx,fhighidx;
  char hdrbuff[DADAHDR_SIZE],oroute[1024],bbbase[4][1024],bbname[4][1024],dadaname[1024],hdrname[1024],ut[32],srcname[1024],dat,mjd[64],*dblk,*blkp0,*blkp1,ra[64],dec[64];
  float freq,bw,flow,fhigh;
  double ts;
  long fsize,sampct,nblkout,UDPsize;
  struct stat filestat;
  struct option longopts[]={
    {"xo", required_argument, NULL, 'W'},
    {"xe", required_argument, NULL, 'X'},
    {"yo", required_argument, NULL, 'Y'},
    {"ye", required_argument, NULL, 'Z'},
    {0, 0, 0, 0}
  };

  optct=0;
  ct=0;
  bs=1;
  freq=-999.999;
  bw=-999.999;
  ndim=1;
  ibg=1;
  ied=-1;
  sampct=0;
  fct=0;
  len=1;
  nblk=4096;
  npol=2;
  ctblk=0;
  flowidx=-1;
  fhighidx=-1;
  //Size of one UDP file
  UDPsize=2147483648;
  // Number of unloaded blks from each UDP file to each dada file (close to 10s)
  //nblkout=2197266;
  nblkout=524288;
  blkp0=(char *)malloc(sizeof(char)*nblk);
  blkp1=(char *)malloc(sizeof(char)*nblk);
  dblk=(char *)malloc(sizeof(char)*nblk*npol);
  strcpy(ra,"00:00:00.00");
  strcpy(dec,"-00:00:00.00");
  memset(hdrbuff,0,DADAHDR_SIZE);
  strcpy(srcname,"Not given");
  strcpy(oroute,"Not given");

  if(argc==1)
    {
      usage(argv[0]);
      exit(0);
    }

  while((arg=getopt_long(argc,argv,"hS:f:b:O:T:N:i:R:D:s:cl:u:",longopts,NULL)) != -1)
    {
      switch(arg)
	{
	case 'W':
	  strcpy(bbbase[0],optarg);
	  break;

	case 'X':
	  strcpy(bbbase[1],optarg);
	  break;

	case 'Y':
	  strcpy(bbbase[2],optarg);
	  break;

	case 'Z':
	  strcpy(bbbase[3],optarg);
	  break;

	case 'S':
	  strcpy(hdrname,optarg);
	  break;

	case 'f':
	  freq=atof(optarg);
	  break;

	case 'b':
	  bw=atof(optarg);
	  break;

	case 'O':
	  strcpy(oroute,optarg);
	  break;

	case 'T':
	  strcpy(ut,optarg);
	  break;

	case 'N':
	  strcpy(srcname,optarg);
	  break;

	case 'i':
	  ibg=atoi(optarg);
	  break;

	case 'h':
	  usage(argv[0]);
	  return 0;

	case 'R':
	  strcpy(ra,optarg);
	  break;

	case 'D':
	  strcpy(dec,optarg);
	  break;

	case 's':
	  bs=atoi(optarg);
	  break;

	case 'c':
	  optct=1;
	  break;

	case 'l':
	  flow=atof(optarg);
	  break;

	case 'u':
	  fhigh=atof(optarg);
	  break;

	case 0:
	  break;

	default:
	  usage(argv[0]);
	  return 0;
	}
    }
    
  printf("Check provided option values...\n");
  // Check provided info
  if(freq<0)
    {
      fprintf(stderr,"Error: Wrong or non-given central frequency.\n");
      exit(0);
    }
  if(bw==-999.999)
    {
      fprintf(stderr,"Error: Wrong or non-given bandwidth.\n");
      exit(0);
    }
  f=open(hdrname,O_RDONLY);
  if(fstat(f,&filestat)<0)
    {
      fprintf(stderr,"Error: Dada header %s cannot access.\n");
      exit(0);
    }
  close(f);
  if(strcmp(srcname,"Not given")==0)
    {
      fprintf(stderr,"Error: Source name not given.\n");
      exit(0);
    }
  if(strcmp(oroute,"Not given")==0)
    {
      fprintf(stderr,"Error: unload route not given.\n");
      exit(0);
    }
  if(strcmp(ut,"Not given")==0)
    {
      fprintf(stderr,"Error: starting UT not given.\n");
      exit(0);
    }
  if(optct == 1 && (flowidx<0 || fhighidx<0))
    {
      fprintf(stderr,"Error: cutting edge not given.\n");
      exit(0);
    }


  // Get MJD from given date 
  date2mjd_ld(ut,mjd);
  printf("MJD: %s\n",mjd);

  // Sampling interval
  ts=1.0/(double)bw/2*ndim;
  
  // Size of an unload file
  fsize=4*nblk*nblkout;
  printf("File size: %ld\n",fsize);

  // Re-calculate output parameters when cutting
  if(optct == 1)
    {
      printf("This mode is not implemented yet.\n");
      exit(0);
    }

  // Read DADA header template
  dadahdr=fopen(hdrname,"rt");
  if(dadahdr==NULL)
    {
      fprintf(stderr,"Error: Cannot open header template file %s.\n",hdrname);
      exit(0);
    }
  fread(hdrbuff,1,DADAHDR_SIZE,dadahdr);
  fclose(dadahdr);

  // Prepare global variables in header
  ascii_header_set(hdrbuff,"UTC_START","%s",ut);
  ascii_header_set(hdrbuff,"SOURCE","%s",srcname);
  ascii_header_set(hdrbuff,"RA","%s",ra);
  ascii_header_set(hdrbuff,"DEC","%s",dec);
  ascii_header_set(hdrbuff,"FREQ","%f",freq);
  ascii_header_set(hdrbuff,"BW","%f",bw*bs);
  ascii_header_set(hdrbuff,"TSAMP","%.18lf",ts);
  ascii_header_set(hdrbuff,"MJD_START","%s",mjd);

  // Specified header values for the first file
  sprintf(dadaname,"%s/%s_%.0f_%016ld.000000.dada",oroute,ut,freq,fsize*fct+UDPsize*(ibg-1)*4);
  ascii_header_set(hdrbuff,"FILE_SIZE","%ld",fsize);
  ascii_header_set(hdrbuff,"FILE_NAME","%s",dadaname);
  ascii_header_set(hdrbuff,"OBS_OFFSET","%ld",fsize*fct+UDPsize*(ibg-1)*4);

  // Check file existence
  printf("Checking available files...\n");
  for(j=0;;j++)
    {
      // Get file status
      for(i=0;i<4;i++)
	{
	  sprintf(bbname[i],"%s_%04i.dat",bbbase[i],ibg+j);
	  f=open(bbname[i],O_RDONLY);
	  if(fstat(f,&filestat)<0) 
	    {
	      ied=ibg+j-1;
	      break;
	    }
	  close(f);
	}
      if(ied>0) break;
    }
  printf("Index starts: %i; Index ends: %i\n",ibg,ied);
  
  odada=fopen(dadaname,"wb");
  fwrite(hdrbuff,1,DADAHDR_SIZE,odada);

  printf("Start data conversion...\n");
  // Main loop
  for(j=ibg;j<=ied;j++)
    {
      printf("Working on index %i...\n",j);
      // Open UDP files
      for(i=0;i<4;i++)
        {
          sprintf(bbname[i],"%s_%04i.dat",bbbase[i],j);
	  bb[i]=fopen(bbname[i],"rb");
	  if(bb[i]==NULL)
	    {
	      fprintf(stderr,"Error: Cannot open %s.\n",bbname[i]);
	      exit(0);
	    }
	}
      do{
	// Copy samples
	fread(blkp0,sizeof(char)*nblk,1,bb[0]);
	fread(blkp1,sizeof(char)*nblk,1,bb[2]);
	if(feof(bb[0])==1 || feof(bb[1])==1 || feof(bb[2])==1 || feof(bb[3])==1) break;
	for(k=0;k<nblk;k++)
	  {
	    dblk[k*npol]=blkp0[k];
	    dblk[k*npol+1]=blkp1[k];
	  }
	fwrite(dblk,sizeof(char)*nblk*npol,1,odada);
	sampct+=nblk;

	// Not implemented below
	// FFT 
	// Select useful channels
	// FFT back
	// Get mean & rms
	// Redigitize
	// Write to dadafile

	fread(blkp0,sizeof(char)*nblk,1,bb[1]);
	fread(blkp1,sizeof(char)*nblk,1,bb[3]);
	for(k=0;k<nblk;k++)
	  {
	    dblk[k*npol]=blkp0[k];
	    dblk[k*npol+1]=blkp1[k];
	  }
	sampct+=nblk;
	fwrite(dblk,sizeof(char)*nblk*npol,1,odada);

	ctblk++;		  
	ct++;
	// If end dada, close and open
	if(ctblk==nblkout)
	  {
	    // Close written file
	    fclose(odada);
	    fct++;
	    printf("%s unloaded.\n",dadaname);

	    // Prepare for the next file
	    //sprintf(dadaname,"%s/%s_%.0f_%016ld.000000.dada",oroute,ut,freq,fsize*fct);
	    sprintf(dadaname,"%s/%s_%.0f_%016ld.000000.dada",oroute,ut,freq,fsize*fct+UDPsize*(ibg-1)*4);
	    ascii_header_set(hdrbuff,"FILE_SIZE","%ld",fsize);
	    ascii_header_set(hdrbuff,"FILE_NAME","%s",dadaname);
	    ascii_header_set(hdrbuff,"OBS_OFFSET","%ld",fsize*fct+UDPsize*(ibg-1)*4);
	    odada=fopen(dadaname,"wb");
	    fwrite(hdrbuff,1,DADAHDR_SIZE,odada);

	    // Reset sample count
	    sampct=0;
	    ctblk=0;
	  }
      }while(1);

      printf("Read %i\n",ct);
      // Close up
      for(i=0;i<4;i++)
	fclose(bb[i]);
      printf("Index %i finished.\n",j);
    }
  fclose(odada);
  printf("%s unloaded.\n",dadaname);
}
