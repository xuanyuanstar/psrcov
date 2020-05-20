//Convert vdif format to dada format
//Convert from 2-bit to 8-bit

//This works for vdif output from ALMA: two files for two polarisations, each has 32 channels;
//A complete sample is 16x2=32bit, so a time sample corresponds to two complete samples (CS).
//Channels are organized as below:
//CS1: ch15 ch14 ... ch01 ch00
//CS2: ch31 ch30 ... ch17 ch16

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "hget.c"
#include "mjd2date.h"
#include "ascii_header.c"
#include "cvrt2to8.c"
#include "ran.c"

//Calculate MJD from number of 6-mon counts and seconds
long double get_mjd(int mon, long sec)
{
  int imjd,yr,mjd_ref;
  long double fmjd;
  
  //Reference MJD of 2000-01-01 UTC
  mjd_ref=51544;
  
  //Number of whole years
  yr=floor(mon/2);

  //Add days from whole years
  imjd=mjd_ref+yr*365+floor((yr+3)/4);

  //Add days from the left part of 6-mon count
  imjd+=(mon-yr*2)*(31+28+31+30+31+30);
  if(mon-yr*2==1 && floor(yr/4)*4==yr) imjd++;

  //Add days from seconds
  imjd+=floor(sec/86400);

  fmjd=(sec-floor(sec/86400)*86400)/86400;

  return((long double)imjd+fmjd);
}

int usage(char *prg_name)
{
  fprintf(stdout,
		  "Convert 2-bit vdif output from ALMA vdif to PSRDADA;\n"
		  "Data configured in two files for individual polarisation;\n"
          "Each file is 2-bit, real sampled, 32 freq chans;\n"
		  "A complete sample is 32 bit, so two CS for a time samp;\n"
		  "Data rate 2 Gbyte/s (include 2 pols);\n"
		  "Code to extract a single channel;\n"
		             "%s [options]\n"
		             " -f   Central frequency of the data (MHz)\n"
		             " -l   Frame size in Byte (without header, by default 8000)\n"
		             " -r   Frame header in Byte (by default 32)\n"
            		 " -i   Input vdif file for p1\n"
		             " -j   Input vdif file for p2\n"
		             " -n   Index of frequency channel to extract\n"
		             " -p   Header table for p1\n"
		             " -q   Header table for p2\n"
		             " -D   Dada file header size (by default 4096)\n"
		             " -B   Bytes for one dada file (by default 2500000000,10s)\n"
		             " -S   Sample data header file\n"
		             " -u   L/U side band (0 for lower, 1 for upper)\n"
		             " -s   Number of seconds to get sample statistics (default 10)\n"
		             " -k   Number of seconds to skip from the beginning when getting statistics (default 10)\n"
		             " -O   Route for output \n"
		  " -h   Available options\n",
		  prg_name);
  exit(0);
}

main(int argc, char *argv[])
{
  FILE *invdif[2],*dada,*hdr,*phdr[2];

  //Default dada header file set up
  int DADAHDR_SIZE=4096;
  
  //Default bytes for a single output dada file
  long B_out=2500000000;

  //Default bytes of a frame (without header)
  int len=8000;

  //Default bytes of a frame header
  int fhdr=32;
  
  char ifile[200], jfile[200],oroute[200], hdrfile[200],phdrfile[200],qhdrfile[200],dadahdr[DADAHDR_SIZE],ut[30],mjd_str[25],filename[200],dat;
  unsigned char *inbuffer[2], *outbuffer[2];
  int arg,j_i,j_j,j_q,j_O,j_S,j_p,n_f,n_cs,mon[2],ctoffset,i,j,k,ifreq,nfchan,B_cs,mon_nxt,n_f_s,bs,nf_stat,dati,n_skip;
  float cw,freq,cfreq,ns_stat,s_skip;
  double mean[2],sq,rms[2];
  long double mjd;
  long int idx[2],sec[2],num[2],offset0,sec_nxt,num_nxt,seed[2];
  time_t t;
  
  j_i=0;
  j_j=0;
  j_O=0;
  j_S=0;
  j_p=0;
  j_q=0;
  freq=0.0;
  ctoffset=0;
  ifreq=-1;
  ns_stat=1.0;
  s_skip=10.0;

  //Hard coded ALMA vdif output
  //Number of channel
  nfchan=32;

  //Bytes per complete sample after 2to8 conversion
  B_cs=16;

  //Bandwidth in MHz
  cw=-62.5;
  
  //Read arguments
  while ((arg=getopt(argc,argv,"hf:l:r:i:j:n:p:q:D:B:S:u:O:")) != -1)
	{
	  switch(arg)
		{
		case 'f':
		  freq=atof(optarg);
		  break;

		case 'l':
		  len=atoi(optarg);
		  break;

		case 'r':
		  fhdr=atoi(optarg);
		  break;
		  
		case 'i':
		  strcpy(ifile,optarg);
		  j_i=1;
		  break;

		case 'j':
		  strcpy(jfile,optarg);
		  j_j=1;
		  break;
		  
		case 'n':
		  ifreq=atoi(optarg);
		  break;

		case 'p':
		  strcpy(phdrfile,optarg);
		  j_p=1;
		  break;

		case 'q':
		  strcpy(qhdrfile,optarg);
		  j_q=1;
		  break;
		  
		case 'D':
		  DADAHDR_SIZE=atoi(optarg);
		  break;

		case 'B':
		  B_out=atoi(optarg);
		  break;

		case 'S':
		  strcpy(hdrfile,optarg);
		  j_S=1;
		  break;

		case 'u':
		  bs=atoi(optarg);
		  break;

		case 'O':
		  strcpy(oroute,optarg);
		  j_O=1;
		  break;

		case 's':
		  ns_stat=atof(optarg);
		  break;

		case 'k':
		  s_skip=atof(optarg);
		  break;
		  
		case 'h':
		  usage(argv[0]);
		  return 0;

		default:
		  usage(argv[0]);
		  return 0;
		}
	}
  
  //Check if arguments are enough to procceed
  if(freq==0.0)
	{
	  printf("Missing info of observing frequency.\n");
	  exit(0);
	}
  
  if(j_i==0)
	{
	  printf("No input file for pol1 provided.\n");
	  exit(0);
	}

  if(j_j==0)
	{
	  printf("No input file for pol2 provided.\n");
	  exit(0);
	}
  
  if(j_S==0)
	{
	  printf("No sample of dada header provided.\n");
	  exit(0);
	}

  if(j_O==0)
	{
	  printf("No output route specified.\n");
	  exit(0);
	}

  if(ifreq<0)
	{
	  printf("No/not valid index of frequency channel.\n");
	  exit(0);
	}

  if(j_p==0)
	{
	  printf("No header info provided for pol1.\n");
	  exit(0);
	}

  if(j_q==0)
	{
	  printf("No header info provided for pol2.\n");
	  exit(0);
	}
  if(bs==0)
	{
	}
  else if(bs==1)
	{
	  cw=-cw;
	}
  else
	{
	  printf("Invalid argument for band sense.\n");
	  exit(0);
	}

  //Get seed for random generator
  srand((unsigned)time(&t));
  seed[0]=0-t;
  seed[1]=0-t-2;
  
  //Number of frames per second data in vdif
  //cw x 2 (real sampled) x nfchan x 2 (bits) / 8 (cvt to byte) / len per frame
  n_f_s=fabs(cw)*1000000*2*nfchan*2/8/len;
  printf("Number of frames per second data: %i.\n",n_f_s);
  
  //Number of frames to read to fill a dada file
  //B_out / 4 (2bit in 8 bit out) / 2 (2 pols in 2 files) / (len / nchan (size of time sample per frame) )
  n_f=B_out/4/2/(len/nfchan);
  printf("Number of frames to read to fill a dada file: %i.\n",n_f);

  //Number of frame to skip
  n_skip=n_f_s*s_skip;
  printf("Number of frames to skip from the beginning to get statistics: %i.\n",n_skip);
  
  //Number of frames to get sample statistics
  nf_stat=ns_stat*n_f_s;
  printf("Number of frames to get statistics: %i.\n",nf_stat);
  
  //Number of complete samples in a frame
  n_cs=len/4;
  
  //Allocate memo for one frame
  for(j=0;j<2;j++)
	{
	  inbuffer[j]=malloc(sizeof(unsigned char)*len);
	  outbuffer[j]=malloc(sizeof(unsigned char)*len*4);
	  memset(inbuffer[j],0,len);
	  memset(outbuffer[j],0,len*4);
	}
  
  //Read sample dada header
  memset(dadahdr,0,DADAHDR_SIZE);
  hdr=fopen(hdrfile,"rt");
  fread(dadahdr,1,DADAHDR_SIZE,hdr);
  fclose(hdr);

  //Read headers file to get time
  phdr[0]=fopen(phdrfile,"rt");
  phdr[1]=fopen(qhdrfile,"rt");
  
  fscanf(phdr[0],"%ld %i %ld %ld",&idx[0],&mon[0],&sec[0],&num[0]);
  fscanf(phdr[1],"%ld %i %ld %ld",&idx[1],&mon[1],&sec[1],&num[1]);

  if(mon[0]!=mon[1])
	{
	  printf("Inconsistent observing month for two pols!\n");
	  exit(0);
	}  
  mon_nxt=mon[0];
  
  //Pick the small sec and num for start
  if(sec[0]<sec[1])
	{
	  sec_nxt=sec[0];
	  num_nxt=num[0];
	}
  else if(sec[0]>sec[1])
	{
	  sec_nxt=sec[1];
	  num_nxt=num[1];
	}
  else
	{
	  sec_nxt=sec[0];
	  if(num[0]<=num[1])
		{
		  num_nxt=num[0];
		}
	  else
		{
		  num_nxt=num[1];
		}
	}
  
  //Offset in bytes at the beginning, for dada file
  offset0=2*fabs(cw)*1000000*2/n_f_s*num_nxt;
  printf("Offset at the beginning in bytes: %ld.\n",offset0);
  
  //Get MJD
  mjd=get_mjd(mon_nxt,sec_nxt);

  //Convert MJD to UT
  mjd2date(mjd,ut);
  printf("UT of start: %s\n",ut);

  //write accurate starting time into a string
  sprintf(mjd_str,"%.16Lf",mjd);

  //Caculate central frequency of the channel
  cfreq=freq+cw*((float)ifreq-15.5);

  //Get sample statistics of each polarisation
  printf("Getting sample statistics from %f s of data, after skipping the first %f...\n",ns_stat,s_skip);  
  for(j=0;j<2;j++)
	{
	  mean[j]=0.0;
	  sq=0.0;
	  rms[j]=0.0;
	  if(j==0)
		invdif[j]=fopen(ifile,"rb");
	  else
		invdif[j]=fopen(jfile,"rb");

	  //Skip the first given length of data
	  for(i=0;i<n_skip;i++)
		fseek(invdif[j],fhdr+len,SEEK_CUR);

	  for(i=0;i<nf_stat;i++)
		{
		  //Skip the frame header
		  fseek(invdif[j],fhdr,SEEK_CUR);

		  //Read data
		  fread(inbuffer[j],1,len,invdif[j]);

		  //Expand from 2-bit to 8-bit
		  convert2to8(outbuffer[j], inbuffer[j], len);

		  //Get values
		  for(k=0;k<len*4;k++)
			{
			  dati=(int)outbuffer[j][k];
			  mean[j]+=(double)dati;
			  sq+=pow((double)dati,2.0);
			}
		}
	  fclose(invdif[j]);
	  mean[j]=mean[j]/nf_stat/len/4;
	  rms[j]=sqrt(sq/nf_stat/len/4-pow(mean[j],2.0));
	  printf("Mean: %lf; rms: %lf\n",mean[j],rms[j]);
	}
  
  //Update dada header
  ascii_header_set(dadahdr,"FILE_SIZE","%ld",B_out);
  ascii_header_set(dadahdr,"UTC_START","%s",ut);
  ascii_header_set(dadahdr,"MJD_START","%s",mjd_str);
  ascii_header_set(dadahdr,"FREQ","%f",cfreq);
  ascii_header_set(dadahdr,"TSAMP","%.16lf",1.0/fabs(cw)/2);
  ascii_header_set(dadahdr,"NDIM","%i",1);
  ascii_header_set(dadahdr,"NPOL","%i",2);
  ascii_header_set(dadahdr,"BW","%f",cw);
  
  //Open files
  invdif[0]=fopen(ifile,"rb");
  invdif[1]=fopen(jfile,"rb");

  //Main loop
  while(feof(phdr[0])!=1 && feof(phdr[1])!=1)
	{
	  //Set filename and byte offset
	  sprintf(filename,"%s_%.01f_%016ld.000000.dada",ut,freq,offset0+B_out*ctoffset);
	  ascii_header_set(dadahdr,"FILE_NAME","%s",filename);
	  ascii_header_set(dadahdr,"OBS_OFFSET","%ld",offset0+B_out*ctoffset);

	  //Open output dada file
	  sprintf(filename,"%s/%s_%.01f_%016ld.000000.dada",oroute,ut,freq,offset0+B_out*ctoffset);
	  dada=fopen(filename,"wb");
	  if(dada==NULL)
		{
		  printf("Could not generate output file.\n");
		  exit(1);
		}
	  //Write header
	  fwrite(dadahdr,1,DADAHDR_SIZE,dada);

	  //Loop over to write content
	  for(i=0;i<n_f;i++)
		{
		  //Treat individual pols
		  for(j=0;j<2;j++)
			{
			  //If the available frame matches the time
			  if(mon[j]==mon_nxt && num[j]==num_nxt && sec[j]==sec_nxt)
				{
				  //Move to the right frame and read the data
				  fseek(invdif[j],(len+fhdr)*idx[j]+fhdr,SEEK_SET);
				  fread(inbuffer[j],1,len,invdif[j]);
				  convert2to8(outbuffer[j], inbuffer[j], len);

				  //Get info of the next available frame
				  fscanf(phdr[j],"%ld %i %ld %ld",&idx[j],&mon[j],&sec[j],&num[j]);
				}
			  //Fill noise
			  else
				{
				  printf("Miss available frame at second %ld and number %ld for pol%i. Fill with random noise.\n",sec_nxt,num_nxt,j);
				  for(k=0;k<len*4;k++)
					{
					  dati=(int)roundf((float)mean[j]+gasdev(&seed[j])*(float)rms[j]);
					  if(dati<0) dati=0;
					  if(dati>255) dati=255;
					  outbuffer[j][k]=(unsigned char)dati;
					}
				}
			}

		  //Write data
          for(k=0;k<n_cs/2;k++)
			{
			  for(j=0;j<2;j++)
				{
				  if(ifreq<16)
					{
					  dat=(char)((int)outbuffer[j][2*k*B_cs+B_cs-ifreq]-128);
					}
				  else
					{
					  dat=(char)((int)outbuffer[j][2*k*B_cs+B_cs*3-ifreq]-128);
					}					  
				  fwrite(&dat,1,1,dada);
				}
			}
		  
		  //If the end of table file, break
		  if(feof(phdr[0])==1 || feof(phdr[1])==1) break;
		  		  
		  //Count the next frame to read
		  num_nxt++;
		  if(num_nxt==n_f_s)
			{
			  num_nxt=0;
			  sec_nxt++;
			}
		}
	  //Close output
	  fclose(dada);
	  ctoffset++;
	  printf("%s created.\n",filename);
	}
	
  //Close and clean up
  for(j=0;j<2;j++)
	{
	  fclose(invdif[j]);
	  fclose(phdr[j]);
	  free(inbuffer[j]);
	  free(outbuffer[j]);
	}
}
