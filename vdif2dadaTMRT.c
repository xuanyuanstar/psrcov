//Convert vdif format to dada format

//This works for vdif output from TMRT: Multi-channels and rotate between polarisations and frequencies;
//Configured as 8 x 64 MHz, sampled in 8 bit, real sample, dual polarisation;
//A complete sample is 8 x 4 x 8 = 256 bit, in inverse order of frequency: ch7_l ch7_r ch6_l ch6_r ....

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <assert.h>
#include "hget.c"
#include "mjd2date.c"
#include "ascii_header.c"

//Convert 2-bit string to 8-bit, taken from vdif2to8
//Arguments are output, input, number of bytes in input
static void convert2to8(unsigned char *dest, const unsigned char *src, int bytes)
{
  /* choose levels such that the ratio of high to low is as close to 3.3359                     
   * as possible to best maintain amplitude scaling.  127.5 is the center of                    
   * the scale (equates to 0).  118.5/35.5 is pretty close to optimal.                          
   */
  const unsigned char levels[4] = {9, 92, 163, 246};
  static int first = 1;
  static unsigned char lut2to8[256][4];   /* mapping from input 8 bits (4 samples) to output 4 8-bit samples */
  int i, o;

  if(first)
	{
	  /* assemble look up table */

	  for(i = 0; i < 256; ++i)
		{
		  int j;

		  for(j = 0; j < 4; ++j)
			{
			  int k;

			  k = (i >> (2*j)) & 0x3;
			  lut2to8[i][j] = levels[k];
			}
		}
	}

  o = 0;
  for(i = 0; i < bytes; ++i)
	{
	  dest[o++] = lut2to8[src[i]][0];
	  dest[o++] = lut2to8[src[i]][1];
	  dest[o++] = lut2to8[src[i]][2];
	  dest[o++] = lut2to8[src[i]][3];
	}
}

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
		  "Convert 8-bit vdif output TMRT Mark6 to PSRDADA.\n"
		  "Data configured as 8 x 64 MHz.\n"
		             "%s [options]\n"
		             " -f   Central frequency of the selected channel (MHz)\n"
		             " -l   VDIF frame size in Byte (without header, by default 8000)\n"
		             " -r   VDIF frame header in Byte (by default 32)\n"
            		     " -i   Input vdif file\n"
		             " -n   Index of frequency channel to extract\n"
		             " -p   Header table\n"
		             " -D   DADA file header size (by default 4096)\n"
		             " -B   Bytes for one dada file (by default 2560000000)\n"
		             " -S   Sample data header file\n"
		             " -O   Route for output \n"
		  " -h   Available options\n",
		  prg_name);
  exit(0);
}

main(int argc, char *argv[])
{
  FILE *invdif,*dada,*hdr,*phdr;

  //Default dada header file set up
  int DADAHDR_SIZE=4096;
  
  //Default bytes for a single output dada file
  long B_out = 2560000000;

  //Default bytes of a frame (without header)
  int len=8192;

  //Default bytes of a frame header
  int fhdr=32;
  
  char ifile[200], oroute[200], hdrfile[200],phdrfile[200],dadahdr[DADAHDR_SIZE],ut[30],mjd_str[25],filename[200],dat;
  unsigned char *inbuffer, *outbuffer;
  int arg,j_i,j_O,j_S,j_p,n_f,n_cs,mon,ctoffset,i,j,k,t,ifreq,nfchan,B_cs,mon_nxt,n_f_s, bs, npol, ndim;
  float bw, freq;
  long double mjd;
  long int idx,sec,num,offset0,sec_nxt,num_nxt;
  
  j_i=0;
  j_O=0;
  j_S=0;
  j_p=0;
  freq=0.0;
  ctoffset=0;
  ifreq=-1;

  //Specific for TMRT vdif output
  bw = 64.0;
  nfchan = 8;
  bs = -1;
  npol = 2;
  ndim = 1;
  B_cs = ndim * npol * nfchan;
  
  //Read arguments
  while ((arg=getopt(argc,argv,"hf:l:r:i:n:p:D:B:S:O:")) != -1)
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
		  
		case 'n':
		  ifreq=atoi(optarg);
		  break;

		case 'p':
		  strcpy(phdrfile,optarg);
		  j_p=1;
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

		case 'O':
		  strcpy(oroute,optarg);
		  j_O=1;
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
	  printf("No input file provided.\n");
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

  if(nfchan<1)
	{
	  printf("No/not valid number of frequency channels.\n");
	  exit(0);
	}

  //Number of frames per second data in vdif
  //bw x 2 (real sampled) x nfchan x 2 (pols) x 1 byte
  n_f_s = bw * 1000000 * 2 * nfchan * 2 / len;
  printf("Number of frames per second data: %i.\n",n_f_s);
  
  //Number of frames to read to fill a dada file
  n_f = B_out / (len / nfchan);
  printf("Number of frames to read to fill a dada file: %i.\n",n_f);

  //Number of complete time samples in a frame
  n_cs = len / B_cs;
  
  //Allocate memo for one frame
  inbuffer=malloc(sizeof(unsigned char)*len);
  outbuffer=malloc(sizeof(unsigned char)*len);
  memset(inbuffer,0,len);
  memset(outbuffer,0,len);
  
  //Read sample dada header
  memset(dadahdr,0,DADAHDR_SIZE);
  hdr=fopen(hdrfile,"rt");
  fread(dadahdr,1,DADAHDR_SIZE,hdr);
  fclose(hdr);

  //Read headers file to get time
  phdr=fopen(phdrfile,"rt");
  fscanf(phdr,"%ld %i %ld %ld",&idx,&mon,&sec,&num);

  mon_nxt=mon;
  sec_nxt=sec;
  num_nxt=num;

  //Offset in bytes at the beginning
  offset0=2*bw*1000000*2/n_f_s*num;
  printf("Offset at the beginning in bytes: %ld.\n",offset0);
  
  //Get MJD
  mjd=get_mjd(mon,sec);

  //Convert MJD to UT
  mjd2date(mjd,ut);
  printf("UT of start: %s\n",ut);

  //write accurate starting time into a string
  sprintf(mjd_str,"%.16Lf",mjd);
  
  //Update header
  ascii_header_set(dadahdr,"FILE_SIZE","%ld",B_out);
  ascii_header_set(dadahdr,"UTC_START","%s",ut);
  ascii_header_set(dadahdr,"MJD_START","%s",mjd_str);
  ascii_header_set(dadahdr,"FREQ","%f",freq);
  ascii_header_set(dadahdr,"TSAMP","%.16lf",1.0/bw/2);
  ascii_header_set(dadahdr,"NDIM","%i",ndim);
  ascii_header_set(dadahdr,"NPOL","%i",npol);
  ascii_header_set(dadahdr,"BW","%f",bw*bs);
  
  //Open files
  invdif=fopen(ifile,"rb");

  //Main loop
  while(feof(phdr)!=1)
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

	  //Loop over frames to write content
	  for(i=0;i<n_f;i++)
		{
		  //If the available frame matches the time
		  if(mon==mon_nxt && num==num_nxt && sec==sec_nxt)
			{
			  //Move to the right frame and read the data
			  fseek(invdif,(len+fhdr)*idx+fhdr,SEEK_SET);
			  fread(inbuffer,1,len,invdif);
			  memcpy(outbuffer, inbuffer, len);

			  //Get info of the next available frame
			  fscanf(phdr,"%ld %i %ld %ld",&idx,&mon,&sec,&num);
			}
		  //Fill zeros
		  else
			{
			  printf("Miss available frame for sec %ld and index %ld. Fill zeros.\n",sec_nxt,num_nxt);
			  memset(outbuffer, 0, len);
			}
		
		  //If the end of table file, break
		  if(feof(phdr)==1) break;

		  //Write data
		  for(k=0;k<n_cs;k++)
			{
			  for(j=0; j<npol*ndim; j++)
			    {
			      dat = outbuffer[k * B_cs + npol * ndim * ifreq + j];
			      fwrite(&dat,1,1,dada);
			    }
			}

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
  fclose(invdif);
  fclose(phdr);
  free(inbuffer);
  free(outbuffer);
}

