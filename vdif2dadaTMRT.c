//Convert vdif format to dada format

//This works for vdif output from TMRT: Multi-channels and rotate between polarisations and frequencies;
//Configured as 8 x 64 MHz, sampled in 8 bit, real sample, dual polarisation;
//A complete sample is 8 x 4 x 8 = 256 bit, in inverse order of frequency: ch7_l ch7_r ch6_l ch6_r ....

#define VDIF_HEADER_BYTES       32
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <assert.h>
#include "hget.c"
#include "ascii_header.c"
#include <time.h>
#include <malloc.h>
#include <stdbool.h>
#include <stdint.h>
#include "vdifio.h"
#include "vdif2psrfits.h"
#include "dec2hms.h"

static uint32_t VDIF_BW = 512; //Bandwidth in MHz

int64_t getVDIFFrameOffset(const vdif_header *headerst, const vdif_header *header, uint32_t fps)
{
  uint32_t day[2],sec[2],num[2];
  int64_t offset;

  // Get epoch
  day[0]=getVDIFFrameMJD(headerst);
  day[1]=getVDIFFrameMJD(header);

  // Get second after epoch
  sec[0]=getVDIFFrameSecond(headerst);
  sec[1]=getVDIFFrameSecond(header);

  // Get number of frame after second
  num[0]=getVDIFFrameNumber(headerst);
  num[1]=getVDIFFrameNumber(header);

  offset=(day[1]-day[0])*86400*fps+(sec[1]-sec[0])*fps+(num[1]-num[0]);

  return offset;
}

int getVDIFFrameInvalid_robust(const vdif_header *header, int framebytes, bool ifverbose)
{
  int f, mjd,sec,num,inval;
  f=getVDIFFrameBytes(header);
  mjd=getVDIFFrameMJD(header);
  sec=getVDIFFrameSecond(header);
  num=getVDIFFrameNumber(header);
  inval=getVDIFFrameInvalid(header);
 
  if(f != framebytes || inval || mjd < 50000 || mjd > 60000 || sec < 0 || num < 0 || num > 125000)
    {
      if (ifverbose)
	fprintf(stderr,"Invalid frame: bytes %i mjd %i sec %i num %i inval %i\n",f,mjd,sec,num,inval);
      return 1;
    }
  else
    return 0;
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
  int fbytes = 8192;

  //Default bytes of a frame header
  int fhdr=32;
  
  char ifile[200], oroute[200], hdrfile[200], dadahdr[DADAHDR_SIZE],ut[30],mjd_str[25],filename[200],dat, vfhdr[VDIF_HEADER_BYTES], vfhdrst[VDIF_HEADER_BYTES];
  unsigned char *outbuffer;
  int arg,j_i,j_O,j_S,j_p,n_f,n_cs,mon,ctoffset,i,j,k,t,ifreq,nfchan,B_cs,mon_nxt,n_f_s, bs, npol, ndim, offset, offset_pre;
  float bw, freq;
  long double mjd;
  long int idx,sec,num,offset0,sec_nxt,num_nxt;
  int fps;
  
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
  fbytes = 8192;
  
  //Read arguments
  while ((arg=getopt(argc,argv,"hf:l:r:i:n:D:B:S:O:")) != -1)
	{
	  switch(arg)
		{
		case 'f':
		  freq=atof(optarg);
		  break;

		case 'l':
		 fbytes=atoi(optarg);
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
  n_f_s = bw * 1000000 * 2 * nfchan * 2 / fbytes;
  printf("Number of frames per second data: %i.\n",n_f_s);
  
  //Number of frames to read to fill a dada file
  n_f = B_out / (fbytes / nfchan);
  printf("Number of frames to read to fill a dada file: %i.\n",n_f);

  //Number of complete time samples in a frame
  n_cs = fbytes / B_cs;

  // Number of frames per second
  fps=1000000 * 2 * VDIF_BW * npol / fbytes;
  
  //Allocate memo for one frame
  outbuffer=malloc(sizeof(unsigned char)*fbytes);
  memset(outbuffer,0,fbytes);
  
  //Read sample dada header
  memset(dadahdr,0,DADAHDR_SIZE);
  hdr=fopen(hdrfile,"rt");
  fread(dadahdr,1,DADAHDR_SIZE,hdr);
  fclose(hdr);

  //Open files
  invdif=fopen(ifile,"rb");

  // Move to the first valid frame
  while (true)
    {
      fread(vfhdr, 1, VDIF_HEADER_BYTES, invdif);
      if(!getVDIFFrameInvalid_robust((const vdif_header *)vfhdr, VDIF_HEADER_BYTES+fbytes, false))
	break;
      printf("Invalid frame. Move to the next.\n");
      fseek(invdif, fbytes, SEEK_CUR);
    }

  // Get MJD and UT from the first valid frame header
  mjd = getVDIFFrameDMJD((const vdif_header *)vfhdr, fps);
  sprintf(mjd_str,"%.16Lf",mjd);
  mjd2date(mjd, ut);
  printf("UT of start: %s\n",ut);
  memcpy(vfhdrst,vfhdr,VDIF_HEADER_BYTES);
  offset_pre = -1;
  fseek(invdif, -VDIF_HEADER_BYTES, SEEK_CUR);

  //Offset in bytes at the beginning
  offset0=2*bw*1000000*2/n_f_s*num;
  printf("Offset at the beginning in bytes: %ld.\n",offset0);
  
  //Update dada header
  ascii_header_set(dadahdr,"FILE_SIZE","%ld",B_out);
  ascii_header_set(dadahdr,"UTC_START","%s",ut);
  ascii_header_set(dadahdr,"MJD_START","%s",mjd_str);
  ascii_header_set(dadahdr,"FREQ","%f",freq);
  ascii_header_set(dadahdr,"TSAMP","%.16lf",1.0/bw/2);
  ascii_header_set(dadahdr,"NDIM","%i",ndim);
  ascii_header_set(dadahdr,"NPOL","%i",npol);
  ascii_header_set(dadahdr,"BW","%f",bw*bs);
  
  //Main loop
  while(feof(invdif) != 1)
    {
      // Set filename and byte offset
      sprintf(filename,"%s_%.01f_%016ld.000000.dada",ut,freq,offset0+B_out*ctoffset);
      ascii_header_set(dadahdr,"FILE_NAME","%s",filename);
      ascii_header_set(dadahdr,"OBS_OFFSET","%ld",offset0+B_out*ctoffset);

      // Open output dada file
      sprintf(filename,"%s/%s_%.01f_%016ld.000000.dada",oroute,ut,freq,offset0+B_out*ctoffset);
      dada=fopen(filename,"wb");
      if(dada==NULL)
	{
	  printf("Could not generate output file.\n");
	  exit(1);
	}
      // Write dada header
      fwrite(dadahdr,1,DADAHDR_SIZE,dada);

      // Loop over VDIF frames to write content
      for(i=0;i<n_f;i++)
	{
	  // Find the next valid frame
	  while (true)
	    {
	      fread(vfhdr, 1, VDIF_HEADER_BYTES, invdif);
	      if(!getVDIFFrameInvalid_robust((const vdif_header *)vfhdr, VDIF_HEADER_BYTES+fbytes, false))
		break;
	      printf("Invalid frame. Move to the next.\n");
	      fseek(invdif, fbytes, SEEK_CUR);
	    }

	  // Get frame offset                                                                                                                                
	  offset = getVDIFFrameOffset((const vdif_header *)vfhdrst, (const vdif_header *)vfhdr, fps);

	  // Gap from the last frames, fill in data buffer with zeros
	  if(offset > offset_pre + 1)
	    {
	      printf("Miss available frame. Fill zeros.\n");
	      memset(outbuffer, 0, fbytes);
	      fseek(invdif,  VDIF_HEADER_BYTES, SEEK_CUR);
	    }
	  // Consecutive
	  else
	    fread(outbuffer, 1, fbytes, invdif);

	  offset_pre++;

	  // Write data
	  for(k=0;k<n_cs;k++)
	    {
	      for(j=0; j<npol*ndim; j++)
		{
		  dat = outbuffer[k * B_cs + npol * ndim * ifreq + j];
		  fwrite(&dat,1,1,dada);
		}
	    }
	  if(feof(invdif) == 1)
	    break;
	}

      // Close output
      fclose(dada);
      ctoffset++;
      printf("%s created.\n",filename);
    }
  
  //Close and clean up
  fclose(invdif);
  free(outbuffer);
}

