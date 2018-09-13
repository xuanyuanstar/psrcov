//Convert UDP to psrfits search mode format in 32-bit float
//Each UDP file corresponds to one offload fits file
//FFT length in unit of nblk=4096 in UDP file

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <malloc.h>
#include <stdbool.h>
#include "psrfits.h"

int usage(char *prg_name)
{
  fprintf(stdout,
           "%s [options]\n"
           " --xo First pol0 file base\n"
           " --xe Second pol0 file base\n"
           " --yo First pol1 file base\n"
           " --ye Second pol1 file base\n"
           " -T   Starting UT in Yr-Mon-Dat-hr:min:sec (e.g., 2017-01-01-01:03:05)\n"
           " -f   Central frequency (MHz)\n"
           " -b   Bandwidth (MHz)\n"
           " -N   Source name\n"
	   " -A   RA (by default 00:00:00)\n"
	   " -C   DEC (by default +00:00:00)\n"
           " -i   Start index (by default 1)\n"
	   " -j   End index (by default the last available)\n"
	   " -n   FFT length factor, need to be power of 2 (by default 1) up to 128\n"
	   " -t   Time sample scrunch factor, need to be power of 2 (by default 1)\n"
	   " -S   Source name\n"
	   " -l   Low-end frequency for unload (MHz)\n"
	   " -u   Up-end frequency for unload (MHz)\n"
	   " -s   Band sense (1 for upper, -1 for lower, by default 1)\n"
	   " -D   Ouput data status (I for Stokes I, C for coherence product, X for pol0 I, Y for pol1 I, by default I)\n"
           " -O   Route for output\n"
	   " -h   Available options\n",
          prg_name);
  exit(0);
}

bool IsPowerofTwo(uint x)
{
  return (x&(x-1))==0;
}

int main(int argc, char *argv[])
{
  FILE *bb[4];
  char oroute[1024],bbbase[4][1024],bbname[4][1024],ut[32],srcname[1024],dstat,ra[16],dec[16];
  int arg,ibg,ied,i,j,k,t,s,npol,nchan,bs,nblk,nsub_ed,ncyc,lf_idx,uf_idx,fd,imjd;
  float freq,bw,lf,uf;
  char *bufp0,*bufp1;
  double fmjd;
  long double ts;
  long UDPsize, UDPsize_ed;
  unsigned int tsf,len;
  struct psrfits pf;
  struct stat filestat;
  struct option longopts[]={
    {"xo", required_argument, NULL, 'W'},
    {"xe", required_argument, NULL, 'X'},
    {"yo", required_argument, NULL, 'Y'},
    {"ye", required_argument, NULL, 'Z'},
    {0, 0, 0, 0}
  };

  freq=-999.999;
  bw=-999.999;
  ibg=1;
  ied=-1;
  dstat='I';
  npol=1;
  len=1;
  tsf=1;
  lf=-999.999;
  uf=-999.999;
  bs=1;
  UDPsize=2147483648; // 2 GB
  UDPsize_ed=UDPsize;
  nblk=4096;
  strcpy(ra,"00:00:00");
  strcpy(dec,"+00:00:00");
  strcpy(srcname,"Not given");
  strcpy(oroute,"Not given");
  strcpy(ut,"Not given");

  if(argc==1)
    {
      usage(argv[0]);
      exit(0);
    }
  
  // Read arguments
  while((arg=getopt_long(argc,argv,"hf:b:O:T:N:t:i:j:l:u:s:D:A:C:n:",longopts,NULL)) != -1)
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

	case 'n':
	  len=atoi(optarg);
	  break;

	case 't':
	  tsf=atoi(optarg);
	  break;

        case 'i':
          ibg=atoi(optarg);
          break;

	case 'j':
	  ied=atoi(optarg);
	  break;

	case 'l':
	  lf=atof(optarg);
	  break;

	case 's':
	  bs=atoi(optarg);
	  break;

	case 'u':
	  uf=atof(optarg);
	  break;

	case 'D':
	  strcpy(&dstat,optarg);
	  if(dstat!='C') npol=1;
	  break;

	case 'A':
	  strcpy(ra,optarg);
	  break;

	case 'C':
	  strcpy(dec,optarg);
	  break;

        case 'h':
          usage(argv[0]);
          return 0;

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
  if(bw<0)
    {
      fprintf(stderr,"Error: Wrong or non-given bandwidth.\n");
      exit(0);
    }
  if(strcmp(srcname,"Not given")==0)
    {
      fprintf(stderr,"Error: Source name not given.\n");
      exit(0);
    }
  if(strcmp(oroute,"Not given")==0)
    {
      fprintf(stderr,"Error: Unload route not given.\n");
      exit(0);
    }
  if(strcmp(ut,"Not given")==0)
    {
      fprintf(stderr,"Error: Starting UT not given.\n");
      exit(0);
    }
  if(dstat!='I' && dstat!='C' && dstat!='X' && dstat!='Y')
    {
      fprintf(stderr,"Not recognized status for output data.\n");
      exit(0);
    }
  if(lf<0)
    {
      fprintf(stderr,"Error: Low-end frequency for unload not given.\n");
      exit(0);
    }
  if(uf<0)
    {
      fprintf(stderr,"Error: Upper-end frequency for unload not given.\n");
      exit(0);
    }
  if(lf<freq-bw/2 || lf>=freq+bw/2 || lf>=uf)
    {
      fprintf(stderr,"Error: Invalid low-end frequency for unload.\n");
      exit(0);
    }
  if(uf>freq+bw/2 || uf<=freq-bw/2)
    {
      fprintf(stderr,"Error: Invalid upper-end frequency for unload.\n");
      exit(0);
    }
  if(bs!=1 && bs!=-1)
    {
      fprintf(stderr,"Error: Invalid band sense.\n");
      exit(0);
    }
  if(tsf<1 || IsPowerofTwo(tsf)!=true )
    {
      fprintf(stderr,"Error: Invalid time scrunch factor.\n");
      exit(0);
    }
  if( len<1 || IsPowerofTwo(len)!=true || len > 128 )
  {
    fprintf(stderr,"Error: Invalid FFT length factor %i.\n",len);
    exit(0);
  }

  // Detection blocks
  float det[nblk*len+1][4],sdet[nblk*len+1][4];

  // Number of read cycles per UDP file
  ncyc=UDPsize/(nblk*len)/tsf;

  // Get MJD from given date
  date2mjd(ut,&imjd,&fmjd);

  // Allocate timeseries buffer
  bufp0=(char *)malloc(sizeof(char)*nblk*len*2);
  bufp1=(char *)malloc(sizeof(char)*nblk*len*2);

  // Get frequency config info
  if(bs == 1)
    {
      lf_idx=(int)((lf-(freq-bw/2))/(bw/(nblk*len)));
      uf_idx=(int)((uf-(freq-bw/2))/(bw/(nblk*len)));
    }
  else if(bs == -1)
    {
      lf_idx=(int)((freq+bw/2-uf)/(bw/(nblk*len)));
      uf_idx=(int)((freq+bw/2-lf)/(bw/(nblk*len)));
    }
  nchan=uf_idx-lf_idx+1;
  printf("Channel indices to keep: %i to %i; Total: %i.\n",lf_idx,uf_idx,nchan);

  // Check UDP file existence
  printf("Check through available UDP files to the end...\n");
  if(ied == -1)
    {
      for(j=0;;j++)
	{     
	  // Get file status for odd & even, two pols
	  for(i=0;i<4;i++)
	    {
	      sprintf(bbname[i],"%s_%04i.dat",bbbase[i],ibg+j);
	      fd=open(bbname[i],O_RDONLY);
	      if(fstat(fd,&filestat)<0)
		{
		  ied=ibg+j-1;
		  break;
		}
	      close(fd);
	    }
	  if(ied>0) break;
	}
    }

  printf("UDP Begin: %04i; End: %04i.\n",ibg,ied);

  // Get minimum size for the last available file
  for(i=0;i<4;i++)
    {
      sprintf(bbname[i],"%s_%04i.dat",bbbase[i],ied);
      fd=open(bbname[i],O_RDONLY);
      fstat(fd,&filestat);
      if(filestat.st_size<UDPsize_ed) UDPsize_ed=filestat.st_size;
    }

  // Set psrfits main header
  printf("Setting up PSRFITS output...\n");
  pf.filenum = 0;           // This is the crucial one to set to initialize things

  //Set values for our hdrinfo structure
  pf.hdr.scanlen = 86400; // in sec
  strcpy(pf.hdr.observer, "A. Eintein");
  strcpy(pf.hdr.telescope, "FAST");
  strcpy(pf.hdr.obs_mode, "SEARCH");
  strcpy(pf.hdr.backend, "ROACH2");
  strcpy(pf.hdr.source, srcname);
  strcpy(pf.hdr.date_obs, ut);
  printf("UT:%s\n",ut);
  strcpy(pf.hdr.poln_type, "LIN");

  //Specify status of output data
  if(dstat=='C')
	strcpy(pf.hdr.poln_order, "AABBCRCI");
  else
	strcpy(pf.hdr.poln_order, "AA+BB");

  strcpy(pf.hdr.track_mode, "TRACK");
  strcpy(pf.hdr.cal_mode, "OFF");
  strcpy(pf.hdr.feed_mode, "FA");
  pf.hdr.dt = (double)(nblk*len*tsf)/bw/1.0e6;
  pf.hdr.fctr = freq-bw/2*bs+bw*bs/(nblk*len)*(uf_idx+lf_idx)/2;
  pf.hdr.BW = (uf_idx-lf_idx+1)*bw/(nblk*len)*bs;
  pf.hdr.nchan = nchan;
  pf.hdr.MJD_epoch = (long double)imjd+(long double)fmjd;
  //pf.hdr.ra2000 = 302.0876876;
  //dec2hms(pf.hdr.ra_str, pf.hdr.ra2000/15.0, 0);
  strcpy(pf.hdr.ra_str,ra);
  //pf.hdr.dec2000 = -3.456987698;
  //dec2hms(pf.hdr.dec_str, pf.hdr.dec2000, 1);
  strcpy(pf.hdr.dec_str,dec);
  pf.hdr.azimuth = 123.123;
  pf.hdr.zenith_ang = 23.0;
  pf.hdr.beam_FWHM = 0.25;
  pf.hdr.start_lst = 10000.0;
  pf.hdr.start_sec = 25000.82736876;
  pf.hdr.start_day = 55000;
  pf.hdr.scan_number = 1;
  pf.hdr.rcvr_polns = 2;
  pf.hdr.onlyI = 0;
  pf.hdr.summed_polns = 0;
  pf.hdr.offset_subint = 0;
  pf.hdr.orig_nchan = pf.hdr.nchan;
  pf.hdr.orig_df = pf.hdr.df = pf.hdr.BW / pf.hdr.nchan;
  pf.hdr.nbits = 32;
  pf.hdr.npol = npol;
  pf.hdr.chan_dm = 0.0;
  pf.hdr.fd_hand = 1;
  pf.hdr.fd_sang = 0;
  pf.hdr.fd_xyph = 0;
  pf.hdr.be_phase = 1;
  pf.hdr.nsblk = 1024;
  pf.rows_per_file = ncyc/pf.hdr.nsblk;  // Need to set this based on PSRFITS_MAXFILELEN 
  pf.hdr.ds_time_fact = 1;
  pf.hdr.ds_freq_fact = 1;
  nsub_ed=UDPsize_ed/(nblk*len*tsf)/pf.hdr.nsblk;
  sprintf(pf.basefilename, "%s/%s",oroute,ut);
  
  psrfits_create(&pf);
  
  //Set values for our subint structure
  pf.sub.tsubint = pf.hdr.nsblk * pf.hdr.dt;
  pf.tot_rows = 0.0;
  pf.sub.offs = (pf.tot_rows + 0.5) * pf.sub.tsubint;
  pf.sub.lst = pf.hdr.start_lst;
  pf.sub.ra = pf.hdr.ra2000;
  pf.sub.dec = pf.hdr.dec2000;
  pf.sub.feed_ang = 0.0;
  pf.sub.pos_ang = 0.0;
  pf.sub.par_ang = 0.0;
  pf.sub.tel_az = pf.hdr.azimuth;
  pf.sub.tel_zen = pf.hdr.zenith_ang;
  pf.sub.bytes_per_subint = (pf.hdr.nbits * pf.hdr.nchan * pf.hdr.npol * pf.hdr.nsblk) / 8;
  pf.sub.FITS_typecode = TBYTE;  // 11 = byte      

  // Create and initialize the subint arrays
  pf.sub.dat_freqs = (float *)malloc(sizeof(float) * pf.hdr.nchan);
  pf.sub.dat_weights = (float *)malloc(sizeof(float) * pf.hdr.nchan);
  for (i = 0 ; i < pf.hdr.nchan ; i++)
	{
	  pf.sub.dat_freqs[i] = pf.hdr.fctr - 0.5 * pf.hdr.BW + 0.5 * pf.hdr.df + i * pf.hdr.df;
	  pf.sub.dat_weights[i] = 1.0; 
	}
  pf.sub.dat_offsets = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
  pf.sub.dat_scales = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
  for (i = 0 ; i < pf.hdr.nchan * pf.hdr.npol ; i++)
	{
	  pf.sub.dat_offsets[i] = 0.0;
	  pf.sub.dat_scales[i] = 1.0;
	}
  
  pf.sub.rawdata = (unsigned char *)malloc(pf.sub.bytes_per_subint);

  printf("Header prepared.\n");

  // main loop over UDP files
  for(j=ibg;j<=ied;j++)
    {
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
      printf("Reading from index %i...\n",j);

      // Loop to create subints
      for(k=0;k<pf.rows_per_file;k++)
	{
	  // Loop to detect time samples
	  for(i=0;i<pf.hdr.nsblk;i++)
	    {
	      // Initialize
	      for(t=0;t<nblk*len+1;t++)
		{
		  sdet[t][0]=0.0;
		  sdet[t][1]=0.0;
		  sdet[t][2]=0.0;
		  sdet[t][3]=0.0;
		}
	      // Loop over samples to scrunch
	      for(t=0;t<tsf;t++)
		{
		  // Read sample block(s)
		  for(s=0;s<len;s++)
		    {
		      fread(bufp0+sizeof(char)*nblk*2*s,sizeof(char)*nblk,1,bb[0]);
		      fread(bufp1+sizeof(char)*nblk*2*s,sizeof(char)*nblk,1,bb[2]);
		      fread(bufp0+sizeof(char)*nblk*2*s+sizeof(char)*nblk,sizeof(char)*nblk,1,bb[1]);
		      fread(bufp1+sizeof(char)*nblk*2*s+sizeof(char)*nblk,sizeof(char)*nblk,1,bb[3]);
		    }

		  // Make detection
		  getUDPDetection(bufp0, bufp1, nblk*len*2, det, dstat);

		  // Time scrunch     
		  for(s=0;s<nblk*len+1;s++)
		    {
		      sdet[s][0]+=det[s][0];
		      sdet[s][1]+=det[s][1];
		      sdet[s][2]+=det[s][2];
		      sdet[s][3]+=det[s][3];
		    }
		}
	      // Value sample blk
	      for(t=lf_idx;t<=uf_idx;t++)
		{
		  if(npol==4)
		    {
		      memcpy(pf.sub.rawdata+i*sizeof(float)*4*nchan+sizeof(float)*(t-lf_idx),&sdet[t][0],sizeof(float));
		      memcpy(pf.sub.rawdata+i*sizeof(float)*4*nchan+sizeof(float)*nchan*1+sizeof(float)*(t-lf_idx),&sdet[t][1],sizeof(float));
		      memcpy(pf.sub.rawdata+i*sizeof(float)*4*nchan+sizeof(float)*nchan*2+sizeof(float)*(t-lf_idx),&sdet[t][2],sizeof(float));
		      memcpy(pf.sub.rawdata+i*sizeof(float)*4*nchan+sizeof(float)*nchan*3+sizeof(float)*(t-lf_idx),&sdet[t][3],sizeof(float));
		    }
		  else if(npol==1)
		    {
		      memcpy(pf.sub.rawdata+i*sizeof(float)*1*nchan+sizeof(float)*(t-lf_idx),&sdet[t][0],sizeof(float));
		    }
		}
	    }

	  // Update offset from Start of subint
	  pf.sub.offs = (pf.tot_rows + 0.5) * pf.sub.tsubint;

	  // Write subint
	  psrfits_write_subint(&pf);
	  printf("Subint %i written.\n",pf.sub.tsubint);

          // Break if it is the last subint to write
          if(j==ied && k==nsub_ed-1) break;
	}
      // Close UDP files 
      for(i=0;i<4;i++)
	fclose(bb[i]);
      printf("UDP index %i done.\n",j);
    }

  // Close the last file and cleanup
  fits_close_file(pf.fptr, &(pf.status));
  free(pf.sub.dat_freqs);
  free(pf.sub.dat_weights);
  free(pf.sub.dat_offsets);
  free(pf.sub.dat_scales);
  free(pf.sub.rawdata);
  free(bufp0);
  free(bufp1);

  printf("Wrote %d subints (%f sec) in %d files.\n",pf.tot_rows, pf.T, pf.filenum);

  return;
}
