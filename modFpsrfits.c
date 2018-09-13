// Modify frequency config in psrfits search mode data

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "psrfits.h"
#include "dec2hms.h"

void usage(char *prg_name)
{
  fprintf(stdout,
		             "%s [options]\n"
		             " -l   Index of channels to extract (low end)\n"
		             " -u   Index of channels to extract (high end)\n"
		             " -f   Frequency scrunch factor (by default 1)\n"
		             " -r   Reverse band sense\n"
		             " -O   Route of the output file(s)\n"
		  " -h   Available options\n",
		  prg_name);
  exit(0);
}

int main(int argc, char *argv[])
{
  FILE *out;
  int fsfac,arg,fbg,fed;

  struct psrfits pf;

  fsfac=1;
  fl=-1;
  fh=-1;
  
  //Read arguments
  while ((arg=getopt(argc,argv,"hfl:fh:f:")) != -1)
	{
	  switch(arg)
		{
		case 'f':
		  fsfac=atoi(optarg);
		  break;

		case 'fl':
		  fl=atoi(optarg);
		  break;
		  
		case 'fh':
		  fh=atof(optarg);
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
  if((fl<0 || fh<0) && (fl != -1 || fh != -1))
	{
	  fprintf(stderr,"ERROR: Invalid channel index provided.\n");
	  exit(0);
	}
  if(fl>fh)
	{
	  fprintf(stderr,"ERROR: Low-end channel index higher than high-end.\n");
	  exit(0);
	}

  //

  
  //Set psrfits main header
  printf("Setting up PSRFITS output...\n");
  pf.filenum = 0;           // This is the crucial one to set to initialize things
  pf.rows_per_file = 200;  // Need to set this based on PSRFITS_MAXFILELEN

  //Set values for our hdrinfo structure
  pf.hdr.scanlen = 86400; // in sec
  strcpy(pf.hdr.observer, "A. Eintein");
  strcpy(pf.hdr.telescope, "ALMA");
  strcpy(pf.hdr.obs_mode, "SEARCH");
  strcpy(pf.hdr.backend, "MARK6");
  strcpy(pf.hdr.source, srcname);
  strcpy(pf.hdr.date_obs, ut);
  strcpy(pf.hdr.poln_type, "LIN");

  //Specify status of output data
  if(dstat=='C')
	strcpy(pf.hdr.poln_order, "AABBCRCI");
  else
	strcpy(pf.hdr.poln_order, "AA+BB");

  strcpy(pf.hdr.track_mode, "TRACK");
  strcpy(pf.hdr.cal_mode, "OFF");
  strcpy(pf.hdr.feed_mode, "FA");
  pf.hdr.dt = spf*tsf/1.0e6;
  pf.hdr.fctr = freq;
  pf.hdr.BW = VDIF_BW*bs;
  pf.hdr.nchan = VDIF_NCHAN;
  pf.hdr.MJD_epoch = mjd;
  pf.hdr.ra2000 = 302.0876876;
  dec2hms(pf.hdr.ra_str, pf.hdr.ra2000/15.0, 0);
  pf.hdr.dec2000 = -3.456987698;
  dec2hms(pf.hdr.dec_str, pf.hdr.dec2000, 1);
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
  pf.hdr.nsblk = 8192;
  pf.hdr.ds_time_fact = 1;
  pf.hdr.ds_freq_fact = 1;
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

  printf("Header prepared. Start to write data...\n");
  adat=fopen("dat.txt","wt");
  
  // Main loop to write subints
  do
	{
	  // Fill time samples in each subint: pf.sub.rawdata
	  for(i=0;i<pf.hdr.nsblk;i++)
		{
		  // Initialize
		  for(k=0;k<VDIF_NCHAN;k++)
			{
			  for(p=0;p<npol;p++)
				sdet[k][p]=0.0;
			} 
		  // Loop over frames
		  for(k=0;k<tsf;k++)
			{
			  //Read frame header
			  fread(vfhdr[0],1,VDIF_HEADER_BYTES,vdif[0]);
			  fread(vfhdr[1],1,VDIF_HEADER_BYTES,vdif[1]);
			  header[0]=(const vdif_header *)vfhdr[0];
			  header[1]=(const vdif_header *)vfhdr[1];

			  //Valid frame
			  if(!getVDIFFrameInvalid(header[0]) && !getVDIFFrameInvalid(header[1]))
				{
				  //Read data in frame
				  fread(buffer[0],1,fbytes,vdif[0]);
				  fread(buffer[1],1,fbytes,vdif[1]);

				  //Get detection
				  getVDIFFrameDetection_32chan(buffer[0],buffer[1],fbytes,det,dstat);
				  
				}
			  //Invalide frame
			  else
				{
				  //Skip the data
				  fseek(vdif[0],fbytes,SEEK_CUR);
				  fseek(vdif[1],fbytes,SEEK_CUR);

				  fprintf(stderr,"Invalid frame detected. Use measured mean & rms to generate fake detection.\n");
				  
				  //Create fake detection with measured mean and rms
				  getVDIFFrameFakeDetection(mean,rms,VDIF_NCHAN,det,seed,fbytes,dstat);
				}
			  
			  //Accumulate value
			  for(j=0;j<VDIF_NCHAN;j++)
				{
				  for(p=0;p<npol;p++)
					sdet[j][p]+=det[j][p];
				}
			  
			  //Break at the end of vdif file
			  if(feof (vdif[0]) || feof (vdif[1])) break;
			}		  
		  //Break when not enough frames to get a sample
		  if(k!=tsf) break;
		  
		  //Write detections in pf.sub.rawdata, in 32-bit float and FPT order (freq, pol, time)
		  for(j=0;j<VDIF_NCHAN;j++)
			{
			  if(npol==4)
				{
				  memcpy(pf.sub.rawdata+i*sizeof(float)*4*VDIF_NCHAN+sizeof(float)*j,&sdet[j][0],sizeof(float));
				  memcpy(pf.sub.rawdata+i*sizeof(float)*4*VDIF_NCHAN+sizeof(float)*VDIF_NCHAN*1+sizeof(float)*j,&sdet[j][1],sizeof(float));
				  memcpy(pf.sub.rawdata+i*sizeof(float)*4*VDIF_NCHAN+sizeof(float)*VDIF_NCHAN*2+sizeof(float)*j,&sdet[j][2],sizeof(float));
				  memcpy(pf.sub.rawdata+i*sizeof(float)*4*VDIF_NCHAN+sizeof(float)*VDIF_NCHAN*3+sizeof(float)*j,&sdet[j][3],sizeof(float));
				}
			  else if(npol==1)
				{
				  memcpy(pf.sub.rawdata+i*sizeof(float)*1*VDIF_NCHAN+sizeof(float)*j,&sdet[j][0],sizeof(float));
				  //if(j==0) fprintf(adat,"%f\n",sdet[j][0]);
				}
			}
		}

	  //Break when subint is not complete
	  if(k!=tsf || i!=pf.hdr.nsblk) break;

	  //Update offset from Start of subint
	  pf.sub.offs = (pf.tot_rows + 0.5) * pf.sub.tsubint;

	  //Write subint
	  psrfits_write_subint(&pf);
	  printf("Subint %i written.\n",pf.sub.tsubint);
	  
	} while(!feof (vdif[0]) && !feof (vdif[1]) && !pf.status && pf.T < pf.hdr.scanlen);
	
  // Close the last file and cleanup
  fits_close_file(pf.fptr, &(pf.status));
  free(pf.sub.dat_freqs);
  free(pf.sub.dat_weights);
  free(pf.sub.dat_offsets);
  free(pf.sub.dat_scales);
  free(pf.sub.rawdata);
  free(buffer[0]);
  free(buffer[1]);
  fclose(vdif[0]);
  fclose(vdif[1]);
  fclose(adat);
  
  printf("Wrote %d subints (%f sec) in %d files.\n",pf.tot_rows, pf.T, pf.filenum);

  return;
}
