//convert nuppi format raw data to DADA format
//extract a single channel

#include <stdio.h>
#include "fitshead.h"
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <assert.h>
#include "hget.c"
#include "mjd2date.c"
#include "ascii_header.c"
#include "srcname_corr.c"

int usage(char *prg_name)
{
  fprintf(stdout,
           "%s [options]\n"
           " -N   Maximum size of nuppi file header (by default 131072)\n"
           " -D   Dada file header size (by default 4096)\n"
           " -B   Required bytes for a single output file (by default 640000000)\n"
           " -b   Subband index to extract, begin from 0\n"
           " -L   List of input files\n"
           " -S   Sample data header file\n"
           " -O   Route for output (by default /data2/kliu/tmp/)\n"
	   " -h   Available options\n",
	  prg_name);
  exit(0);
}

main(int argc, char *argv[])
{
  FILE *fraw,*dadahdr_spl,*list,*outdada;

  //default nuppi&dada header file set up
  int MAX_HEADER_SIZE=1024*128;
  int DADAHDR_SIZE=4096;
  char outroute[200];
  strcpy(outroute,"/data2/kliu/tmp/");

  //default requested data bytes for a single output file
  //640000000 for 10s integration time
  long B_out=640000000;

  int bdidx=-1,arg;
  char *listfile=0,*dadahdrsamp=0;

  //read in arguments
  while ((arg=getopt(argc,argv,"hN:D:B:b:L:S:O:")) != -1)
    {
      switch(arg)
	{
	case 'N':
	  MAX_HEADER_SIZE=atoi(optarg);
	  break;

	case 'D':
	  DADAHDR_SIZE=atoi(optarg);
	  break;

	case 'B':
	  B_out=(long)atoi(optarg);
	  break;

	case 'b':
	  bdidx=atoi(optarg);
	  break;

	case 'L':
	  listfile=optarg;
	  break;

	case 'S':
	  dadahdrsamp=optarg;
	  break;

        case 'O':
          strcpy(outroute,optarg);          
          break;

	case 'h':
	  usage(argv[0]);
	  return 0;

	default:
	  usage(argv[0]);
	  return 0;
	}
    }

  //check if provided info is enough
  if(bdidx==-1)
    {
      printf("Extracted band index must be give.\n");
      exit(0);
    }

  if(listfile==0)
    {
      printf("File list must be given.\n");
      exit(0);
    }

  if(dadahdrsamp==0)
    {
      printf("Dada sample file must be given.\n");
      exit(0);
    }

  int hdrlength,blocksize,imjd,smjd,n_bit,n_pol,n_band,blocksize_chan,blocksize_p1,blocksize_p2,ct_block,ct_outfile,i,pktidx_step,pktidx,pktsize,pktidx_pre,overlap,opkt,obyte_chan;
  int jd;
  fpos_t fileposi;
  long fileoffset;
  char hdr_buffer[MAX_HEADER_SIZE],src_name[20],src_name_new[20],ra[15],dec[15],cmd[200],datafilename[200],filebasename[50],outname[50],ut[30],dadahdr[DADAHDR_SIZE],mjd_str[25],*block,*block_p1,*block_p2,buf;
  float freq,t_samp,freq_sub,bw;
  long double mjd,fmjd;

  //initialize counters
  fileoffset=0;
  blocksize_p1=0;
  blocksize_p2=0;
  ct_outfile=0;

  list=fopen(listfile,"rt");
  if(list==NULL)
    {
      printf("List file %s not found.\n",listfile);
      exit(0);
    }

  fscanf(list,"%s",datafilename);

  fraw=fopen(datafilename,"rt");
  if(fraw==NULL)
    {
      printf("Data file %s not found.\n",datafilename);
      exit(0);
    }

  //read data environment from the first datafile

  //initialize header buffer  
  memset(hdr_buffer,0,MAX_HEADER_SIZE);

  //read data of maximum header length                   
  fread(hdr_buffer,1,MAX_HEADER_SIZE,fraw);

  /*-------Obtain header information-------*/
  //get header length                                  
  hdrlength = gethlength(hdr_buffer);

  //Get block size without header   
  hgeti4(hdr_buffer, "BLOCSIZE", &blocksize);

  //Get packet size
  hgeti4(hdr_buffer, "PKTSIZE", &pktsize);

  //Get number of bands 
  hgeti4(hdr_buffer,"OBSNCHAN",&n_band);
  if(bdidx>=n_band)
    {
      printf("Required channel index exceeds number of bands.\n");
      exit(0);
    }

  //Get number of polarisation samples
  hgeti4(hdr_buffer,"NPOL",&n_pol);

  //Get number of sampling bits
  hgeti4(hdr_buffer,"NBITS",&n_bit);

  //Get overlap between blocks, in samples per channel
  hgeti4(hdr_buffer, "OVERLAP", &overlap);

  //Get the packet index of the new block    
  hgeti4(hdr_buffer, "PKTIDX", &pktidx);

  //Number of packets within a block, not overlapped
  obyte_chan=overlap*n_pol*n_bit/8;
  opkt=overlap*n_band*n_pol*n_bit/8/pktsize;
  pktidx_step=blocksize/pktsize-opkt;
  pktidx_pre=-pktidx_step;
  blocksize_chan=blocksize/n_band;
  if(B_out<blocksize_chan)
    {
      printf("Output file size too small.\n");
      exit(0);
    }

  //get source name                                                     
  hgets(hdr_buffer,"SRC_NAME",20,src_name);

  //get starting MJD integer 
  hgeti4(hdr_buffer,"STT_IMJD",&imjd);

  //get starting MJD second                                             
  hgeti4(hdr_buffer,"STT_SMJD",&smjd);

  //convert starting MJD to UT, based on the fact that starting time from an integer second     
  mjd=(long double)imjd+(long double)smjd/86400.0;
  mjd2date(mjd,ut);

  //write accurate starting time into a string                      
  sprintf(mjd_str,"%i.%016.0Lf",imjd,(long double)smjd/86400*1.0e16);

  //get RA                                                              
  hgets(hdr_buffer,"RA_STR",15,ra);

  //get DEC                                                             
  hgets(hdr_buffer,"DEC_STR",15,dec);

  //get central frequency                                               
  hgetr4(hdr_buffer,"OBSFREQ",&freq);

  //get length of full bandwidth
  hgetr4(hdr_buffer,"OBSBW",&bw);

  //get sampling time, and convert into microsecond                     
  hgetr4(hdr_buffer,"TBIN",&t_samp);
  t_samp*=1.0e6;

  //get file basic name without extension                               
  hgets(hdr_buffer,"BASENAME",50,filebasename);

  /*-----------------------------------*/

  //go back to the beginning of data file 
  fseek(fraw, 0, SEEK_SET);

  //calculate central frequency of the output file
  freq_sub=freq-bw/2+bw/n_band/2+bw*bdidx/n_band;

  //read sample DADA format header into memory
  memset(dadahdr,0,DADAHDR_SIZE);
  dadahdr_spl=fopen(dadahdrsamp,"rt");
  fread(dadahdr,1,DADAHDR_SIZE,dadahdr_spl);
  fclose(dadahdr_spl);

  //allocate memo for channel block
  block=malloc((blocksize_chan-obyte_chan)*sizeof(char));

  while(feof(list)==0)
  {
    //calculate expected fileoffset for the next output file
    fileoffset=B_out*ct_outfile;

    //open an output file, file name: UT(obs start)_freq_offset.dada
    sprintf(outname,"%s%s_%.1f_%016ld.000000.dada",outroute,ut,freq_sub,fileoffset);
    outdada=fopen(outname,"wb");

    //modify sample header saved in memory
    ascii_header_set(dadahdr,"FILE_SIZE","%ld",B_out);
    ascii_header_set(dadahdr,"UTC_START","%s",ut);
    ascii_header_set(dadahdr,"MJD_START","%s",mjd_str);
    //modify source name to match dada format requirement
    srcname_corr(src_name,src_name_new);
    ascii_header_set(dadahdr,"SOURCE","%s",src_name_new);
    ascii_header_set(dadahdr,"RA","%s",ra);
    ascii_header_set(dadahdr,"DEC","%s",dec);
    ascii_header_set(dadahdr,"FREQ","%f",freq_sub);
    ascii_header_set(dadahdr,"BW","%f",bw/n_band);
    ascii_header_set(dadahdr,"TSAMP","%f",t_samp);
    ascii_header_set(dadahdr,"FILE_NAME","%s",filebasename);
    ascii_header_set(dadahdr,"OBS_OFFSET","%ld",fileoffset);

    //Write header into output file
    fwrite(dadahdr,1,DADAHDR_SIZE,outdada);

    //If there is left over from previously read block, write in
    if(blocksize_p2!=0) 
      {
	fwrite(block_p2,1,blocksize_p2,outdada);
	free(block_p2);
      }

    //Compute how many integer number of blocks to read next and the leftover
    ct_block=floor((B_out-blocksize_p2)/(blocksize_chan-obyte_chan));
    printf("Integer number of block to read: %i\n",ct_block);
    blocksize_p1=B_out-(blocksize_chan-obyte_chan)*ct_block-blocksize_p2;
    blocksize_p2=blocksize_chan-obyte_chan-blocksize_p1;
    printf("Fraction at the tail: %d\n",blocksize_p1);
    printf("Leftover for the next: %d\n",blocksize_p2);

    //write in data
    for(i=0;i<=ct_block;i++)
      {
	//Read in one more char to see if reach the end of nuppi file
	fgetc(fraw);

	//if the end of the nuppi file, switch to the next
	if(feof(fraw)!=0)
	  {
	    printf("%s finished.\n",datafilename);
	    fclose(fraw);
	    fscanf(list,"%s",datafilename);

	    //If no more data to read
	    if(feof(list)!=0)
	      {
		printf("Data file list finished.\n");
		free(block);
		fclose(list);
		fclose(outdada);
		exit(0);
	      }

	    fraw=fopen(datafilename,"rt");
	    if(fraw==NULL)
	      {
		printf("Data file %s not found.\n",datafilename);
		exit(0);
	      }
	  }
	else
	  //Switch back by one char
	  {
	    fseek(fraw,-sizeof(char),SEEK_CUR);
	  }

	//Read block header into buffer
	memset(hdr_buffer,0,MAX_HEADER_SIZE);
	fread(hdr_buffer,1,MAX_HEADER_SIZE,fraw);

	//Get the packet index of the new block
	hgeti4(hdr_buffer, "PKTIDX", &pktidx);

	//Switch back to the beginning of block
	fseek(fraw, -MAX_HEADER_SIZE, SEEK_CUR);

	//If bloc missed
	if(pktidx!=pktidx_pre+pktidx_step)
	  {
	    printf("%i\n",pktidx);
	    printf("Bloc missed at PKTIDX=%i. Supplement with bloc of zeros.\n",pktidx_pre+pktidx_step);

	    //If the last block to read, write in fractional part and save the left
	    if(i==ct_block)
	      {
		//Allocate memo for part to read and part to leave   
		block_p1=malloc(blocksize_p1*sizeof(char));
		block_p2=malloc(blocksize_p2*sizeof(char));

		//Fill the memo with zero
		memset(block_p1,0,blocksize_p1);
		memset(block_p2,0,blocksize_p2);

		//Write out tail
		fwrite(block_p1,1,blocksize_p1,outdada);
		free(block_p1);
	      }
	    //Write entire bloc
	    else
	      {
		//Fill the memo with zero
		memset(block,0,blocksize_chan-obyte_chan);

		//Write out channel data
		fwrite(block,1,blocksize_chan-obyte_chan,outdada);
	      }
	    pktidx_pre+=pktidx_step;
	  }
	//If bloc not missed
	else{
	  //If the last block to read, write in fractional part and save the left
	  if(i==ct_block)
	    {
	      //Allocate memo for part to read and part to leave
	      block_p1=malloc(blocksize_p1*sizeof(char));
	      block_p2=malloc(blocksize_p2*sizeof(char));

	      //Skip the header and find the write channel
	      fseek(fraw,hdrlength+blocksize_chan*bdidx,SEEK_CUR);

	      //Read in channel data into two parts
	      fread(block_p1,1,blocksize_p1,fraw);
	      fread(block_p2,1,blocksize_p2,fraw);

	      //Switch to the end of the block
	      fseek(fraw,-blocksize_chan*bdidx-(blocksize_chan-obyte_chan)+blocksize,SEEK_CUR);

	      //Write out tail
	      fwrite(block_p1,1,blocksize_p1,outdada);
	      free(block_p1);
	    }
	  //Write in entire block
	  else
	    {
	      //Skip the header and find the write channel
	      fseek(fraw,hdrlength+blocksize_chan*bdidx,SEEK_CUR);

	      //Read in channel data
	      fread(block,1,blocksize_chan-obyte_chan,fraw);

	      //Switch to the end of the block
	      fseek(fraw,-blocksize_chan*bdidx-(blocksize_chan-obyte_chan)+blocksize,SEEK_CUR);

	      //write out channel data
	      fwrite(block,1,blocksize_chan-obyte_chan,outdada);
	    }
	  //Update the previous packet index
	  pktidx_pre+=pktidx_step;
	}
      }
    printf("%s created.\n",outname);
    fclose(outdada);
    ct_outfile++;
  }
}
