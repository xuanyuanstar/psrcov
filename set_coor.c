/* fix_psrfits.c */
#include <stdio.h>
#include <stdlib.h>
#include <fitsio.h>
#include <getopt.h>
#include <string.h>

int usage(char *prg_name)
{
  fprintf(stdout,
           "%s [options]\n"
	   " -f   File name\n"
	   " -n   Jname\n"
           " -R   RA in hr:min:sec\n"
           " -D   DEC in hr:min:sec\n"
 	   " -h   Available options\n",
          prg_name);
  exit(0);
}

int main(int argc, char *argv[]) {

  int arg;
  char filename[4096],ra[64],dec[64],jname[64];

  // Read arguments
  if(argc==1)
    {
      usage(argv[0]);
      exit(0);
    }

  while((arg=getopt(argc,argv,"hf:n:R:D:")) != -1)
    {
      switch(arg)
        {
	case 'f':
	  strcpy(filename,optarg);
	  break;

	case 'n':
	  strcpy(jname,optarg);
	  break;

        case 'R':
          strcpy(ra,optarg);
          break;

        case 'D':
          strcpy(dec,optarg);
          break;

	case 0:
          break;

        default:
          usage(argv[0]);
          return 0;
	}
    }

    int status=0;
    fitsfile *f;

    /* open fits file */
    status=0;
    fits_open_file(&f, filename, READWRITE, &status);

    /* Set coordinates */
    fits_movabs_hdu(f, 1, NULL, &status);
    fits_update_key(f, TSTRING, "SRC_NAME", jname, NULL, &status);
    fits_update_key(f, TSTRING, "RA", ra, NULL, &status);
    fits_update_key(f, TSTRING, "DEC", dec, NULL, &status);
    fits_update_key(f, TSTRING, "STT_CRD1", ra, NULL, &status);
    fits_update_key(f, TSTRING, "STP_CRD1", ra, NULL, &status);
    fits_update_key(f, TSTRING, "STT_CRD2", dec, NULL, &status);
    fits_update_key(f, TSTRING, "STP_CRD2", dec, NULL, &status);
        
    /* write it out */
    fits_close_file(f, &status);

    exit(0);
}
