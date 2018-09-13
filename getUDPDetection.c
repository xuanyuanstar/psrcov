#include <malloc.h>
#include <complex.h>
#include <fftw3.h>
void getDetection(float p0r, float p0i, float p1r, float p1i, float *det, char dstat);

void getUDPDetection(const char *src_p0, const char *src_p1, int bbytes, float det[][4], char dstat)
{
  fftwf_complex *out_p0,*out_p1;
  fftwf_plan pl0,pl1;
  float dets[4],*in_p0,*in_p1;
  int i,j,k,Nts;
  unsigned char *buff8[2];

  //Memo for FFT for two pols                                              
  out_p0 = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*bbytes/2+1);
  out_p1 = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*bbytes/2+1);
  in_p0 = (float *) fftwf_malloc(sizeof(float)*bbytes);
  in_p1 = (float *) fftwf_malloc(sizeof(float)*bbytes);

  // Initialize
  for(i=0;i<bbytes;i++)
    {
      in_p0[i]=(float)((int)src_p0[i]);
      in_p1[i]=(float)((int)src_p1[i]);
    }

  // Perform FFT
  pl0 = fftwf_plan_dft_r2c_1d(bbytes, in_p0, out_p0, FFTW_ESTIMATE);
  pl1 = fftwf_plan_dft_r2c_1d(bbytes, in_p1, out_p1, FFTW_ESTIMATE);
  fftwf_execute(pl0);
  fftwf_execute(pl1);

  // Make detection
  for(i=0;i<bbytes/2+1;i++)
    {
      // Detection each FFT spetrum
      getDetection(creal(out_p0[i]),cimag(out_p0[i]),creal(out_p1[i]),cimag(out_p1[i]),dets,dstat);
      for(j=0;j<4;j++)
	det[i][j]=dets[j];
    }

  //Free up memo
  fftwf_free(in_p0);
  fftwf_free(in_p1);
  fftwf_free(out_p0);
  fftwf_free(out_p1);
  fftwf_destroy_plan(pl0);
  fftwf_destroy_plan(pl1);
}
