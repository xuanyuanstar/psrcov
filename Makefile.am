AM_CFLAGS = -I@top_srcdir@/ @CFITSIO_CFLAGS@ @FFTW_CFLAGS@

bin_PROGRAMS= vdif2psrfitsALMA vdif2psrfitsPico UDP2psrfits set_coor UDP2dada19BEAM UDP2dadaUWB nuppi2dada vdif2dadaALMA vdif2dadaEB
lib_LTLIBRARIES=libVDIF.la

libVDIF_la_SOURCES = dec2hms.c downsample.c polyco.c vdifio.c write_psrfits.c cvrt2to8.c mjd2date.c getVDIFFrameDetection.c getUDPDetection.c date2mjd.c date2mjd_ld.c ascii_header.c
libVDIF_la_LIBADD = @CFITSIO_LIBS@ @FFTW_LIBS@ 

vdif2psrfitsPico_SOURCES = vdif2psrfitsPico.c
vdif2psrfitsPico_LDADD = libVDIF.la @CFITSIO_LIBS@ @FFTW_LIBS@ -lfftw3f_threads 

vdif2psrfitsALMA_SOURCES = vdif2psrfitsALMA.c
vdif2psrfitsALMA_LDADD = libVDIF.la @CFITSIO_LIBS@ @FFTW_LIBS@ -lfftw3f_threads

vdif2dadaALMA_SOURCES = vdif2dadaALMA.c
vdif2dadaALMA_LDADD = libVDIF.la

vdif2dadaEB_SOURCES = vdif2dadaEB.c
vdif2dadaEB_LDADD = libVDIF.la

UDP2psrfits_SOURCES = UDP2psrfits.c
UDP2psrfits_LDADD = libVDIF.la @CFITSIO_LIBS@ @FFTW_LIBS@

set_coor_SOURCES = set_coor.c
set_coor_LDADD = @CFITSIO_LIBS@

UDP2dada19BEAM_SOURCES = UDP2dada19BEAM.c

UDP2dadaUWB_SOURCES = UDP2dadaUWB.c

nuppi2dada_SOURCES = nuppi2dada.c
nuppi2dada_LDADD = libVDIF.la

AM_CPPFLAGS = -DPSRFITS_TEMPLATE_DIR='"/cluster/pulsar/kliu/Soft/psrcov"'

ACLOCAL_AMFLAGS = -I config
