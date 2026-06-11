# Edit Makefile.am to specify the route for psrfits template file
# To compile, run:
./bootstrap

./configure --prefix=install_root

make & make install
