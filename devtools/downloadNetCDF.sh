#!/bin/bash

# install into here:
PREFIX=/opt/


HDF5VERSION=hdf5-1.8.17
NETCDF4VERSION=netcdf-4.4.1

echo $0: Downloads $NETCDF4VERSIO and $HDF5VERSION \(to /tmp/\) and configures them to be installed in $PREFIX
echo For updated version change this script, but do not use a > hdf5 1.10.x  version.
echo install dirs are: $PREFIX/$NETCDF4VERSION-parallel and $PREFIX/$HDF5VERSION-parallel

read -p "To continue press [Enter], [n|ctrl-c] to abort: "

if [ "$REPLY" != "" ]; then
   exit
fi

cd /tmp/

if [ ! -e $HDF5VERSION.tar.gz ]; then
    wget  https://support.hdfgroup.org/ftp/HDF5/current/src/$HDF5VERSION.tar.gz
fi

if [ ! -e $NETCDF4VERSION.tar.gz ]; then
  wget ftp://ftp.unidata.ucar.edu/pub/netcdf/$NETCDF4VERSION.tar.gz
fi

tar xfz $HDF5VERSION.tar.gz
tar xfz $NETCDF4VERSION.tar.gz

cd $HDF5VERSION

CC=mpicc ./configure  --prefix=/$PREFIX/$HDF5VERSION-parallel/ --enable-parallel

make -j3  
#make check 

echo
echo

echo $0:  Next step needs sudo pw, which will probably automatically also be used for the next. If in doubt hit Ctrl-C and abort, do make install manually, then look in this script and ./configure and make netcdf afterwards with subsequent make install.

echo echo$0: 'sudo make install'
sudo make install

cd ../$NETCDF4VERSION

CC=mpicc CPPFLAGS=-I/$PREFIX/$HDF5VERSION-parallel/include/ LDFLAGS=-L/$PREFIX/$HDF5VERSION-parallel/lib/ ./configure --enable-parallel-tests --prefix=/$PREFIX/$NETCDF4VERSION-parallel

make  -j3
#make check

echo $0: calling make install on netcdf
sudo make install


cd ..

echo
echo
echo $0: Do not forget to add export LD_LIBRARY_PATH=/$PREFIX/$NETCDF4VERSION-parallel/lib:\$LD_LIBRARY_PATH
echo $0: you need to cleanup downloaded files and unpacked dirs files in /tmp yourself.

