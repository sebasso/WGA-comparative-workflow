#!/bin/bash
#assuming brew or similar installed gcc to PATH /usr/local/bin/
# && dependent on automake and autoconf  install with brew install automake
# NOTE: gcc,cc and cxx must be exported so this script find their location.
#export CC=/usr/local/Cellar/gcc@4.9/4.9.4/bin/gcc-4.9
#export CXX=/usr/local/Cellar/gcc@4.9/4.9.4/bin/c++-4.9

OS=`uname`

if [ "$OS" == "Darwin" ];
then
printf "\nOS: $OS\n"
#gccpath=`which gcc`
printf "\n gccpath: \n$gcc\n"
exec=`echo $gcc | grep -o /usr/local/bin/gcc-'[0-9]'."[0-9]"`
printf "exec:\n $exec\n"
  if [ ! -z "$exec" ];
  then
    printf "valid gcc compiler found \n installing parsnp with $gcc .."

    cd muscle/
    ./autogen.sh
    ./configure --prefix=`pwd` CXXFLAGS='-fopenmp'
    make install
    cd ..
    ./autogen.sh
    ./configure --prefix=`pwd` CXXFLAGS='-fopenmp'
    make install

  elif [ "$gccpath" == *"/usr/bin/gcc"* ];
  then
    printf "OS is using clang compiler, need to use gcc for open \n openmp flag compilation"
    printf "\n either by installing homebrew with: \n /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)" \n and brew search gcc || brew install gcc-yourversion \n Manually\n https://gcc.gnu.org/wiki/InstallingGCC"
  else
    printf "unknown compiler \n aborting ...\n"
  fi
elif [ "$OS" == "Linux" ];
then
    cd muscle/
    ./autogen.sh
    ./configure --prefix=`pwd`
    make install
    cd ..
    ./autogen.sh
    ./configure
    make install
else
  printf "unknown OS try running ./autogen.sh in muscle dir and parsnp dir to check for compatibility"
fi
