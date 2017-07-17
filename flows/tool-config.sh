#!/bin/bash

compile_parsnp(){
  cd $parsnp_path/bin
  # Does parsnp executeables exists && is  right OS exec -> ASSUMING other executables work
  compiled=0
  # compiles if executable is wrong format || doesn't exist.
  if [ "$OS" == "Darwin" ];
  then
    parsnpcompiled=$(file parsnp | grep -o "Mach-O 64-bit executable x86_64")
      if [[ ! $parsnpcompiled == "Mach-O 64-bit executable x86_64" ]]; then
        >&2 printf "parsnp NOT compiled or not compiled CORRECTLY!\n compiling ....."
        bash build_parsnp_nix.sh 2&>1 build_parsnp.std
        if [[ $? -ne 0 ]]; then
          >&2 printf "Compilation FAILED, exit code: $?"
          exit 1
        fi
        compiled=1
      fi
  elif [ "$OS" == "Linux" ];
  then
    parsnpcompiled=$(file parsnp | grep -o "ELF 64-bit LSB executable")
      if [[ ! $parsnpcompiled == "ELF 64-bit LSB executable" ]]; then
        >&2 printf "parsnp NOT compiled or not compiled CORRECTLY!\n compiling $parsnpcompiled....."
        bash build_parsnp_nix.sh 2&>1 build_parsnp.std
        if [[ $? -ne 0 ]]; then
          >&2 printf "Compilation FAILED, exit code: $?"
          exit 1
        fi
        compiled=1
      fi
  fi
  cd ../
}


# tool runs definitions:

run_ksnp(){
  cd $kSNP_path
  #if previous run failed, multiple entries in these files will give problems, (tsch doesnt support functions, or trap clean_function EXIT should be used)
  ./kSNP3 -in $input_files -outdir $ksnp_output -k 13 -CPU $CPUS \
  -kchooser "1" -ML -path $kSNP_path  > $ksnp_output/stdout  2> $ksnp_output/stderr
}

run_parsnp(){
    compile_parsnp
    printf "Parsnp path: $parsnp_path Compiled: $compiled (0=no,1=yes)\n" > $parsnp_output/stdout
    python ./Parsnp.py $parsnp_path -r $ref -d $genome_path -p $CPUS -o $parsnp_output \
    > $parsnp_output/stdout  2> $parsnp_output/stderr
}


#fetching all tool definitions
#typeset -F gives function names only, -f gives whole function of functions
all_tools_functions="typeset -F | grep -o '\<run_.*\>'"
declare -a tool_array
while read -r line; do
   tool_array+=("$line")
done <<< "`eval $all_tools_functions`"  # inputs file to be read line by line
