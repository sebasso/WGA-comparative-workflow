#!/bin/bash

exit_module(){
  printf "\nExiting ["
  for i in {1..70}; do sleep 0.003; printf "#"; done
  printf "]\n"; exit;
}


kill_subshells(){
  printf " ----------------------------------- \n"
  printf " Killing subshells\n\n"
  printf '%s\t' "${tool_names[@]}"
  printf '%s\t' "${pids[@]}"
  printf "\n"
  for pid in ${pids[@]};
  do
    printf "killing pid: $pid\n"
    #kill all subprocesses of tools
    ps -o pid,ppid | grep $pid | tail -n+2 | awk -F  " " '{print $1}' | xargs kill -15
    kill -15 $pid
  done
  printf " \n\nKilling done\n"
  printf " ----------------------------------- \n"
  exit_module
}


cleanup(){
  printf "cleaning up: \n"
  printf "exit status: $?\n"
  if [ ! -z $run_specific ]; then
    rm -rf $run_specific
  fi

  if [ -d $parsnp_output ]; then
    rm -rf $parsnp_output
  fi
  #ksnp
  if [ -d $ksnp_output ]; then
    rm -rf $ksnp_output
  fi
  if [ -f $kSNP_path/fasta_list ]; then
   rm $kSNP_path/fasta_list
  fi
  if [ -f $kSNP_path/annotate_genomes ]; then
   rm $kSNP_path/annotate_genomes
  fi
  if [ -f $kSNP_path/in_list ]; then
   rm $kSNP_path/in_list
  fi

  kill_subshells
}
