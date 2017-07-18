#!/bin/bash

exit_module(){
  printf "\nExiting ["
  for i in {1..70}; do sleep 0.003; printf "#"; done
  printf "]\n"; exit;
}


kill_subshells(){
  >&2 printf " ----------------------------------- \n"
  >&2 printf '%s\t' "${tool_names[@]}"
  >&2 printf '%s\t' "${pids[@]}"
  >&2 printf "\n"
  for pid in ${pids[@]};
  do
    >&2 printf "killing pid: $pid\n"
    #kill all subprocesses of tools
    ps -o pid,ppid | grep $pid | tail -n+2 | awk -F  " " '{print $1}' | xargs kill -15
    kill -15 $pid
  done
  >&2 printf " Killing done\n"
  >&2 printf " ----------------------------------- \n"
  exit_module
}


cleanup(){
  >&2 printf "cleaning up: \n"
  exit_status=$?
  >&2 printf "exit status: $?\n"
  if [[ $exit_status -ne 0 ]]; then
    #save stdout/stderr
    mkdir $parsnp_path/$NOW/logs
    mkdir $kSNP_path/$NOW/logs

    mv $parsnp_output/stderr  $parsnp_path/$NOW/logs
    mv $parsnp_output/stdout  $parsnp_path/$NOW/logs

    mv $ksnp_output/stderr $kSNP_path/$NOW/logs
    mv $ksnp_output/stdout $kSNP_path/$NOW/logs
  fi

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
