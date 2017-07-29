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
}


cleanup(){
  >&2 printf "cleaning up: \n"
  exit_status=$?
  >&2 printf "exit status: $?\n"
  #remove subprocesses from directories first or rm can fail if subprocesses is in an writing phase
  kill_subshells
  #if [[ $exit_status -ne 0 ]]; then
    #save stdout/stderr

  crashed = logs/$NOW$name
  mkdir -p $parsnp_path/logs/$NOW
  mkdir -p $kSNP_path/logs/$NOW
  # save logs
  mv $parsnp_output/stderr $parsnp_path/$crashed
  mv $parsnp_output/stdout $parsnp_path/$crashed
  mv $ksnp_output/stderr $kSNP_path/$crashed
  mv $ksnp_output/stdout $kSNP_path/$crashed
  #Save files
  mv $ksnp_output/ksnp.tree $kSNP_path/$crashed
  mv $ksnp_output/kSNP_SNPs_POS_formatted.tsv $kSNP_path/$crashed
  mv $ksnp_output/SNPs_all $kSNP_path/$crashed
  mv $ksnp_output/SNPs_all_matrix.fasta $kSNP_path/$crashed
  mv $parsnp_output/parsnp.tree $parsnp_path/$crashed
  mv $parsnp_output/parsnp_SNPs_POS_formatted.tsv $parsnp_path/$crashed
  mv $parsnp_output/parsnp.ggr $parsnp_path/$crashed

  if [ ! -z $run_specific ]; then
    rm -rf $run_specific
  fi

  cleanup_junk
}

cleanup_junk(){
  if [ -d $parsnp_output ]; then
    rm -rf $parsnp_output
  fi

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
  exit_module
}
