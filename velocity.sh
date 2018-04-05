#!/bin/bash 
#######
#######
#SBATCH -N 1
#SBATCH -J dcn
#SBATCH -t 48:00:00
#SBATCH -p share
#SBATCH -o dcn.out
#SBATCH -e dcn.err
#######
#

ulimit -s unlimited
#
. /etc/cmod/sh.init
#
./dcn.x
#
