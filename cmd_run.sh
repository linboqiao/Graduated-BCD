make clean;make

## run PFBCD on shared memory distributed system
yhrun -J PFBCD-E -p th_ft1 -N 1 -n 1 --wait=0 --mincpus=32 --mem=16000 -v ./main
