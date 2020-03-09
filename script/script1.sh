#$ -S /bin/bash
#$ -N thimble_Nt_2
#$ -l h_rt=10:00:00
#$ -l h_vmem=800m
#$ -l h_fsize=500m
#$ -q run64new

export OMP_NUM_THREADS=8

i=2

cp /lpt/jquarroz/python/U1_python/param.py .
cp /lpt/jquarroz/python/U1_python/$2.py .

mkdir ./data

python3 $2.py 1 $1

cp ./data/accept_p_Nt_${i}_B_${1}.txt /lpt/jquarroz/python/U1_python/data/
cp ./data/accept_q_Nt_${i}_B_${1}.txt /lpt/jquarroz/python/U1_python/data/
cp ./data/expect_Nt_${i}_B_${1}.txt /lpt/jquarroz/python/U1_python/data/
cp ./data/x4_Nt_${i}_B_${1}.txt /lpt/jquarroz/python/U1_python/data/
cp ./data/y4_Nt_${i}_B_${1}.txt /lpt/jquarroz/python/U1_python/data/

cp ./param.py /lpt/jquarroz/python/U1_python/data/
