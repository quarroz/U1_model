#$ -S /bin/bash
#$ -N thimble_error
#$ -l h_rt=10:00:00
#$ -l h_vmem=800m
#$ -l h_fsize=500m
#$ -q run64new

export OMP_NUM_THREADS=8

Nt=20
beta=10

touch data_${1}

cp /lpt/jquarroz/python/U1/param.py .
cp /lpt/jquarroz/python/U1/U1.py .

mkdir ./data

python3 U1.py 2 10 > data_${1}.txt

#cp ./data/expect_Nt_${Nt}_B_${beta}.txt /lpt/jquarroz/python/U1/data/expectation_${1}_False.txt
cp ./data_${1} /lpt/jquarroz/python/U1/data/
#cp ./param.py /lpt/jquarroz/python/U1/data/
