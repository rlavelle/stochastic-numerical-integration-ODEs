mkdir protein-trials-noclb1-noclb3 protein-trials-noclb3 protein-trials-noclb1-noclb3-noclb4 protein-trials-noclb1 protein-trials-noclb1-noclb4 protein-trials-noclb4 protein-trials-noclb3-noclb4 protein-trials-wt

python3 gillespie.py 15 /home/cvanoeve/stochastic-numerical-integration-ODEs/protein-trials-noclb1-noclb3 0 0 0 0 0.2 0.1 &
python3 gillespie.py 15 /home/cvanoeve/stochastic-numerical-integration-ODEs/protein-trials-noclb3 0.002 0.2 0 0 0.2 0.1 &
python3 gillespie.py 15 /home/cvanoeve/stochastic-numerical-integration-ODEs/protein-trials-noclb1-noclb3-noclb4 0 0 0 0 0 0 &
python3 gillespie.py 15 /home/cvanoeve/stochastic-numerical-integration-ODEs/protein-trials-noclb1 0 0 0.002 0.5 0.2 0.1 &
python3 gillespie.py 15 /home/cvanoeve/stochastic-numerical-integration-ODEs/protein-trials-noclb1-noclb4 0 0 0.002 0.5 0 0 &
python3 gillespie.py 15 /home/cvanoeve/stochastic-numerical-integration-ODEs/protein-trials-noclb4 0.002 0.2 0.002 0.5 0 0 &
python3 gillespie.py 15 /home/cvanoeve/stochastic-numerical-integration-ODEs/protein-trials-noclb3-noclb4 0.002 0.2 0 0 0 0 &
python3 gillespie.py 15 /home/cvanoeve/stochastic-numerical-integration-ODEs/protein-trials-wt 0.002 0.2 0.002 0.5 0.2 0.1 &
