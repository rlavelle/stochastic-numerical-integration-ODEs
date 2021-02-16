# stochastic-numerical-integration-ODEs

- Rowan Lavelle <rowlavel@iu.edu>
- Curtis Goss <cvanoeve@iu.edu>

This project is not meant to be ran on a regular work computer, however it is runnable. It will take up 5 cores, and output multiple gigabytes of files. The work computer friendly version is the notebooks, which are used for testing before creating the scripts which run the stochastic gillespie process in parallel.

The goal of this project is to create a quantitative understanding of the meiotic termination network by looking at the counts of certain proteins in a cell while it is going through its phases. The stochastic simulation of this process is coupled differential equations to control birth and death rates of different proteins using the gillespie process.

We compare a stochastic model to a deterministic model to create more real life scenarios. The second half of this project is done in the bio lab, where proteins are florescently marked so their sizes can be recorded. This data in turn is what helps create the constants and values in the coupled differential equation system, which inturn controls the stochastic process. 
