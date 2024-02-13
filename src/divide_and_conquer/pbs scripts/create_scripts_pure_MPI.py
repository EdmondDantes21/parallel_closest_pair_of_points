def write_to_file(nodes, cores, n):
    filename = "pure_mpi_dac_" + str(nodes) + "_" + str(cores) + "_" + str(n) + ".sh"
    f = open(filename, "a")
    f.write("#!/bin/bash\n\n")
    f.write("#PBS -l select=" + str(nodes) + ":ncpus=" + str(cores) + ":mem=4gb\n\n")
    f.write("#set max execution time\n")
    f.write("#PBS -l walltime=0:15:00\n\n")
    f.write("#PBS -q short_cpuQ\n\n")
    f.write("mpiexec -n " + str(cores * nodes) + " ./pure_mpi_closest_pair " + str(2**n))


nodes = [1,1,1,1,2,4]
cores = [2,4,8,16,16,16]
N = [20,21,22,23,24,25]

for i in range(6):
    for j in range(6):
        write_to_file(nodes[i], cores[i], N[j])
        