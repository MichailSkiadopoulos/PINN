JOB_NUMBER = [6, ]

for i in JOB_NUMBER:
    
    nodes = 1
    nCPUs = int(3)
    memPernode = int(88)
    precision = "double" 
    
    fh=open('Job-%d.pbs' %(i),'w')
    
    fh.write("#!/bin/bash\n")
    fh.write("#SBATCH -N %d\n" %(nodes))
    fh.write("#SBATCH --tasks-per-node=1\n")
    fh.write("#SBATCH --cpus-per-task=%d\n" %(nCPUs))
    fh.write("#SBATCH --mem-per-cpu=%dGB\n" %(memPernode))
    fh.write("#SBATCH --constraint='sc&icelake'\n") #Its standard core & icelake
    fh.write("#SBATCH --account=open\n")
    fh.write("#SBATCH --partition=open\n")
    fh.write("#SBATCH --time=48:00:00\n")
    fh.write("#SBATCH -o STDOUT.out\n")
    fh.write("#SBATCH --mail-type=end   # send email when job ends\n")
    fh.write("#SBATCH --mail-type=fail  # send email if job fails\n")
    fh.write("#SBATCH --mail-user=mbs6742@psu.edu\n")
    fh.write("\n")
    fh.write('echo "Job started on `hostname` at `date`"\n')
    fh.write("\n")
    fh.write("# Load Module\n")
    fh.write("module purge\n")
    fh.write("module load abaqus\n")    
    fh.write("\n")     
    fh.write("# Directory from which the job is submitted\n")
    fh.write("cd /storage/home/mbs6742/work/Wrong_sims3/\n")
    fh.write("\n")    
    fh.write("# Run job\n")
    fh.write("abaqus job=Job-%d  {} interactive ask_delete=OFF\n".format(precision) %(i))
    fh.write("\n")
    fh.write('echo "Job Ended at `date`"\n')
    
    fh.close()
