#!/bin/csh
#SBATCH --time=00:10:00  #Walltime in hours:minutes:seconds
#SBATCH --ntasks=1  #Number of Processors needed 
#SBATCH --nodes=1   #Number of cluster nodes needed
#SBATCH -o FDTD.%j.out  #Name of stdout, using the job number (%j) 
#SBATCH -e FDTD.%j.err  #Name of stdout error file using the job number
#SBATCH --mail-type=FAIL,BEGIN,END  #Flags on when to send notification emails
#SBATCH --mail-user=<youruniversityid>@utah.edu  #E-mail to send notifications to

#SBATCH --account=ECE5340 #charge the account set up for this class
#SBATCH --partition=kingspeak #you can use kingspeak or lonepeak

# Set the WORKING directory (this is the directory from which you submitted this file to 
# the scheduler)
setenv WORKDIR $SLURM_SUBMIT_DIR
echo "Set Working Directory"

# Set the SCRACTCH directory (this is the directory where your code will run)
setenv SCRDIR /scratch/kingspeak/serial/$USER/$SLURM_JOBID

# Create the SCRATCH directory 
mkdir -p $SCRDIR

# Copy your executable from the WORKING directory to the SCRATCH directory
cp $WORKDIR/a.out $SCRDIR

# Move over to the SCRATCH directory
cd $SCRDIR
echo "Set Scratch Directory"

# Set the EscXE
setenv EXE $SCRDIR/a.out
echo "Set EXE"
echo "Set number of processors"

# Load the appropriate modules
module purge
module load intel
module load impi
echo "Loaded Modules"

setenv I_MPI_FABRICS shm:ofa

echo "Start of calculations: `date`"
# Run the code
mpirun -bootstrap slurm -np $SLURM_NTASKS $EXE  >& out.$SLURM_JOBID.out
echo "End of Calculations: `date`"

# Copy any output from the SCRATCH directory to the WORKING directory for easy access
cp -r $SCRDIR/*.dat $WORKDIR
cp -r $SCRDIR/*.out $WORKDIR

#rm -Rf $SCRDIR
