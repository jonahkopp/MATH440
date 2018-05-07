Running the program:

First, change the message.txt file to be any message as long as it is under 2000 characters in length

On Sayer's machine:

1. Navigate to final_proj directory containing all project files such as main.f90 and mod_file.f90

2. Compile the program by typing "make main_exe" without the quotes

3. Execute by typing "mpiexec -machinefile lab_ma_file -np P ./main_exe" without quotes, where P is the number of processing cores

4. The message has now been encoded, decoded, and written to file output.txt. Open by typing "emacs output.txt" or replace emacs with any text editor



On MIO:

1. Get onto MIO and navigate to the same directory as above steps

2. Type "module purge", hit Enter, type "module load StdEnv", hit Enter, type "module load impi/intel/latest", hit Enter (all without quotes) to load only necessary modules

3. Type "mpiifort -c mod_file.f90" without quotes to compile the module file so that step 5 will work properly

4. Type "emacs final_proj.slurm" without quotes to open the batch file, then change number of MIO nodes, number of tasks per node, and number of total tasks
   to whatever you'd like, so long as number of tasks per node is 16 or less and number of total tasks = number of nodes * number of tasks per node 

5. Type "mpif90 -fopenmp mod_file.f90 main.f90 -o main_exe" without quotes to compile the main program and create main_exe

6. Type "sbatch final_proj.slurm" without quotes to send the job to MIO

7. Once job is finished, open con.out using a text editor to view the runtime and open the output.txt file to view the original message which has been encoded and decoded
