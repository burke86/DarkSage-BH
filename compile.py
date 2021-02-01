"""
Run this python wrapper to update hard-coded variables in the Dark Sage codebase, pointing to a relevant parameter file, then `make'.
"""

from __future__ import print_function
import subprocess
import os
import multiprocessing

# Initialisation
hard_variable_names = ["N_BINS", "N_AGE_BINS"]
wd = "code/"
og_file = wd + "core_allvars.h"
old_file = og_file
temp_count = 0

# Rename raw_input() / input() to accomodate both Python 2 and Python 3
try:
    cmnd_input = raw_input
except NameError:
    cmnd_input = input

# Provide input for parameter file name.  If nothing provided, will skip updating the hard parameters and just compile as is.
par_file = cmnd_input("Provide path to a parameter file to update hard parameters (or hit enter to skip)\n")

if len(par_file)>0:

    # Read parameter file to get relevant parameter values, checking for common reasons  why the fed name might not work  
    successful = False
    attempt_count = 0
    while not successful:
        if attempt_count==1 and par_file[-4:]!='.par':
            w = par_file.find('.')
            if w==-1:
                par_file += '.par' # If user is too lazy to add '.par' to the end
            else:
                par_file = par_file[:w] + '.par' # If user mistyped '.par' for something else (e.g. '.py')
            
        if attempt_count==2 and par_file[:6]!='input/': 
            par_file = 'input/' + par_file # If user forgot to specify the input directory
        
        if attempt_count==3:
            print("Could not locate the parameter file")
            quit()
        
        try:
            f = open(par_file, 'r')
            successful = True
        except IOError:
            attempt_count += 1
            successful = False

        
    hard_variable_values = []
    lines = f.readlines()
    for name in hard_variable_names:
        lstr = len(name)
        for line in lines:
            words = line.split()
            try: 
                first_word = words[0]
            except IndexError:
                continue
            if first_word==name:
                hard_variable_values += [words[1]]
                break
    
    # Create temporary files, updating the next hard parameter with each iteration
    for var, name in zip(hard_variable_values, hard_variable_names):
        temp_count += 1
        temp_file = wd + "core_allvars_temp_"+str(temp_count)+".h"
        os.system("sed \'s/^#define "+name+".*/#define "+name+" "+var+"/\' "+old_file+ " >"+temp_file)
        old_file = temp_file
        
    subprocess.call(['mv', temp_file, og_file]) # replace the original file with the latest iteration

    # remove temporary files
    for i in range(1,temp_count):
        temp_file = wd + "core_allvars_temp_"+str(i)+".h"
        subprocess.call(['rm', temp_file])

# clean and compile Dark Sage
subprocess.call(['make', 'clean'])
subprocess.call(['make'])

if len(par_file)>0:
    print("\nIf compilation successful, run with the following command:")
    print("./darksage", par_file)
    print("\nOr run in parallel with:")
    print("mpirun -np ", multiprocessing.cpu_count(), "./darksage", par_file, "\n")
