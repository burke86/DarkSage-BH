#
#    ICRAR - International Centre for Radio Astronomy Research
#    (c) UWA - The University of Western Australia, 2019
#    Copyright by UWA (in the framework of the ICRAR)
#    All rights reserved
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#    MA 02111-1307  USA
#

import logging
import multiprocessing
import os
import shutil
import subprocess
import tempfile
import time

import numpy as np

import common


logger = logging.getLogger(__name__)

def count_jobs(job_name):
    """Returns how many jobs with self.jobs_name are currently queued or running"""

    try:
        out, err, code = common.exec_command("squeue")
    except OSError:
        raise RuntimeError("Couldn't run squeue, is it installed?")

    if code:
        raise RuntimeError("squeue failed with code %d: stdout: %s, stderr: %s" % (code, out, err))

    lines_with_jobname = [l for l in out.splitlines() if job_name in l]
    return len(lines_with_jobname)

def _exec_shark(msg, cmdline):
    logger.info('%s with command line: %s', msg, subprocess.list2cmdline(cmdline))
    out, err, code = common.exec_command(cmdline)
    if code != 0:
        logger.error('Error while executing %s (exit code %d):\n' +
                     'stdout:\n%s\nstderr:\n%s', cmdline[0], code,
                     common.b2s(out), common.b2s(err))
        raise RuntimeError('%s error' % cmdline[0])


def _to_shark_options(particle, space):
    """Given `particle` in `space` return an iterable with the corresponding
    shark options settings corresponding to that particle"""
    for value, name, is_log in zip(particle, space['name'], space['is_log']):
        if is_log:
            value = 10 ** value
        yield '%s=%s' % (name, value)

def _evaluate(constraint, stat_test, modeldir, subvols):
    y_obs, y_mod, err = constraint.get_data(modeldir, subvols)
    return stat_test(y_obs, y_mod, err)

count = 0
def run_shark_hpc(particles, *args):
    """
    - Handler function for running PSO on Shark on a SLURM based cluster.
    - Swarm size and number of iterations need to be set within the script for now
    - Function needs the relative path to a Shark config file under the -c option
    - For now the subprocess call within must be altered if you are changing shark submit options
    - To find appropriate memory allocations peruse the initial output lines of each particle
    """

    global count

    opts, space, subvols, statTest = args

    # Prepare the file that will be used by the shark submission scripts
    # to determine which values shark will be run for. We put a final \n so the
    # final line gets properly counted by wc (used by shark-submit)
    shark_options = [
        ' '.join(['-o "%s"' % option for option in _to_shark_options(particle, space)])
        for particle in particles
    ]
    positions_fname = tempfile.mktemp('particle_positions.txt')
    logger.info('Creating particle positions file at %s', positions_fname)
    with open(positions_fname, 'wt') as f:
        f.write('\n'.join(shark_options) + '\n')

    # Submit the execution of multiple shark instances, one for each particle
    job_name = 'PSOSMF_%d' % count
    shark_output_base = os.path.join(opts.outdir, job_name)
    cmdline = ['./shark-submit', '-S', opts.shark_binary, '-w', opts.walltime,
               '-n', job_name, '-O', shark_output_base, '-E', positions_fname,
               '-V', ' '.join(map(str, subvols))]
    if opts.account:
        cmdline += ['-a', opts.account]
    if opts.queue:
        cmdline += ['-Q', opts.queue]
    if opts.nodes:
        cmdline += ['-N', str(opts.nodes)]
    else:
        cmdline += ['-m', opts.memory, '-c', str(opts.cpus)]
    cmdline.append(opts.config)
    _exec_shark('Queueing PSO particles', cmdline)

    # Actually wait for the jobs to finish...
    while count_jobs(job_name) > 0:
        time.sleep(10)

    ss = len(particles)
    fx = np.zeros([ss, 3])
    for i in range(ss):
        _, simu, model, _ = common.read_configuration(opts.config)
        particle_outdir = os.path.join(shark_output_base, str(i))
        modeldir = common.get_shark_output_dir(particle_outdir, simu, model)
        fx[i] = [_evaluate(c, statTest, modeldir, subvols) for c in opts.constraints]
        if not opts.keep:
            shutil.rmtree(particle_outdir)

    fx = np.sum(fx, 1)
    logger.info('Particles %r evaluated to %r', particles, fx)

    # this global count just tracks the number of iterations so they can be saved to different files
    count += 1

    return fx

# this is the 'func' that is passed into the pso routine in pso.py
# this doesn't require the new pso.py file, can just run on the original pyswarm function
# only the hpc module above needs the modified pso function, as it has parallelization updates
# the vanilla version only does one particle at a time, whereas the hpc version takes multiple particles at once
def run_shark(particle, *args):

# space is the thing containing the parameter values for the model
    opts, space, subvols, statTest = args

    pid = multiprocessing.current_process().pid
    shark_output_base = os.path.join(opts.outdir, 'output_%d' % pid)
    _, simu, model, _ = common.read_configuration(opts.config)
    modeldir = common.get_shark_output_dir(shark_output_base, simu, model)

# here is where things are executed.  I will need to modify my parameter file here and then run Dark Sage
    cmdline = [opts.shark_binary, opts.config,
               '-o', 'execution.output_directory=%s' % shark_output_base,
               '-o', 'execution.simulation_batches=%s' % ' '.join(map(str, subvols))]
    for option in _to_shark_options(particle, space):
        cmdline += ['-o', option]
    _exec_shark('Executing shark instance', cmdline)

    total = sum(_evaluate(c, statTest, modeldir, subvols) for c in opts.constraints)
    logger.info('Particle %r evaluated to %f', particle, total)

    if not opts.keep:
        shutil.rmtree(shark_output_base)

    return total



def run_darksage(particle, *args):
    
    # space is the thing containing the parameter values for the model
    opts, space, subvols, statTest = args
    
    # create/clear directory for temporary Dark Sage output
    spid = str(multiprocessing.current_process().pid)
#    modeldir = '/Users/adam/DarkSage/autocalibration/DS_output_'+spid+'/'
    modeldir = '/fred/oz245/DarkSage_output/MTNG/autocalibration/02/DS_output_' + spid + '/'
    if not os.path.exists(modeldir): os.makedirs(modeldir)
    if os.path.isfile(modeldir+'model_z0.000_0'): subprocess.call(['rm', modeldir+'model*'])

    # copy template parameter file and edit accordingly
    slash = len(opts.config) - opts.config[::-1].find('/') - 1
#    temp_filename = '/Users/adam/DarkSage/autocalibration/' + opts.config[slash+1:-4] + '_' + spid + '_temp.par'
    temp_filename = '/fred/oz245/DarkSage_output/MTNG/autocalibration/02/' + opts.config[slash+1:-4] + '_' + spid + '_temp.par'
    #
    f = open(opts.config).readlines()
    s = open(temp_filename, 'w')
    #
    Np = len(space['name'])
    Ndone = 0
    for l, line in enumerate(f):
        if line[:9] == 'OutputDir': f[l] = 'OutputDir              '+modeldir+'\n'
        for p in range(Np):
            if line[:len(space['name'][p])] == space['name'][p]: 
                f[l] = space['name'][p]+'          '+str(round(particle[p],5))+'\n'
                Ndone += 1
        if Ndone==Np: break
    #
    s.writelines(f)
    s.close()

    # currently assuming that serial in parallel here is the best
    print(temp_filename)
#    cmdline = ['mpirun', '-np', '8', opts.shark_binary, temp_filename]
    cmdline = [opts.shark_binary, temp_filename]
    _exec_shark('Executing Dark Sage instance', cmdline)

    # changed to now take the weighted sum of logs of the (summed) reduced chi^2 of each constraint
    total = 10**sum(np.log10(np.sum(_evaluate(c, statTest, modeldir, subvols))*c.weight) for c in opts.constraints)
    logger.info('Particle %r evaluated to %f', particle, total)

    shutil.rmtree(modeldir)
    return total
