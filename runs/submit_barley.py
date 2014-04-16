#!/usr/bin/python

import sys
import os
import subprocess as subp
import shlex
import time

qname = "main.q"

def ensure_dir_exists(dir):
    if (not os.path.exists(dir)):
        os.makedirs(dir)

def runcmd(cmd, input=None):
    cmd = shlex.split(cmd)
    #print "cmd is:", cmd
    #print "input is:", input
    p = subp.Popen(cmd, stdin=subp.PIPE, stdout=subp.PIPE, stderr=subp.PIPE)
    # uncomment next line to see exactly what is piped to Popen
    # p.stdin = sys.stderr

    (stdoutdata, stderrdata) = p.communicate(input=input)
    #print "stdoutdata:", stdoutdata[:-1]
    #print "stderrdata:", stderrdata[:-1]
    return(stdoutdata[:-1], stderrdata[:-1])

def qsub_launch(dir, prefix, line, qname=qname, jobname = dir, qsub_extras=""):
    # run with qsub

    if (1): # for Sun Grid Engine (SGE) queue system
        cmd = "qsub"
        if (qname != ""):
            #cmd = cmd + " -w w -q %s -l qname=%s" % (qname, qname.split("@")[0])
            cmd = cmd + " -w w -q %s" % (qname)
        if (qsub_extras != ""):
            cmd = cmd + " "+qsub_extras
        print "qsub cmd is:", cmd

        base = os.getenv("PWD")+"/"+dir
        input=""
        input=input+"#!/bin/bash\n"
        input=input+"cd %s\n" % (base) # follows symlinks properly with "PWD"
        input=input+"#$ -N %s\n" % (jobname)
        if (1): # change to a 0 if you want no .out file
            input=input+"#$ -o %s/outs/%s.out\n" % (base, prefix)
        else:
            input=input+"#$ -o /dev/null\n"
        input=input+"#$ -e %s/errs/%s.err\n" % (base, prefix)
        # write cmd line for one run
        input=input+line+"\n"

        # log some job info when job is done
        input=input+"echo ============ Job information follows ============\n"
        #input=input+"export\n"
        input=input+"echo QUEUE: $QUEUE\n"
        input=input+"echo HOSTNAME: $HOSTNAME\n"
        input=input+"echo NSLOTS: $NSLOTS\n"
        input=input+"echo PWD: $PWD\n"
        input=input+"echo JOB_NAME: $JOB_NAME\n"
        input=input+"echo JOB_ID: $JOB_ID\n"
        input=input+"echo RESTARTED: $RESTARTED\n"
        input=input+"echo Your command line was: "+line+"\n"
        input=input+"echo =================================================\n"

    else: # for PBS queue system
        cmd = "qsub"
        if (qname != ""):
            cmd = cmd + " -q "+qname
        if (qsub_extras != ""):
            cmd = cmd + " "+qsub_extras
        print "qsub cmd is:", cmd

        input=""
        input=input+"#!/bin/bash\n"
        input=input+"#PBS -N "+jobname+"\n"
        input=input+"#PBS -o outs/%s.out\n" % prefix
        input=input+"#PBS -e errs/%s.err\n" % prefix
        input=input+"#PBS -l select=1:ncpus=1\n"

        fulldir = os.getenv("PWD")+"/"+dir
        input=input+"cd "+fulldir+"\n"
        input=input+line+"\n"

        # log some job info when job is done
        input=input+"echo ============ Job information follows ============\n"
        input=input+"echo PBS_O_HOST: $PBS_O_HOST\n"
        input=input+"echo PBS_O_WORKDIR: $PBS_O_WORKDIR\n"
        input=input+"echo PBS_JOBID: $PBS_JOBID\n"
        input=input+"echo PBS_JOBNAME: $PBS_JOBNAME\n"
        input=input+"echo PBS_QUEUE: $PBS_QUEUE\n"
        input=input+"echo PBS_NODEFILE: $PBS_NODEFILE\n"
        input=input+"echo Your command line was: "+line+"\n"

        input=input+"echo =================================================\n"

    (stdoutdata, stderrdata) = runcmd(cmd, input)
    print "stdoutdata:", stdoutdata
    print "stderrdata:", stderrdata

    if (0): # retry if got an error submitting
        count = 0
        while (len(stderrdata) != 0):
            print "Uh oh! Got something on stderr... sleep and retry!", count
            # e.g., maybe we got a msg on stderr that the queue is full
            count = count + 1
            time.sleep(60)
            (stdoutdata, stderrdata) = runcmd(cmd, input)

            # but don't retry forever
            if (count == 120):
                print "hit count of", count, "exiting."
                sys.exit()


###

if __name__ == "__main__":
    # parse command line
    if (len(sys.argv) != 2):
        print "usage: submit.py <dir>"
        sys.exit()
    # dir is the name of a particular set of runs
    # e.g., I had a set called "predprey2"
    dir = sys.argv[1]

    # This makes sure <workdir>/<dir> exists
    # All output files go into <dir>/*
    # stdout and stderr for into <dir>/outs and <dir>/errs, respectively
    ensure_dir_exists(dir)
    ensure_dir_exists(dir+"/outs")
    ensure_dir_exists(dir+"/errs")

    # do everything in <dir>
    os.chdir(dir)

    # Maybe you want to put your own loops for parameter ranges, in
    # this example we just have one loop 'idx'.
    for idx in range(0, 10): # 10 replicates
        # 'prefix' is the prefix of the output files. You may want
        # to change this to reflect all the parameters of that
        # run.
        prefix = "r"+str(idx)

        # 'line' (below) is the command that will be run, e.g.,
        # your executable with some parameters specific to a given
        # run.
        #
        # You will want to put your own command line parameters
        # here to specify the parameters for this particular
        # run. The example command below just copies (with dd) the
        # system file /proc/loadavg to an output file called
        # <prefix>.log inside the <dir> directory.
        #
        # You may also want to pass <prefix> to your executable so
        # that it can write its output log files with that prefix.
        #
        # The stdout and stderr of each run are captured in
        # <dir>/outs and <dir>/errs.
        #
        # Run 'submit.py foobar' and then look in <dir>,
        # <dir>/errs, and <dir>/outs, to see the results.
        #
        # Note that your line will be executed in <dir>, so if you
        # have an executable, you should refer to it as
        # "../<executable>" or by absolute path.
        line = "time /bin/dd if=/proc/loadavg of="+prefix+".log"

        print "dir:", dir, "prefix:", prefix, "line:", line
        qsub_launch(dir=dir, prefix=prefix, line=line, qname=qname)

