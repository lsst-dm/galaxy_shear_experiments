import argparse
import glob
import os
import shutil
import signal
import sys
import time
import subprocess


import lsst.pex.config as pexConfig
import lsst.afw.table

from great3sims import constants, run
from lsst.galaxy_shear.shearConfig import RunShearConfig
from lsst.galaxy_shear.analyzeShearTest import runAnal

def waitforpids(pidlist, waituntil=0):
    """
    input a list of pids which we want to wait on
    return when all are done, and none has returned an error
    """

    while len(pidlist) > waituntil:
        (waitpid,retcode) = os.wait()    
        print "pid %d is done, retcode %d\n" % (waitpid, retcode)
        pidlist.remove(waitpid)
        if retcode != 0:
            for killpid in pidlist:
                os.kill(killpid,signal.SIGQUIT);
            raise StandardError("Forked process failed for pid %d\n" % waitpid)

def runcmd(args,env=None,stdoutname=None,stderrname=None,append=True):
    """
    runcmd is used to make external calls.  
    args is a list of what would be space separated arguments, with the convention that
    args[0] is the name of the command.  An argument of space separated names should not be quoted

    If stdoutname is set, it is assumed that you want to output to a file of that name.
    If stdouthandle is set, it is assumed that you want to output to that handle.
    Otherwise, this command received stdout back through a pipe and returns it.

    On error, this program will throw StandardError with a string containing the stderr return
    from the calling program.  Throws may also occur directly from process.POpen (e.g., OSError)            

    This program will also return a cmd style facsimile of the command and arguments 
    """
    
    cmdstring = args[0]
    for i in range(1,len(args)):
        if args[i].find(' ') >= 0: cmdstring += " '%s'" % args[i]
        else: cmdstring += " %s" % args[i]
    if stdoutname != None:
        cmdstring += " > %s" % stdoutname
    infostring = '%s: running "%s"\n' % (time.asctime(time.localtime()),cmdstring)

    errcode = 0
    # actually do it
    if stdoutname != None:
        if append:
            stdout_handle = open(stdoutname,"a")
        else:
            stdout_handle = open(stdoutname,"w")
    else:
        stdout_handle = None 

    if stderrname != None:
        if append:
            stderr_handle = open(stderrname,"a")
        else:
            stderr_handle = open(stderrname,"w")
    else:
        stderr_handle = None

    #    Call the requested command, wait for its return, then get any return pipe output
    proc = subprocess.Popen(args,stdout=stdout_handle,stderr=stderr_handle,close_fds=False,env=env)
    errcode = proc.wait()

    if stdoutname != None: stdout_handle.close()
    if stderrname != None: stderr_handle.close()

    #    The calling program will not throw on failure, but we should to make it compatible with
    #    Python programming practices
    if errcode:
        errstring = 'Command "%s" failed with code %d\n' % (cmdstring,errcode) 
        raise StandardError(errstring)
    return

def runShear(baseDirs, test=None, forks=1, clobber=1, great3=False, galsim=False, meas=False, anal=False):
    #   open the run_params file for this run and extract what we need.
    
    pidlist = []
    if great3:
        for base in baseDirs:
                config = RunShearConfig()
                config.load(os.path.join(base, "shear.config"))
                great3_dir = base + "/" + config.exp_type + "/ground/constant"
                if not os.path.isdir(great3_dir):
                    os.makedirs(great3_dir)
                baseArgs = base.split("_")
                if  baseArgs[0][0] == 'f' and baseArgs[0][1:].isdigit() and len(baseArgs) == 3:
                    config.filter = int(baseArgs[0][1:])
                    config.seeing = float(baseArgs[1])
                    config.shear_value = float(baseArgs[2])
                    config.psf_dir = os.path.join(config.psf_lib_dir,
                                        "f%d_%s"%(config.filter, config.seeing))
                    if not os.path.isdir(os.path.join(great3_dir, "psfs")):
                        os.symlink(config.psf_dir, os.path.join(great3_dir, "psfs")) 
                else:
                    os.symlink(os.path.abspath(config.psf_dir), os.path.join(great3_dir, "psfs")) 
                fout = open(os.path.join(great3_dir, "_mapper"), "w")
                fout.write("lsst.obs.great3.Great3Mapper")
                fout.close()
                constants.image_size_deg = config.fov
                constants.nrows = config.ndims
                constants.ncols = config.ndims
                constants.n_subfields =  config.n_subfields
                constants.xsize["ground"][True] = config.galaxy_stamp_size
                constants.xsize["ground"][False] = config.galaxy_stamp_size
                constants.n_subfields_per_field["constant"][True] = 1
                constants.subfield_grid_subsampling = 1
                constants.n_deep_subfields = 0
                constants.deep_frac = 0.0
                subfield_max = constants.n_subfields + constants.n_deep_subfields - 1
                if forks > 1:
                    if len(pidlist) >= forks:
                        waitforpids(pidlist, waitutil=forks-1)
                    pid = os.fork()
                    if pid:
                        # we are the parent
                        pidlist.append(pid)
                    else:
                        # we are the child
                        run(base, gal_dir=config.gal_dir, steps=["metaparameters", "catalogs", "config",],
                            experiments=[config.exp_type], obs_type="ground", shear_type=["constant"], 
                            draw_psf_src = '%s/control/ground/constant/psfs/psfs.index'%base, 
                            subfield_max=subfield_max,
                            shear_value=config.shear_value,
                            shear_angle=config.shear_angle
                        )
                        sys.exit(0)
                else:
                    run(base, gal_dir=config.gal_dir, steps=["metaparameters", "catalogs", "config",],
                        experiments=[config.exp_type], obs_type="ground", shear_type=["constant"], 
                        draw_psf_src = '%s/control/ground/constant/psfs/psfs.index'%base, 
                        subfield_max=subfield_max,
                        shear_value=config.shear_value,
                        shear_angle=config.shear_angle
                    )

    waitforpids(pidlist)

    if galsim:
        for base in baseDirs:
                config = RunShearConfig()
                config.load(os.path.join(base, "shear.config"))
                cwd = os.getcwd()
                os.chdir(base)
                print "output.noclobber=%s"%(not clobber)
                runcmd(("galsim", "cgc.yaml", "output.noclobber=%s"%(not clobber)),
                        stdoutname="galsim.stdout", append=False)
                runcmd(("galsim", "cgc_psf.yaml", "output.noclobber=%s"%(not clobber)),
                        stdoutname="galsim.stdout", append=True)
                runcmd(("galsim", "cgc_star_test.yaml", "output.noclobber=%s"%(not clobber)),
                        stdoutname="galsim.stdout", append=True)
                os.chdir(cwd)

    if meas:
        for base in baseDirs:
            config = RunShearConfig()
            config.load(os.path.join(base, "shear.config"))
            tempPath = os.path.join(base, "processShearTest.py") 
            shutil.copy("processShearTest.py", tempPath)
            fout = open(tempPath, "a")
            fout.write("root.galaxyStampSize = %d\n"%config.galaxy_stamp_size)
            fout.write("root.measPlugin = '%s'\n"%config.shape_field)
            fout.write("root.noClobber = %s\n"%(not clobber))
            if not test is None and test[0] == 'n':
                nGrow = int(test[1:])
                fout.write('root.measurement.plugins["modelfit_CModel"].region.nGrowFootprint=%d\n'%nGrow)
            if not test is None:
                fout.write('root.test="%s"\n'%test)
            if not test is None and test[0] == 's':
                footprintSize = int(test[1:])
                fout.write('root.measurement.plugins["modelfit_CModel"].region.nInitialRadii=0\n')
                fout.write('root.footprintSize=%d\n'%footprintSize)
            if not test is None:
                fout.write('root.test="%s"\n'%test)
            fout.close()
            out = os.path.join(base, config.exp_type + "/ground/constant")
            runcmd(("processShearTest.py", out, "-j", str(forks), "--configfile=%s"%(tempPath),
                        "--id", "subfield=0..%d"%(config.n_subfields-1), "epoch=0..%d"%(config.n_epochs-1),
                        "--output", out
                        ), stdoutname="%s/meas.stdout"%base, append=False)
    if anal:
        if test is None:
            outfile = "anal.fits"
        else: 
            outfile = "anal_%s.fits"%test
        for base in baseDirs:
                config = RunShearConfig()
                config.load(os.path.join(base, "shear.config"))
                config.psf_dir = os.path.join(config.psf_lib_dir, "f%d_%s"%(config.filter, config.seeing))
                if os.path.isfile(os.path.join(base, outfile)) and not clobber:
                    continue
                if forks > 1:
                    if len(pidlist) >= forks:
                        waitforpids(pidlist, waitutil=forks-1)
                    pid = os.fork()
                    if pid:
                        # we are the parent
                        pidlist.append(pid)
                    else:
                        # we are the child
                        runAnal(base, outfile, config, test=test)
                        sys.exit(0)
                else:
                        runAnal(base, outfile, config, test=test)
        if forks > 1:
            waitforpids(pidlist)
        outCat = None

        for base in baseDirs:
            sourceCat = lsst.afw.table.BaseCatalog.readFits(os.path.join(base, "anal_subfield.fits"))
            if outCat is None:
                outCat = lsst.afw.table.BaseCatalog(sourceCat.getSchema())
            print "catalog for %s has %d rows"%(base, len(sourceCat))
            outCat.extend(sourceCat)
        print "outCat: ", len(outCat)
        if not outCat is None: 
            outCat.writeFits(os.path.join("sum_"+outfile))
            print "output catalog %d rows"%len(outCat)
            
if __name__ == "__main__":
    """
    This main program reads data for multiple runs in filter, seeing, and shear value, and

    filter                runs with the specified filter
    seeing_values         runs with the specified seeing
    shear_values          runs with the specified shear values
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--base", help="directory containing shear.config", type = str)
    parser.add_argument("-v", "--shear_values", help="comma separated shear values",
        type = str, default=None)
    parser.add_argument("-f", "--filter", help="number of the filter for PhoSim", type = int)
    parser.add_argument("-s", "--seeing_values", help="comma separated seeing values",
        type = str, default=None)
    parser.add_argument("-c", "--clobber", type=int, help="set to 0 to preserve existing outputs",
        default=1)
    parser.add_argument("-t", "--test", help="run on specific data subset")
    parser.add_argument("-p", "--processes", help="Number of forks, max", type = int, default=1)
    parser.add_argument("-3", "--great3", action='store_true', help="run great3sims")
    parser.add_argument("-g", "--galsim", action='store_true', help="run galsim")
    parser.add_argument("-m", "--meas", action='store_true', help="run measurement algorithm")
    parser.add_argument("-a", "--anal", action='store_true', help="run analysis program")

    args = parser.parse_args()

    #   open the config file for this run
    config = RunShearConfig()
    if args.base is None:
        config.load("shear.config")
    else:
        config.load(os.path.join(args.base, "shear.config"))
    if not args.filter is None:
        config.filter = args.filter

    if args.base is None and args.shear_values is None:
        print "Must specify either a base directory or a set of shear values."
        sys.exit(1)

    baseDirs = []
    #  Create directories named with filter, seeing_values, and shear_values.
    if args.base is None:
        seeing_values = args.seeing_values
        if seeing_values is None:
            seeing_values = str(config.seeing)

        for seeing in seeing_values.split(','):
            config.seeing = float(seeing)
            if args.shear_values == "all":
                for dir in glob.glob("f%d_%s_.[0-9]*"%(config.filter, seeing)):
                    baseDirs.append(dir)
            else:    
                for shear_value in args.shear_values.split(","):
                    base = 'f%d_%s_%s'%(config.filter, seeing, shear_value)
                    baseDirs.append(base)

        for base in baseDirs:
            config.shear_value = float(base.split("_")[2])
            config.seeing = float(base.split("_")[1])
            if not os.path.isdir(base):
                os.makedirs(base)
            if not os.path.isfile(os.path.join(base, "shear.config")):
                config.save(os.path.join(base, "shear.config"))

    #  If the base directory is specified, assume that a single seeing and shear
    #  value have been specified, or that they are in the default config
    #  Otherwise, an error will occur.
    else:
        if not args.shear_values is None:
            config.shear_value = float(args.shear_values)
        if not args.seeing_values is None:
            config.seeing_value = float(args.seeing_values)
        if not os.path.isdir(args.base):
            os.makedirs(args.base)
        if not os.path.isfile(os.path.join(base, "shear.config")):
            config.save(os.path.join(base, "shear.config"))
        baseDirs.append(args.base)

    runShear(baseDirs, clobber=args.clobber, forks=args.processes, great3=args.great3,
             galsim=args.galsim, meas=args.meas, anal=args.anal, test=args.test)
