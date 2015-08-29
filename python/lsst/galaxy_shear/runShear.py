from great3sims import constants, run
import argparse
import os
import shutil
import signal
import sys
import time
import subprocess
import lsst.galaxy_shear.analyzeShearTest as analyzeShearTest
import lsst.afw.table
#    input a list of pids which we want to wait on
#    return when all are done, and none has returned an error

def waitforpids(pidlist, waituntil=0):
    while len(pidlist) > waituntil:
        (waitpid,retcode) = os.wait()    
        print "pid %d is done, retcode %d\n" % (waitpid, retcode)
        pidlist.remove(waitpid)
        if retcode != 0:
            for killpid in pidlist:
                os.kill(killpid,signal.SIGQUIT);
            raise StandardError("Forked process failed for pid %d\n" % waitpid)

#    runcmd is used to make external calls.  
#    args is a list of what would be space separated arguments, with the convention that
#    args[0] is the name of the command.  An argument of space separated names should not be quoted
#
#    If stdoutname is set, it is assumed that you want to output to a file of that name.
#    If stdouthandle is set, it is assumed that you want to output to that handle.
#    Otherwise, this command received stdout back through a pipe and returns it.
#
#    On error, this program will throw StandardError with a string containing the stderr return
#    from the calling program.  Throws may also occur directly from process.POpen (e.g., OSError)            
#
#    This program will also return a cmd style facsimile of the command and arguments 

def runcmd(args,env=None,stdoutname=None,stderrname=None,append=True):
    
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

def runShear(baseDirs, test=None, forks=1, clobber=True, great3=False, galsim=False, meas=False, anal=False, nGrow=None, footprintSize=None):
    #   open the run_params file for this run and extract what we need.
    fin = open("run_params", 'r')
    params = dict()
    for line in fin:
        if line.isspace():
            continue
        inp = line.split()
        if inp[0][0] == '#':
            continue
        params[inp[0]] = inp[1]
    
    n_subfields = int(params["n_subfields"])
    n_epochs = int(params["n_epochs"])
    ndims = int(params["ndims"])
    galaxy_stamp_size = int(params["galaxy_stamp_size"])
    gal_dir = os.path.abspath(params["gal_dir"])
    fov = float(params["fov"])
    exp_type = params["exp_type"]
    output_dir = params["output_dir"]
    if "psf_dir" in params.keys():
        psf_dir = os.path.abspath(params["psf_dir"])
    if "shear_value" in params.keys():
        shear_value = params["shear_value"]
    shape_field = params["shape_field"] 
    #  Do initial directory setup
    pidlist = []
    if great3:
        for base in baseDirs:
                if os.path.isdir(base + "/control"):
                     shutil.rmtree(base + "/control")
                os.makedirs(base + "/control/ground/constant")
                shutil.copy("run_params", os.path.join(base, "run_params"))
                baseArgs = base.split("_")
                if  baseArgs[0][0] == 'f' and baseArgs[0][1:].isdigit() and len(baseArgs) == 3:
                    filter = int(baseArgs[0][1:])
                    seeing = float(baseArgs[1])
                    shear_value = float(baseArgs[2])
                    psf_dir = os.path.join(params["psf_lib_dir"], "f%d_%s"%(filter, seeing))
                    fout = open(os.path.join(base, "run_params"), "a")
                    fout.write("filter %d\n"%filter)
                    fout.write("seeing %s\n"%seeing)
                    fout.write("shear_value %s\n"%shear_value)
                    fout.write("psf_dir %s\n"%psf_dir)
                    fout.close()

                os.symlink(psf_dir, base + "/control/ground/constant/psfs") 
                fout = open("%s/%s/_mapper"%(base, output_dir), "w")
                fout.write("lsst.obs.great3.Great3Mapper")
                fout.close()
                shear_value=float(shear_value)
                constants.image_size_deg = fov
                constants.nrows = ndims
                constants.ncols = ndims
                constants.n_subfields =  n_subfields
                constants.xsize["ground"][True] = galaxy_stamp_size
                constants.xsize["ground"][False] = galaxy_stamp_size
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
                        run(base, gal_dir=gal_dir, steps=["metaparameters", "catalogs", "config",],
                            experiments=[exp_type], obs_type="ground", shear_type=["constant"], 
                            draw_psf_src = '%s/control/ground/constant/psfs/psfs.index'%base, subfield_max=subfield_max,
                            shear_value=shear_value, shear_angle=60)
                        sys.exit(0)
                else:
                        run(base, gal_dir=gal_dir, steps=["metaparameters", "catalogs", "config",],
                            experiments=[exp_type], obs_type="ground", shear_type=["constant"], 
                            draw_psf_src = '%s/control/ground/constant/psfs/psfs.index'%base, subfield_max=subfield_max,
                            shear_value=shear_value, shear_angle=60)
        waitforpids(pidlist)

    if galsim:
        for base in baseDirs:
                cwd = os.getcwd()
                os.chdir(base)
                #runcmd(("galsim", "cgc.yaml", "output.nproc=%d"%forks, "image.nproc=%d"%forks),
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
            tries = 0
            shutil.copy("processShearTest.py", "temp.py")
            fout = open("temp.py", "a")
            fout.write("root.galaxyStampSize = %d\n"%galaxy_stamp_size)
            fout.write("root.measPlugin = '%s'\n"%shape_field)
            fout.write("root.noClobber = %s\n"%(not clobber))
            if not nGrow is None:
                fout.write('root.measurement.plugins["modelfit_CModel"].region.nGrowFootprint=%d\n'%nGrow)
            if not test is None:
                fout.write('root.test="%s"\n'%test)
            fout.close()
            out = os.path.join(base, output_dir)
            while tries < 5:
                try:
                    runcmd(("processShearTest.py", out, "-j", str(forks), "--configfile=temp.py",
                        "--id", "subfield=0..%d"%(n_subfields-1), "epoch=0..%d"%(n_epochs-1),
                        "--output", out
                        ), stdoutname="%s/meas.stdout"%base, append=False)
                    break
                except:
                    tries = tries + 1
                    continue
    if anal:
        if test is None:
            outfile = "anal.fits"
        else: 
            outfile = "anal_%s.fits"%test
        for base in baseDirs:
                baseArgs = base.split("_")
                if  baseArgs[0][0] == 'f' and baseArgs[0][1:].isdigit() and len(baseArgs) == 3:
                    filter = int(baseArgs[0][1:])
                    seeing = float(baseArgs[1])
                    shear_value = float(baseArgs[2])
                    psf_dir = os.path.join(params["psf_lib_dir"], "f%d_%s"%(filter, seeing))
                    params["filter"] = filter
                    params["seeing"] = seeing
                    params["shear_value"] = shear_value
                    params["psf_dir"] = psf_dir
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
                        analyzeShearTest.runAnal(base, os.path.join(base, outfile), os.path.join(base, "sum_"+outfile), params, test=test)
                        sys.exit(0)
                else:
                        analyzeShearTest.runAnal(base, os.path.join(base, outfile), os.path.join(base, "sum_"+outfile), params, test=test)
        if forks > 1:
            waitforpids(pidlist)
        outCat = None

        for base in baseDirs:
            sourceCat = lsst.afw.table.BaseCatalog.readFits(os.path.join(base, "sum_"+outfile))
            if outCat is None:
                outCat = lsst.afw.table.BaseCatalog(sourceCat.getSchema())
            outCat.append(sourceCat[0])
        print "outCat: ", len(outCat)
        if not outCat is None: 
            outCat.writeFits(os.path.join("sum_"+outfile))
            
            
        

if __name__ == "__main__":
    """
    This main program runs a single directory containing a run_params

    basedir               runs from the specified directory
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("base_dir", help="name of the directory to run from", type = str)
    parser.add_argument("-3", "--great3", action='store_true', help="run great3sims")
    parser.add_argument("-g", "--galsim", action='store_true', help="run galsim")
    parser.add_argument("-m", "--meas", action='store_true', help="run measurement algorithm")
    parser.add_argument("-a", "--anal", action='store_true', help="run analysis program")
    parser.add_argument("-t", "--test", type=str, help="Name of the test", default=None)
    parser.add_argument("-n", "--nGrow", type=int, help="Grow footprints", default=None)
    parser.add_argument("-s", "--footprintSize", type=int, help="set footprint to a square of this size", default=None)
    parser.add_argument("-c", "--clobber", type=int, help="Delete or keeps existing files?", default=1)
    args = parser.parse_args()

    #  Test names automatically set the corresponding argument
    if args.test[0] == "n" and args.nGrow is None:
        nGrow = int(args.test[1:])
    else:
        nGrow = args.nGrow
    if args.test[0] == "s" and args.footprintSize is None:
        footprintSize = int(args.test[1:])
    else:
        footprintSize = args.footprintSize
        
    runShear((args.base_dir,), clobber=args.clobber, great3=args.great3, galsim=args.galsim, meas=args.meas, anal=args.anal, test=args.test, nGrow=nGrow, footprintSize=footprintSize)
