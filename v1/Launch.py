#!/usr/bin/env python 
# encoding: utf-8

"""
Processing of NMR data-sets

NPK - the python front end for processing NMR data-sets
=======================================================

USAGE
    i)    interactive mode
        NPK
                to start an interact session
    or
    ii)   script execution
        NPK -e script.py
                to execute in NPK the script : script.py
    or
    iii)  NPK help
                to show this documentation
    or
    iv)   inline commande
        NPK action input_datafile output_datafile { processing_parameterfile } { additional parameters }
        where
            action is the processing to apply
            input_datafile and output_datafile are gifa binary files
            processing_parameterfile contains all the needed parameters for the processing
                is optionnal, if omited, default values will be systematically used
                however, must be present (even empty) if additional parameters are to be used.

    All python commands, as well as NPK commands are available in the four types of usage.

    Note :
    i) a regular python 2.1 console, type ^D to exit (^Z in Windows)

    ii) a regular python 2.1 interpretor

    iv) is the standard way of invoking NMR processing



IN-LINE COMMAND mode:
====================

NPK action input_datafile output_datafile { processing_parameterfile } { additional parameters }

ACTION
    possible actions are :
        FT1D         - perfoms 1D Fourier Analysis
                       see details in Process1D documentation
        MaxEnt1D     - perfoms 1D Maximum Entropy Analysis
                       see details in Process1D documentation
        FT2D         - perfoms 2D Fourier Analysis
                       see details in Process2D documentation
        MaxEnt2D     - perfoms 2D Maximum Entropy Analysis
                       see details in Process2D documentation
        PreDOSY2D    - prepares for DOSY analysis
                       see details in Process2D documentation
        DOSY2D       - performs a DOSY analysis
                       see details in ProcessDosy documentation
        PDOSY2D      - same as DOSY2D; but will start a parrallel process instead
                       see details in batch_dosy documentation
        FT3D         - perfoms 2D Fourier Analysis
                       see details in Process3D documentation
        ImportBruker - Import from Bruker XWinNMR/TopSpin format, a file or a set of files
        ImportVarian - Import from Varian VNMR format, a file or a set of files
        CHECK        - no action, permits to chek syntax
        help         - this help file

NMR BINARY FILES
    NPK reads only Gifa files.
    These files are read and created by the various versions of the Gifa program
        Gifa v4, see http://abcis.cbs.cnrs.fr
    and Gifa v5, see http://www.NMRtec.com
    The programm NMRNoteBook, developped by NMRtec is also able to export and import Gifa files
    
    File extensions :
    Back in Gifa v4 time, no extension was used to specify the NMR Gifa binary files.
    With Gifa V5 a new convention has been set-up for the file extensions :
            time domain     freq domain
    1D          .gf1            .gs1
    2D          .gf2            .gs2
    3D          .gf3            .gs3
    These extensions are now enforced by NPK
    It should be noted that no provision is made to specify a file in an intermediate state (processed in F2 but not in F1)
    it is up to you name these files depending on what you want to do with them afterwards.

    If you want to import NMR binary files from you spectrometer, you can
    - either use one of the Gifa v4 importers found on the WEB site (the util sub directory)
    - use the provided python importers : ImportBruker, ImportVarian
    - use any other importer

PARAMETER FILES

    The NMR processing is controled by parameters described in the respective documentation.
    Each and every parameter has a default value which is not too silly, so usually, giving no parameter may give you some result.
    However, you will have to modify the processing results by changing these parameters.

    Parameters are given either through a parameter file, or by inline commands.
    parameter files have the .gtb extention
    parameter files contains parameters in the form of a property list :

        #some comment
        param1=value string
        param2=value string    ... etc ...

    Addtional parameters can be given to the processing, independently to the parameter file
    by directly giving the parameter on the calling line, in the form additionnal_param=value
    
    WINDOWS users:  = is used as a separator in windows, so you will have to type "additionnal_param=value"
    
    additional parameterfiles are also associated to the binary data files
        input_datafile.gtb
            should be present near the input_datafile
            and contains details on the data not found in the binary gifa file format.
        output_datafile.gtb
            is created by the processing and contains processing details

EXAMPLES
========
    NPK FT1D input.gf1 output.gs1 

        launch the processing of a 1D FID, called input.gf1, and create the file output.gs1
        if the file is present, it is used as the processing parameter file, otherwise default parameters are used

    NPK FT2D input.gf2 output.gs2 processing.gtb

        launch the processing of a 2D experiment, called input.gf2, and create the 2D spectrum output.gs2
        the file processing.gtb is used as the processing parameter file.
        Note that this file may be empty.

    NPK DOSY2D Post_F2.gs2 DOSY.gs2 processing.gtb audittrail=test1.html  me_details=1  me_size=512

        launch the processing of a DOSY, on a dataset called Post_F2.gs2, and create the 2D DOSY spectrum DOSY.gs2
        the file processing.gtb is used as the processing parameter file.
        additionnal parameters are used :
            audittrail=test1.html overwrite the default audittrail file name
            me_details=1  me_size=512   overwrite the size of the Inverse Laplace processing defined in the processing.gtb file

STARTING NPK  & POSSIBLE ERRORS
    NPK requires a certain number of environment files to be present.
    $JAVA and $NPK_HOME are defined in the NPK shell script.

    The following sequence of event takes place when you start NPK 
    - a java virtual machine (JVM) must be available (found in $JAVA)
    - the jython code is launched in the JVM (found in $NPK_HOME/Java)
    - additional jython libraries are loaded (found in $NPK_HOME/python/Lib)
    - NPK code is initialized
        - loading the binary libraries (found in $NPK_HOME/Java)
        - loading and compiling the NPK python library (found in $NPK_HOME/python/NPK)
    - command line is parsed and processing starts ( found in $NPK_HOME/python/NPK/Launch.py )
    
See also:
    some additional code is provided in the Console code which is loaded when entering the interactive mode :
    
    - status()     : displays the internal state of the Kernel
    - import Graph : opens a minimum graphic interface to the 1D buffer (resize to redraw)


"""

from __future__ import print_function

__author__ = "Marc A. Delsuc <delsuc@igbmc.fr>"
__date__ = "Oct 2009"


import sys
import time
import re
import os
#from NPK import *                              #Commented when moved to v1 folder
from Generic import *                           #Added when copied in v1 folder
import Param

def status():
    """
    print a summary of the internal state of the kernel
    """
    d = get_dim()
    if d == 1:
        d1 = "- current working buffer"
    else:
        d1 = ""
    if d == 2:
        d2 = "- current working buffer"
    else:
        d2 = ""
    if d == 3:
        d3 = "- current working buffer"
    else:
        d3 = ""
    dim(1)
    com_max()
    report = """
DIM 1 %s
=====
   buffer size : %i - itype %i
   values from : %f to %f
Spectral width : %f
%i peak(s) in database

"""%(d1,get_si1_1d(),get_itype_1d(),geta_max(2),geta_max(1),get_specw_1d(),get_npk1d())

    dim(2)
    com_max()
    report = report+"""
DIM 2 %s
=====
    buffer sizes : %i x %i - itype %i
     values from : %f to %f
Spectral widthes : %f x %f
%i peak(s) in database

"""%(d2,get_si1_2d(),get_si2_2d(),get_itype_2d(),geta_max(2),geta_max(1),get_specw_1_2d(),get_specw_2_2d(),get_npk2d())
    dim(3)
    com_max()
    report = report+"""
DIM 3 %s
=====
    buffer sizes : %i x %i x %i - itype %i
     values from : %f to %f
Spectral widthes : %f x %f x %f
%i peak(s) in database

"""%(d3,get_si1_3d(),get_si2_3d(), get_si3_3d(),get_itype_3d(),geta_max(2),geta_max(1),get_specw_1_3d(),get_specw_2_3d(),get_specw_3_3d(),get_npk3d())
    print("\n"+report)
    dim(d)



def ImportBruker(audit,inputfilename,outputfilename):
    """Import Bruker data-sets
    
    inputfilename can be either a file (fid or ser) or a directory
    """
    import Bruker
    import os.path as op
    if (op.isfile(inputfilename)):
        audittrail( audit, "phase", "Importation from Bruker file")
        Bruker.Import(audit,inputfilename,outputfilename)
    elif (op.isdir(inputfilename)):
        if (outputfilename != "-"):
            print("WARNING output file name is ignored when importing series")
        audittrail( audit, "phase", "Importation Bruker files from directory")
        audittrail( audit, "text", "source directory"," ",inputfilename)
        Bruker.Import(audit,inputfilename,"-")
    else:
        raise inputfilename+" : file not found or not available"

def ImportVarian(audit,inputfilename,outputfilename):
    """Import Varian data-sets
    
    inputfilename can be either a file (fid) or a directory
    """
    import Varian
    import os.path as op
    if (op.isfile(inputfilename)):
        audittrail( audit, "phase", "Importation from Varian file")
        Varian.Import(audit,inputfilename,outputfilename)
    elif (op.isdir(inputfilename)):
        if (outputfilename != "-"):
            print("WARNING output file name is ignored when importing series")
        audittrail( audit, "phase", "Importation Varian files from directory")
        audittrail( audit, "text", "source directory"," ",inputfilename)
        Varian.Import(audit,inputfilename,"-")
    else:
        raise filename+" : file not found or not available"

def PDOSY2D(phase,audit,p_in,inputfilename,outputfilename):
    """
    calls the routine to realise a processing in parrallel
    phase : additional parameter LAUNCH / FINAL
    """
    import batch_dosy
    try:
        nchunk = p_in["N_JOBS"]
    except:
        nchunk = batch_dosy.N_JOBS
    # launch phase
    if phase == "LAUNCH":
        audittrail( audit, "phase", "LAUNCHING phase")
        audittrail( audit, "text", "Job details",
            "type of jobs",batch_dosy.BATCH,
            "nbre of jobs",nchunk,
            )
    if phase == "FINAL":
        audittrail( audit, "phase", "FINAL phase")
    batch_dosy.process(phase,inputfilename,outputfilename,L_paramfilename,nchunk)

    
def actionparser(action, audit,p_in,f_in,f_out,inputfilename,outputfilename):
            """ action parser utility"""
            if (action == "FT1D"):
                FT1D(audit,p_in,f_in,f_out,inputfilename,outputfilename)
            elif (action == "MaxEnt1D"):
                MaxEnt1D(audit,p_in,f_in,f_out,inputfilename,outputfilename)
            elif (action == "FT2D"):
                FT2D(audit,p_in,f_in,f_out,inputfilename,outputfilename)
            elif (action == "MaxEnt1D"):
                print("Sorry, MaxEnt1D Not Fully implemented yet")
            elif (action == "MaxEnt2D"):
                MaxEnt2D(audit,p_in,f_in,f_out,inputfilename,outputfilename)
            elif (action == "PreDOSY2D"):
                PreDosy2D(audit,p_in,f_in,f_out,inputfilename,outputfilename)
            elif (action == "DOSY2D"):
                Dosy2D(audit,p_in,f_in,f_out,inputfilename,outputfilename)
            elif (action == "PDOSY2D"):
                PDOSY2D("LAUNCH",audit,inputfilename,outputfilename)
            elif (action == "PDOSYCollect"):
                PDOSY2D("FINAL",audit,inputfilename,outputfilename)
            elif (action == "ImportBruker"):
                ImportBruker(audit,inputfilename,outputfilename)
            elif (action == "ImportVarian"):
                ImportVarian(audit,inputfilename,outputfilename)
            elif (action == "FT3D"):
                FT3D(audit,p_in,f_in,f_out,inputfilename,outputfilename)
            elif (action == "CHECK" or action == "check" ):
                print("CHECK command")
                print(("action: %s  input: %s  output: %s  param: %s"%(action,inputfilename,outputfilename,paramfilename)))
                sys.exit()
                
            else:
                print("Action -" + action + "- unknown\n")
                print("Possible actions : help CHECK FT1D FT2D FT3D PreDOSY2D DOSY2D PDOSY2D MaxEnt1D MaxEnt2D ImportBruker ImportVarian\n")

def auditname(filename,action):
        """ Creates standard name for audit trail
        """
        dir = os.path.dirname(filename)
        base = os.path.basename(filename)
        if base == ".": # happens !
            base = "current"
        return  os.path.join(dir,"NPK_audit_" + action + "_" + base + ".html")

if __name__ =="__main__":
# WARNING, because of deficiency in jython, we do not define a main() here.
# this implies that every thing here is run in the global name scope.
# so you might have interferences between name defined here, and name defined in external scripts
# so all variable are named L_xx to reduce risks
    #    sys.argv.pop(0)
    L_argv=[]
    for i in sys.argv:
        L_argv.append(i)
    # pop the first arg, to remove the name of this script
    L_argv.pop(0)
    L_narg = len(L_argv)

    # this permit to help the lay user to exit when an error occurs
    sys.ps1="*** an error occured in NPK; ( Type ^D (unix) or ^Z (windows) to exit) >"
    
    # NO ARG case, => starts interactive session
    if ( L_narg == 0 ):
        sys.ps1="NPK>"
        sys.ps2="NPK..."
        print("""
        NPK interactive shell
        Type ^D (unix) or ^Z (windows) to exit
        """)
    else:

    # get 1st arg
        L_arg1=L_argv[0]

    # if -e script - starts script
        if ( L_arg1 == '-e' ):
            try:
                L_script = L_argv[1]
                L_argv.pop(0)
                sys.argv.pop(0)     # remove two entries from sys.argv, to keep same arg offset wether called from NPK or from python
                sys.argv.pop(0)
            except:
                print("script name missing")
                sys.exit(1)
            L_auditfile = auditname(L_script,"script")
            L_audit = auditinitial(L_auditfile,"NPK processing script : "+L_script)
            L_time = time.clock()
            try:
                execfile(L_script)
# script is run here in the same global name scope, so it access the variable defined here, to reduce risks, we preset then with L_
            except:
                sys.excepthook(sys.exc_info()[0],sys.exc_info()[1],sys.exc_info()[2])
                print("\nscript "+L_script+"  failed.  exiting..")
                sys.exit(1)
            L_time = time.clock()-L_time
            audittrail( L_audit, "text", "total elapsed time", "time",str(L_time)+" sec")
            audittrail( L_audit, "close")
            L_audit.close()

            print("total elapsed time", "time:" ,str(L_time)+" sec")
            sys.exit(0)


    # if HELP, show doc
        if ( L_arg1 == "help" or L_arg1 == "HELP" or L_arg1 == "--help" ):
            print(__doc__)
            sys.exit(0)

    # else, parse action line
        try:
            L_action = L_argv.pop(0)
            L_inputfilename = L_argv.pop(0)
            L_outputfilename = L_argv.pop(0)
            try:
                L_paramfilename = L_argv.pop(0)
            except:
                L_paramfilename="-"
                print("no parameter file")
        except:
            print("\nwrong number of arguments")
            print("\nArguments :")
            try:
                    print(("action: "+str(L_action)))
                    print(("inputfilename: "+str(L_inputfilename)))
                    print(("outputfilename: "+str(L_outputfilename)))
                    print(("paramfilename: "+str(L_paramfilename)))
            except:
                    pass
            print("\ntry >NPK help\n")
            print("wrong number of arguments\n")
            sys.exit(1)

        #start timer
        L_time = time.clock()

        L_auditfile = auditname(L_inputfilename,L_action)

        # before opening, check that the auditfile key is not present, and use if it's there
        for L_i in range(1,len(L_argv)):
                if L_argv[L_i].startswith("auditfile="):
                    L_l = re.split(r'(?<!\\)=', L_argv[i],1)   #matches = but not \=
                    L_auditfile = re.sub(r'\\=', '=',L_l[1])      # replaces \= by =

        print("Audit file is : ", L_auditfile)
        L_audit = auditinitial(L_auditfile,"NPK processing of "+str(L_action))
        L_f_out={}

#        try:        # this try - finally is to insure that the audit file will be closed !
        if (1==1):  # if you want to debug if the NPK code, you'd rather comment the try above, as well as the finally further down, and use this if

            # load parameters
            if L_paramfilename == "-":
                audittrail( L_audit, "text", "open default processing parameters")
                L_p_in={}
            else:
                audittrail( L_audit, "phase", "Loading parameters") 
                L_p_in = dict_load(L_paramfilename)                # load p_in as a regular dict, will be copied to the NPKParam when building default values
                audittrail( L_audit, "text", "open processing parameters","parameter file",L_paramfilename)

            # gather additionnal parameters on the calling line, in the form key=value
            while L_argv:
                L_k = L_argv.pop(0)
                try:
                    (L_dkey,L_fval) = Param.parse(L_k)
                except:
                    raise L_k+': error in parameter syntax, should be key=value'
                print("Additional property :", L_dkey, "=", L_fval)
                try:            # check key existence
                    L_val=L_p_in[L_dkey]
                except:
                    pass
                else:
                    print("WARNING, key -",L_dkey,"- defined twice")
                    print("         previous value :",L_val)
                L_p_in[L_dkey]=L_fval
                audittrail( L_audit, "text", "additional parameter", L_dkey,L_fval)

            audittrail( L_audit, "text", "open input file parameters")
            L_f_in=Param.NPKParam()
            try:
                L_f_in.load(L_inputfilename+".gtb")
            except:
                print("Warning, Input data file property empty")
            
            actionparser(L_action, L_audit, L_p_in, L_f_in, L_f_out, L_inputfilename, L_outputfilename)  # launch the job !

#        finally:    # finally, dump parameters and exit
            audittrail( L_audit, "phase", "End of processing")
            if (L_f_out != {}):
                try:
                    L_f_out.dump(L_outputfilename+".gtb")
                    audittrail( L_audit, "text", "write output file parameters", "file name", L_outputfilename+".gtb")
                except:
                    pass
            L_time = time.clock()-L_time
            audittrail( L_audit, "text", "total elapsed time", "time",str(L_time)+" sec")
            audittrail( L_audit, "close")
            L_audit.close()

            print("total elapsed time", "time:" ,str(L_time)+" sec")
            sys.exit()

