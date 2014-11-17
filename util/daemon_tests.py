import daemon
import time
import hgapi as hg
import os, sys
import shutil as sh
from datetime import datetime
import glob
from subprocess import Popen, PIPE

'''
Daemon automatic NPKv2 Tests.
Installation:
    install "python-daemon" from Pypi.
    install "hgapi" with easy_install.
Running:
    For launching the program open a Terminal and write "python daemon_NPKv2_Tests.py start"
    For stopping, write "python daemon_NPKv2_Tests.py stop"
    Apparently you can run it directly from Textmate.
Result:
    You obtain every time for the hour defined by the variable 'hour_for_report' a report for the automatic Tests
    at the address defined by 'path_root'
    
'''
                                                    
class automatic_tests(object):
    '''
    '''
    def __init__(self):
        self.hour_for_report = 16                                                                       # hour at which we want the Tests to be done.
        self.date_old = ''                                                                              # initializing                                                                    
        self.path_root = '/Users/chiron/Desktop/encours/Divers/daemonPython'                            # Path for stocking the report
        self.path_draft = os.path.join(self.path_root, 'draft')                                         # address of the repository once downloaded.
        self.address_repository = 'https://Yo@bitbucket.org/NPK_v2/draft'                               # address of the repository to test.
        self.delay_after_clone = 15                                                                     # delay after hg clone. 

    def prepare(self):
        try:
            print "erasing old NPKv2"
            sh.rmtree(self.path_draft)                                                                  # removes the previous NPKv2 directory if still there.
        except:
            print "no existing NPKv2 directory here."
        os.mkdir(self.path_draft)
        ##### Clone
        hg.hg_clone(self.address_repository, self.path_draft)                                           # Retrieves automatically NPKv2
        time.sleep(self.delay_after_clone)                                                              # waiting to be sure to have all
        os.chdir(self.path_draft)                                                                       # Getting positionned in the NPKv2 directory for executing the file Tests.py
    
    def Test_and_Report(self):
        ##### Makes the Tests                                                       
        (self.sout, self.serr) = Popen(["python", "autom_Tests.py"], stdout = PIPE).communicate()
        ##### Writes a report for stoud
        self.report_stout = os.path.join(self.path_root, "Tests_report.txt")                            # address for the report
        try:
            os.remove(self.report_stout)
        except:
            print "no existing report_stout"
        with open(self.report_stout, "w") as test_report:
            test_report.write(self.sout)
        self.date_old = self.date
        
    def do_NPKv2_tests(self):
        '''
        Loads and performs Tests.py on NPKv2.
        '''
        while True:
            os.chdir(self.path_root)                                                                                    # get positionned in the directory for loading a new NPKv2 directory.
            hour = 0                                                      
            while hour < self.hour_for_report or hour > self.hour_for_report:                                           # check if the hour is OK
                hour = int(datetime.strftime(datetime.now(), '%H'))
                time.sleep(10)                                                                                          # waiting to avoid the processor to work too much.
            date = datetime.strftime(datetime.now(), '%Y-%m-%d')
            if self.date != self.date_old:                                                                              # if the date is new
                self.prepare()
                self.Test_and_Report()

def run():
    at = automatic_tests()
    with daemon.DaemonContext():
        at.do_NPKv2_tests()

if __name__ == "__main__":
    run()