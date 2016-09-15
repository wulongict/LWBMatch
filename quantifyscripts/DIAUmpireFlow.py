#!/usr/bin/python
# Python 2.7 Required
# Workflow: DIA-Umpire  --> Comet --> Xinteract
# Input:
# An mzXML file
# A Fasta file in the same folder of mzXML file
#
# Output:
# The mgf file of three quality tiers by DIA-Umpire;
# The pep.xml files of identifications by Comet;
# The ipro.pep.xml files by xinteract
#
#
# -----------------------
# Author: L. WU
# Date: Aug. 18, 2016
# ------------------------

from indexmzXML import *
from CometSearch import *
from Xinteract import *
from WorkflowController import *
from FastaReader import *
import glob


class DB:
    def __init__(self, wfc):
        """
        :type wfc: WorkflowController
        :param wfc: the work flow param
        """
        print '[Info] Creating Database...'
        self.wfc = wfc

    def Run(self, fastapath, wp):
        dbname = self.get_target_decoy_db_from(fastapath, wp)
        return dbname

    def create_target_decoy_db(self, fastafile, wp):
        print '[Info] Creating decoy by perl script'
        [outputpath, outputfile] = os.path.split(fastafile)
        outputpathfile = os.path.join(outputpath, outputfile[0:-6] + 'TargetDecoy.fasta')
        cmdline = 'perl ' + wp.mkdecoyDBperl + ' -f ' + fastafile + ' -o ' + outputpathfile

        if os.path.isfile(outputpathfile) and not self.wfc.run_db_generate:
            return outputpathfile

        fr = FastaReader(fastafile)
        # remove reverse sequence started with ">rev_"
        fr.rewrite()

        print '[Info] >' + cmdline
        os.system(cmdline)
        return outputpathfile

    def get_target_decoy_db_from(self, fastapath, wp):
        print 'Looking for Target Decoy database'
        dbnamelist = glob.glob(os.path.join(fastapath, '*.fasta'))
        if len(dbnamelist) == 1:
            return self.create_target_decoy_db(dbnamelist[0], wp)
        else:
            for eachdbname in dbnamelist:
                if 'TargetDecoy' in eachdbname:
                    return eachdbname
            print 'Could not find TargetDecoy database, please check your mzXML path and make sure the DB or TD-DB is with the mzXML file'
            exit()


class DIAUmpire:
    def __init__(self, wfc, wp):
        print '[Info] Running DIAUmpire...'
        self.wfc = wfc
        self.wp = wp

    def Run(self, wp, inputmzXML):
        print '[Info] Running DIA Umpire: mzXML -> MGF...'

        # in case that the mzXML is not indexed correctly:
        cnvt = Converter(self.wfc,wp)
        cnvt.check_mzxml(inputmzXML)

        mgflist = [inputmzXML[0:-6] + '_Q1.mgf', inputmzXML[0:-6] + '_Q2.mgf', inputmzXML[0:-6] + '_Q3.mgf']
        # self.database = self.get_target_decoy_db_from()

        if len(glob.glob(mgflist[0])) != 0 and self.wfc.run_dia_umpire == False:
            return mgflist
        # self.ResetReRun()

        diaumpireconfigurefile = os.path.join(wp.DIAUmpirePath, 'diaumpire.se_params')
        diaumpirecmd = 'java -jar -Xmx16G ' + os.path.join(wp.DIAUmpirePath,
                                                           'DIA_Umpire_SE.jar') + ' ' + inputmzXML + ' ' + diaumpireconfigurefile
        print '[Info] >' + diaumpirecmd
        os.system(diaumpirecmd)
        return mgflist


def run_dia_workflow(inputmzXML, wfc, wp):
    """
    :param wfc: work flow controller
    @type wfc: WorkflowController
    """
    workflowname = 'DIAUmpire -->-- Comet -->-- Xinteract Workflow'
    print '[Info] ' + workflowname + ' is running...'

    diaumpire = DIAUmpire(wfc, wp)
    MGFList = diaumpire.Run(wp, inputmzXML)
    # self.RunDIAUmpire()
    cometse = comet(wfc, wp)
    fastapath = os.path.split(inputmzXML)[0]
    db = DB(wfc)
    pepxmllist = cometse.Run(wp, MGFList, db.Run(fastapath, wp))

    xint = Xinteract(wfc, wp)
    xint.Run(pepxmllist)
    # xint.RunAll(pepxmllist)
    # RunXinteract()
    print '[Info] DIAUmpire -->-- Comet -->-- Xinteract Workflow finished'


if __name__ == '__main__':
    print 'Test WorkFlow'
    mzXML = r'H:\lwu\testdata\t1.mzXML'
    # wp = PipelineConfig()
    code_path = os.path.split(sys.argv[0])[0]
    wp = PipelineConfig()
    wp.load_from_config_file(os.path.join(code_path, 'RTQuant.ini'))
    wfc = WorkflowController()
    run_dia_workflow(mzXML, wfc, wp)
