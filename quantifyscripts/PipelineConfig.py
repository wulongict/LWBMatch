#!/usr/bin/python
# Python 2.7 Required
# Configuration of the entire workflow; keeps paths to different softwares.
# This script can not run alone.
#
#
# -----------------------
# Author: L. WU
# Date: Aug. 18, 2016
# ------------------------


import os

import platform
import ConfigParser


class PipelineConfig:
    def __init__(self):
        self.systemname = platform.system()
        self.search_parent_ion_tolerance_Th_low = 3.0  # or 0.3
        self.search_fragmentation_ion_tolerance_Th_low = 1.004
        # if self.systemname == 'Linux':
        #     self.binpath = r'/data/wulong/bitbucket/rt-alignment/'
        #     self.lwbmatchpath = r'/data/wulong/bitbucket/rt-alignment/Release'
        #     self.OpenMSpath = r'/usr/local/bin/'
        #     self.mkdecoyDBperl = os.path.join(r'/tools/bin/', 'createDecoyDatabase-interleave.pl')
        #     self.DIAUmpirePath = r'/data/wulong/code/DIAUmpire/'
        #     self.TPPpath = r'/tools/tpp-4.6.3/bin'
        #     self.msconvertpath = self.TPPpath
        #     self.cometpath = r'/tools/bin/'
        #     self.extname = ''
        #     self.recallprecision_binpath = os.path.join(self.binpath, 'Release')
        #     self.cometname = 'comet'
        #
        #
        #
        # else:
        #     self.binpath = r'F:\SpecBackUp\wulongspec\bin'
        #     self.lwbmatchpath = self.binpath
        #     self.OpenMSpath = r'C:\Program Files\OpenMS-1.11\bin'
        #     self.mkdecoyDBperl = os.path.join(self.binpath, 'createDecoyDatabase-interleave.pl')
        #     self.DIAUmpirePath = r'F:\SpecBackUp\wulongspec\DIAUmpire\DIA-Umpire_v1_284'
        #     self.TPPpath = r'C:\Inetpub\tpp-bin'
        #     self.msconvertpath = r'C:\Program Files\ProteoWizard\ProteoWizard 3.0.7389'
        #     self.extname = '.exe'
        #     self.cometpath = self.binpath
        #     self.recallprecision_binpath = self.binpath
        #     self.cometname = 'comet'

    def Print(self):
        print self.__dict__

    def output_as_configure_file(self, outputfile):

        config = ConfigParser.RawConfigParser()

        config.add_section('basic settings')
        config.set('basic settings', 'TPPpath', self.TPPpath)
        config.set('basic settings', 'msconvertpath', self.msconvertpath)
        config.set('basic settings', 'extname', self.extname)
        config.set('basic settings', 'cometpath', self.cometpath)
        config.set('basic settings', 'recallprecision_binpath', self.binpath)
        config.set('basic settings', 'binpath', self.binpath)
        config.set('basic settings', 'lwbmatchpath', self.lwbmatchpath)
        config.set('basic settings', 'OpenMSpath', self.OpenMSpath)
        config.set('basic settings', 'mkdecoyDBperl', self.mkdecoyDBperl)
        config.set('basic settings', 'DIAUmpirepath', self.DIAUmpirePath)
        config.set('basic settings', 'systemname', self.systemname)
        config.set('basic settings', 'cometname', self.cometname)

        print 'start output -------------', outputfile

        with open(outputfile, 'wb') as configfile:
            config.write(configfile)

    def load_from_config_file(self, configfile):
        '''
        :param configfile:
        :return:

        config = ConfigParser.RawConfigParser()
        config.read('example.cfg')
        '''

        config = ConfigParser.RawConfigParser()
        config.read(configfile)

        print 'TPPpath is ', config.get('basic settings', 'TPPpath')

        # Set the params into the variables

        self.binpath = config.get('basic settings', 'binpath')
        self.lwbmatchpath = config.get('basic settings', 'lwbmatchpath')
        self.OpenMSpath = config.get('basic settings', 'OpenMSpath')
        self.mkdecoyDBperl = config.get('basic settings', 'mkdecoyDBperl')
        self.DIAUmpirePath = config.get('basic settings', 'DIAUmpirePath')
        self.TPPpath = config.get('basic settings', 'TPPpath')
        self.msconvertpath = config.get('basic settings', 'msconvertpath')
        self.cometpath = config.get('basic settings', 'cometpath')
        self.extname = config.get('basic settings', 'extname')
        self.recallprecision_binpath = config.get('basic settings', 'recallprecision_binpath')
        self.cometname = config.get('basic settings', 'cometname')

        self.configs = config
        print '[Info] Parameters are as following: '
        print self.__dict__

        return config
