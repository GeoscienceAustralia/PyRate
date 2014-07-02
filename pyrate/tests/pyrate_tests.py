'''
Tests end-to-end runs of the PyRate workflow module.

Created on 17/09/2012
@author: Ben Davies, NCI
'''

import os
from os.path import join
import shutil

import unittest
from pyrate import pyrate


# testing constants
BASE_DIR = '/tmp/pyrate/workflow'
BASE_TIF_DIR = join(BASE_DIR, 'tif')
BASE_CFG_FILE = join(BASE_DIR, 'pyrate_workflow_test.conf')

TEST_CORE = join(os.environ['PYRATEPATH'], 'tests', 'sydney_test')

CURRENT_DIR = os.getcwd()



class PyRateTests(unittest.TestCase):
	# initialise and run workflow from the class setup
	# unit tests verify the different steps have completed

	@classmethod
	def setUpClass(cls):
		try:
			# assume outputs/working dir are part of same tree
			# link the sydney data in
			if not os.path.exists(BASE_DIR):
				os.makedirs(BASE_DIR)

			# link to config file and source data
			if not os.path.exists(BASE_CFG_FILE):
				orig_cfg = join(TEST_CORE, 'pyrate_workflow_test.conf')
				os.symlink(orig_cfg, BASE_CFG_FILE)

			if not os.path.exists(BASE_TIF_DIR):
				orig_tif = join(TEST_CORE, 'tif')
				os.symlink(orig_tif, BASE_TIF_DIR)

			os.chdir(BASE_DIR)
			pyrate.main(BASE_CFG_FILE, verbose=False)
		except:
			# get out to avoid paths busting other tests
			os.chdir(CURRENT_DIR)
			raise

	@classmethod
	def tearDownClass(cls):
		shutil.rmtree(BASE_DIR)


	def test_initial_setup(self):
		self.assertTrue(os.path.exists(BASE_DIR))
		self.assertTrue(os.path.exists(BASE_TIF_DIR))
		self.assertTrue(os.path.exists(BASE_CFG_FILE))



if __name__ == "__main__":
	unittest.main()

