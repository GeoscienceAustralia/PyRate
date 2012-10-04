'''
Created on 17/09/2012
@author: bpd900
'''

import os, unittest

import main



class PyRateTests(unittest.TestCase):

	def test_pyrate_main(self):
		cwd = os.getcwd()
		td = "../../tests/sydney_test"

		os.chdir(td)

		main.main(verbose=True)

		os.chdir(cwd)


if __name__ == "__main__":
	unittest.main()
