# File:    TaxonomyUtilsTests.py
# Author:  J. Westbrook
# Date:    9-Mar-2019
# Version: 0.001
#
# Update:
#
##
"""
Tests for extraction, supplementing and packaging dictionary metadata for schema construction.

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import glob
import logging
import os
import time
import unittest

from rcsb.utils.io import __version__
from rcsb.utils.taxonomy.TaxonomyUtils import TaxonomyUtils

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))


logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s')
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class TaxonomyUtilsTests(unittest.TestCase):

    def setUp(self):
        self.__verbose = True
        #
        self.__workPath = os.path.join(HERE, 'test-output')
        fpL = glob.glob(os.path.join(self.__workPath, "*.pic"))
        if fpL:
            for fp in fpL:
                os.remove(fp)

        self.__urlTarget = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
        #
        self.__startTime = time.time()
        logger.debug("Running tests on version %s" % __version__)
        logger.debug("Starting %s at %s" % (self.id(),
                                            time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

    def tearDown(self):
        endTime = time.time()
        logger.debug("Completed %s at %s (%.4f seconds)" % (self.id(),
                                                            time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                            endTime - self.__startTime))

    def testAccessTaxonomyData(self):
        """ Test the case read and write taxonomy resource file
        """
        try:
            taxId = 9606
            tU = TaxonomyUtils(taxDirPath=self.__workPath)
            tL = tU.getLineage(taxId)
            logger.debug("Lineage for %d (%d): %r" % (taxId, len(tL), tL))
            self.assertGreaterEqual(len(tL), 30)
            #
            ok = tU.isEukaryota(taxId)
            self.assertTrue(ok)

            ok = tU.isVirus(taxId)
            self.assertFalse(ok)

            ok = tU.isArchaea(taxId)
            self.assertFalse(ok)
            ok = tU.isOther(taxId)
            self.assertFalse(ok)
            ok = tU.isBacteria(taxId)
            self.assertFalse(ok)
            ok = tU.isUnclassified(taxId)
            self.assertFalse(ok)
            #
            sn = tU.getScientificName(taxId)
            logger.debug("Scientific name (%d): %s" % (taxId, sn))
            self.assertGreater(len(sn), 10)
            cnL = tU.getCommonNames(taxId)
            self.assertGreaterEqual(len(cnL), 2)
            logger.debug("Common names (%d): %r" % (taxId, cnL))

        except Exception as e:
            logger.exception("Failing with %s" % str(e))
            self.fail()

    def testLineageTaxonomyData(self):
        """ Test the case read taxonomy resource file
        """
        try:
            taxId = 9606
            tU = TaxonomyUtils(taxDirPath=self.__workPath)
            tL = tU.getLineageWithNames(taxId)
            logger.debug("tL(%d) %r" % (len(tL), tL))
            self.assertGreaterEqual(len(tL), 60)
        except Exception as e:
            logger.exception("Failing with %s" % str(e))
            self.fail()


def utilReadSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(TaxonomyUtilsTests("testAccessTaxonomyData"))
    suiteSelect.addTest(TaxonomyUtilsTests("testLineageTaxonomyData"))
    return suiteSelect


if __name__ == '__main__':
    #
    if True:
        mySuite = utilReadSuite()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
