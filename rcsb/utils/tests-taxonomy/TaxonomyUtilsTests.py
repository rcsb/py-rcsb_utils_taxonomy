# File:    TaxonomyUtilsTests.py
# Author:  J. Westbrook
# Date:    9-Mar-2019
# Version: 0.001
#
# Update:
#   25-Mar-2019  jdw add test for merged taxons
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
        self.__useCache = True
        #
        self.__workPath = os.path.join(HERE, 'test-output')
        if not self.__useCache:
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
            psn = tU.getParentScientificName(taxId)
            logger.info("Parent scientific name %s" % psn)
        except Exception as e:
            logger.exception("Failing with %s" % str(e))
            self.fail()

    def testLineageTaxonomyData(self):
        """ Test the case read taxonomy resource file
        """
        try:
            taxId = 1354003
            tU = TaxonomyUtils(taxDirPath=self.__workPath)
            taxId = tU.getMergedTaxId(taxId)
            sn = tU.getScientificName(taxId)
            logger.debug("Scientific name (%d): %s" % (taxId, sn))

            tL = tU.getLineage(taxId)
            logger.info("tL(%d) %r" % (len(tL), tL))
            #
            tL = tU.getLineageWithNames(taxId)
            logger.info("tL(%d) %r" % (len(tL), tL))
            self.assertGreaterEqual(len(tL), 32)
            psn = tU.getParentScientificName(taxId)
            logger.info("Parent scientific name %s" % psn)
        except Exception as e:
            logger.exception("Failing with %s" % str(e))
            self.fail()

    def testMissingTaxIds(self):
        """
        """
        missingTaxIds = [100641, 10559, 10624, 10625, 10626, 10629,
                         10633, 10634, 107412, 1077152, 10870, 10871,
                         10930, 10969, 110206, 11161, 111792, 11207,
                         11498, 11531, 11593, 11598, 11599, 11604,
                         11613, 11619, 11630, 1172, 11741, 126825,
                         129818, 135726, 136134, 145481, 1520516,
                         152465, 160970, 1620918, 1637813, 1826778,
                         184226, 195041, 195949, 197185, 213401,
                         216978, 2244, 2274, 241811, 255776, 260809,
                         265180, 269483, 27317, 2775, 28451,
                         28575, 289331, 29163, 29165, 295351, 3053,
                         311449, 314269, 316397, 3263, 33166, 34362,
                         35758, 36924, 373547, 37705, 378181, 380086,
                         412977, 420151, 423445, 434928, 45596, 45996,
                         46607, 469587, 49008, 4905, 4906, 4963, 5085,
                         535911, 5386, 53957, 5457, 5532, 558164,
                         563191, 5724, 591, 592, 5937, 601, 602,
                         60708, 634178, 6352, 6419, 6437, 654421,
                         66077, 6621, 6633, 6634, 6644, 668994, 67004,
                         69803, 6988, 73900, 7404, 74103, 74577,
                         79679, 82835, 8836, 887, 89769, 92694, 92829,
                         92830, 95649, 96564, 96635, 97708,
                         999953]
        try:
            tU = TaxonomyUtils(taxDirPath=self.__workPath)
            count = 0
            for taxId in missingTaxIds:
                if tU.getMergedTaxId(taxId) == taxId:
                    logger.info("Missing taxid %d" % taxId)
                    count += 1
            self.assertEqual(count, 0)
        except Exception as e:
            logger.exception("Failing with %s" % str(e))
            self.fail()

    def testExportTree(self):
        """ Test export taxonomy data in a particular data structure.
        """
        try:
            useCache = False
            if not useCache:
                cachePath = os.path.join(self.__workPath, 'NCBI')
                if not os.access(cachePath, os.W_OK):
                    try:
                        os.makedirs(cachePath, 0o755)
                    except Exception:
                        pass
                ##
                fpL = glob.glob(os.path.join(cachePath, "*"))
                if fpL:
                    for fp in fpL:
                        os.remove(fp)
            dL = []
            tU = TaxonomyUtils(taxDirPath=cachePath, useCache=useCache)
            dL = tU.exportNodeList()
            logger.info("Node list length %d" % len(dL))
            logger.info("Node list %r" % dL[:20])
        except Exception as e:
            logger.exception("Failing with %s" % str(e))
            self.fail()


def utilReadSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(TaxonomyUtilsTests("testAccessTaxonomyData"))
    suiteSelect.addTest(TaxonomyUtilsTests("testLineageTaxonomyData"))
    suiteSelect.addTest(TaxonomyUtilsTests("testMissingTaxIds"))
    return suiteSelect


def utilTreeSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(TaxonomyUtilsTests("testExportTree"))
    return suiteSelect


if __name__ == '__main__':
    #
    if True:
        mySuite = utilReadSuite()
        unittest.TextTestRunner(verbosity=2).run(mySuite)

    if True:
        mySuite = utilTreeSuite()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
