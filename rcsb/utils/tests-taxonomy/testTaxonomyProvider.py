# File:    testTaxonomyProvider.py
# Author:  J. Westbrook
# Date:    9-Mar-2019
# Version: 0.001
#
# Update:
#   25-Mar-2019  jdw add test for merged taxons
#    3-Jan-2022  dwp remove taxID 255776 from list of missingTaxIds, which is now present in the merged.dmp and nodes.dmp NCBI source data files
#                    (from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz);
#                    Also reduce minimum required list lengths in testLineageTaxonomySpecial
#   30-Jan-2023  dwp reduce taxonomy lineage data length threshold from 32 to 30
#    9-Dec-2024  dwp Update test to account for "kingdom" classification level, after having been added to lineage by NCBI (e.g., Pseudomonadati)
#                    (See browser at: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Root)
#   15-Jul-2025  dwp reduce testLineageTaxonomySpecial threshold from 30 to 27 (returned count dropped to 28 around 15-Jun-2025)
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

from importlib.metadata import version as get_package_version
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.taxonomy.TaxonomyProvider import TaxonomyProvider

__version__ = get_package_version("rcsb.utils.taxonomy")

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))


logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)
# logger.setLevel(logging.DEBUG)


class TaxonomyProviderTests(unittest.TestCase):
    def setUp(self):
        self.__verbose = True
        self.__useCache = True
        #
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__cachePathFb = os.path.join(HERE, "test-output", "CACHE-fallback")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__missingNamePath = os.path.join(self.__dataPath, "missingSrcNames.json")
        #
        if not self.__useCache:
            fpL = glob.glob(os.path.join(self.__cachePath, "*.pic"))
            if fpL:
                for fp in fpL:
                    os.remove(fp)

        self.__urlTarget = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
        #
        self.__startTime = time.time()
        logger.debug("Running tests on version %s", __version__)
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testBuildCache(self):
        """Test building the cache from NCBI resource"""
        tU = TaxonomyProvider(cachePath=self.__cachePath, useCache=False)
        ok = tU.testCache()
        self.assertTrue(ok)

    def testAANameLookup(self):
        #
        tU = TaxonomyProvider(cachePath=self.__cachePath, useCache=True)
        ok = tU.testCache()
        #
        if not ok:
            tU = TaxonomyProvider(cachePath=self.__cachePath, useCache=False)
            ok = tU.testCache()
            self.assertTrue(ok)
        #
        mU = MarshalUtil(workPath=self.__cachePath)
        nD = mU.doImport(self.__missingNamePath, fmt="json")
        for nm in nD:
            taxId = tU.getTaxId(nm)
            if not taxId:
                logger.debug("Unknown source name %s (%r)", nm, nD[nm])

    def testFallback(self):
        """Test taxonomy data fallback."""
        try:
            tU = TaxonomyProvider(cachePath=self.__cachePathFb, ncbiTaxonomyUrl=self.__urlTarget, useCache=False)
            ok = tU.testCache()
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testAccessTaxonomyData(self):
        """Test the case read and write taxonomy resource file"""
        try:
            taxId = 9606
            tU = TaxonomyProvider(cachePath=self.__cachePath, useCache=True)

            tL = tU.getLineage(taxId)
            logger.debug("Lineage for %d (%d): %r", taxId, len(tL), tL)
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
            logger.debug("Scientific name (%d): %s", taxId, sn)
            self.assertGreater(len(sn), 10)
            cnL = tU.getCommonNames(taxId)
            self.assertGreaterEqual(len(cnL), 1)
            logger.debug("Common names (%d): %r", taxId, cnL)
            psn = tU.getParentScientificName(taxId)
            logger.debug("Parent scientific name %s", psn)
            #
            taxId = 37486
            sn = tU.getScientificName(taxId)
            logger.debug("Scientific name (%d): %s", taxId, sn)
            self.assertGreater(len(sn), 10)
            cnL = tU.getCommonNames(taxId)
            self.assertGreaterEqual(len(cnL), 1)
            logger.debug("Common names (%d): %r", taxId, cnL)
            psn = tU.getParentScientificName(taxId)
            logger.debug("Parent scientific name %s", psn)
            logger.debug("Lineage names %r", tU.getLineageWithNames(taxId))
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testLineageTaxonomyData(self):
        """Test the case taxonomy lineage"""
        try:
            # Escherichia coli #1/H766
            taxId = 1354003
            tU = TaxonomyProvider(cachePath=self.__cachePath, useCache=True)
            taxId = tU.getMergedTaxId(taxId)
            sn = tU.getScientificName(taxId)
            logger.debug("Scientific name (%d): %s", taxId, sn)

            tL = tU.getLineage(taxId)
            logger.debug("tL(%d) %r", len(tL), tL)
            #
            tL = tU.getLineageWithNames(taxId)
            logger.debug("tL(%d) %r", len(tL), tL)
            self.assertGreaterEqual(len(tL), 30)
            psn = tU.getParentScientificName(taxId)
            logger.debug("Parent scientific name %s", psn)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testLineageTaxonomySpecial(self):
        """Test the case special virus"""
        try:
            # Severe acute respiratory syndrome coronavirus 2
            taxId = 2697049
            tU = TaxonomyProvider(cachePath=self.__cachePath, useCache=True)
            taxId = tU.getMergedTaxId(taxId)
            sn = tU.getScientificName(taxId)
            logger.info("Scientific name (%d): %s", taxId, sn)

            tL = tU.getLineage(taxId)
            logger.info("tL(%d) %r", len(tL), tL)
            #
            tL = tU.getLineageWithNames(taxId)
            logger.debug("tL(%d) %r", len(tL), tL)
            self.assertGreaterEqual(len(tL), 27)
            psn = tU.getParentScientificName(taxId)
            logger.debug("Parent scientific name %s", psn)
            cnL = tU.getCommonNames(taxId)
            self.assertGreaterEqual(len(cnL), 8)
            logger.debug("Common names (%d): %r", taxId, cnL)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMissingTaxIds(self):
        """ """
        missingTaxIds = [
            100641,
            10559,
            10624,
            10625,
            10626,
            10629,
            10633,
            10634,
            107412,
            1077152,
            10870,
            10871,
            10930,
            10969,
            110206,
            11161,
            111792,
            11207,
            11498,
            11531,
            11593,
            11598,
            11599,
            11604,
            11613,
            11619,
            11630,
            1172,
            11741,
            126825,
            129818,
            135726,
            136134,
            145481,
            1520516,
            152465,
            160970,
            1620918,
            1637813,
            1826778,
            184226,
            195041,
            195949,
            197185,
            213401,
            216978,
            2244,
            2274,
            241811,
            260809,
            265180,
            269483,
            2775,
            28451,
            28575,
            289331,
            29163,
            29165,
            295351,
            3053,
            311449,
            314269,
            316397,
            3263,
            33166,
            34362,
            35758,
            36924,
            373547,
            37705,
            378181,
            380086,
            412977,
            420151,
            423445,
            434928,
            45596,
            45996,
            46607,
            469587,
            49008,
            4905,
            4906,
            4963,
            5085,
            535911,
            5386,
            53957,
            5457,
            5532,
            558164,
            563191,
            5724,
            591,
            592,
            5937,
            601,
            602,
            60708,
            634178,
            6352,
            6419,
            6437,
            654421,
            66077,
            6621,
            6633,
            6634,
            6644,
            668994,
            67004,
            69803,
            6988,
            73900,
            7404,
            74103,
            74577,
            79679,
            82835,
            8836,
            887,
            89769,
            92694,
            92829,
            92830,
            95649,
            96564,
            96635,
            97708,
            999953,
        ]
        try:
            tU = TaxonomyProvider(cachePath=self.__cachePath, useCache=True)
            count = 0
            for taxId in missingTaxIds:
                if tU.getMergedTaxId(taxId) == taxId:
                    logger.debug("Missing taxid %d", taxId)
                    count += 1
            self.assertEqual(count, 0)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testExportTree(self):
        """Test export taxonomy data in a particular data structure."""
        try:
            dL = []
            tU = TaxonomyProvider(cachePath=self.__cachePath, useCache=True)
            dL = tU.exportNodeList()
            logger.debug("Node list length %d", len(dL))
            logger.debug("Node list %r", dL[:20])
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testGraphOps(self):
        """Test graph operations."""
        try:
            tU = TaxonomyProvider(cachePath=self.__cachePath, useCache=True)
            #
            startTime = time.time()
            lcTaxId = tU.getLowestCommonAncestor(63221, 741158)
            logger.info("LCA %r in (%.4f seconds)", lcTaxId, time.time() - startTime)
            self.assertEqual(9606, lcTaxId)
            #
            #
            startTime = time.time()
            lcTaxId = tU.getLowestCommonAncestor(9606, 9606)
            logger.info("LCA %r in (%.4f seconds)", lcTaxId, time.time() - startTime)
            self.assertEqual(9606, lcTaxId)
            #
            startTime = time.time()
            lcTaxId = tU.getLowestCommonAncestor(9606, 741158)
            logger.info("LCA %r in (%.4f seconds)", lcTaxId, time.time() - startTime)
            self.assertEqual(9606, lcTaxId)
            #
            startTime = time.time()
            lcTaxId = tU.getLowestCommonAncestor(866768, 2569093)
            logger.info("LCA %r in (%.4f seconds)", lcTaxId, time.time() - startTime)
            self.assertEqual(3379134, lcTaxId)
            #
            startTime = time.time()
            lcTaxId = tU.getLowestCommonAncestor(866768, 91061)
            logger.info("LCA %r in (%.4f seconds)", lcTaxId, time.time() - startTime)
            self.assertEqual(2, lcTaxId)
            #
            startTime = time.time()
            txL = [(63221, 741158), (9606, 9606), (9606, 741158), (866768, 2569093)]
            lcTaxIdPairD = tU.getLowestCommonAncestors(txL)
            logger.info("List LCA %r in (%.4f seconds)", lcTaxIdPairD, time.time() - startTime)
            #
            for taxPair in txL:
                startTime = time.time()
                status, lcTaxId, rank = tU.compareTaxons(taxPair[0], taxPair[1])
                logger.info("status %r lca %r rank %r in (%.4f seconds)", status, lcTaxId, rank, time.time() - startTime)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def utilReadSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(TaxonomyProviderTests("testBuildCache"))
    suiteSelect.addTest(TaxonomyProviderTests("testAANameLookup"))
    suiteSelect.addTest(TaxonomyProviderTests("testFallback"))
    suiteSelect.addTest(TaxonomyProviderTests("testAccessTaxonomyData"))
    suiteSelect.addTest(TaxonomyProviderTests("testLineageTaxonomyData"))
    suiteSelect.addTest(TaxonomyProviderTests("testLineageTaxonomySpecial"))
    suiteSelect.addTest(TaxonomyProviderTests("testMissingTaxIds"))
    return suiteSelect


def utilTreeSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(TaxonomyProviderTests("testExportTree"))
    suiteSelect.addTest(TaxonomyProviderTests("testGraphOps"))
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite = utilReadSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)

    mySuite = utilTreeSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
