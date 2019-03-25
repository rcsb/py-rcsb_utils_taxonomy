##
# File: TaxonomyUtils.py
# Author:  J. Westbrook
# Date:    9-Mar-2019
# Version: 0.001
#
# Updates:
# 23-Mar-2019 jdw make cache file names python version specific
# 24-Mar-2019 jdw add leaf node to taxonomy lineage
# 25-Mar-2019 jdw add tests for merged taxons and method getMergedTaxId()
##

import logging
import os.path
import sys

from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class TaxonomyUtils(object):

    def __init__(self, **kwargs):
        """
        """
        self.__taxDirPath = kwargs.get("taxDirPath", '.')
        useCache = kwargs.get("useCache", True)
        #
        self.__urlTarget = kwargs.get("ncbiTaxonomyUrl", "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz")
        #
        self.__mU = MarshalUtil(workPath=self.__taxDirPath)
        #
        self.__nodeD = {}
        self.__nameD = {}
        self.__mergeD = {}
        #
        self.__nameD, self.__nodeD, self.__mergeD = self.__reload(self.__urlTarget, self.__taxDirPath, useCache=useCache)

    def getMergedTaxId(self, taxId):
        try:
            taxId = self.__mergeD[taxId] if taxId in self.__mergeD else taxId
        except Exception:
            pass
        return taxId

    def getScientificName(self, taxId):
        try:
            taxId = self.__mergeD[taxId] if taxId in self.__mergeD else taxId
            return self.__nameD[int(taxId)]['sn']
        except Exception:
            pass
        return None

    def getCommonNames(self, taxId):
        try:
            taxId = self.__mergeD[taxId] if taxId in self.__mergeD else taxId
            return list(set(self.__nameD[int(taxId)]['cn']))
        except Exception:
            pass
        return None

    def getParentTaxid(self, taxId):
        try:
            taxId = self.__mergeD[taxId] if taxId in self.__mergeD else taxId
            return self.__nodeD[int(taxId)][0]
        except Exception:
            pass
        return None

    def getLineage(self, taxId):
        pList = []
        try:
            taxId = self.__mergeD[taxId] if taxId in self.__mergeD else taxId
            pList.append(taxId)
            pt = self.getParentTaxid(taxId)
            while ((pt is not None) and (pt != 1)):
                pList.append(pt)
                pt = self.getParentTaxid(pt)
        except Exception as e:
            logger.exception("Failing with %s" % str(e))
        #
        pList.reverse()
        return pList

    def getLineageWithNames(self, taxId):
        rL = []
        try:
            taxId = self.__mergeD[taxId] if taxId in self.__mergeD else taxId
            pTaxIdL = self.getLineage(taxId)
            for ii, pTaxId in enumerate(pTaxIdL, 1):
                nmL = [self.getScientificName(pTaxId)]
                cnL = self.getCommonNames(pTaxId)
                if cnL:
                    nmL.extend(cnL)
                nmL = set(nmL)
                for nm in nmL:
                    rL.append((ii, pTaxId, nm))
        except Exception as e:
            logger.exception("Failing with %s" % str(e))
        #
        return rL

    def getParentList(self, taxId):
        """ Return a list of tuples containing taxid & scientific name for
            parents of the input taxid.
        """
        taxId = self.__mergeD[taxId] if taxId in self.__mergeD else taxId
        pList = []
        pt = self.getParentTaxid(taxId)
        while ((pt is not None) and (pt != '1')):
            pList.append((pt, self.getScientificName(pt)))
            pt = self.getParentTaxid(pt)
        return pList
    #
    ##

    def getLineageScientificNames(self, taxId):
        nmL = []
        try:
            taxId = self.__mergeD[taxId] if taxId in self.__mergeD else taxId
            lineage = self.getLineage(taxId)
            nmL = [self.getScientificName(taxid) for taxid in lineage]
        except Exception as e:
            logger.exception("Failing for taxId %r with %s" % (taxId, str(e)))
        return nmL

    def __getLineageD(self, taxId):
        pD = {}
        try:
            taxId = self.__mergeD[taxId] if taxId in self.__mergeD else taxId
            pD[taxId] = True
            pt = self.getParentTaxid(taxId)
            while ((pt is not None) and (pt != 1)):
                pD[pt] = True
                pt = self.getParentTaxid(pt)
        except Exception as e:
            logger.exception("Failing with %s" % str(e))
        return pD

    def isBacteria(self, taxId):
        try:
            lineage = self.__getLineageD(taxId)
            return True if 2 in lineage else False
        except Exception as e:
            logger.exception("Failing for taxId %r with %s" % (taxId, str(e)))
        return False

    def isEukaryota(self, taxId):
        try:
            taxId = self.__mergeD[taxId] if taxId in self.__mergeD else taxId
            lineage = self.__getLineageD(taxId)
            return True if 2759 in lineage else False
        except Exception as e:
            logger.exception("Failing for taxId %r with %s" % (taxId, str(e)))
        return False

    def isVirus(self, taxId):
        try:
            taxId = self.__mergeD[taxId] if taxId in self.__mergeD else taxId
            lineage = self.__getLineageD(taxId)
            return True if 10239 in lineage else False
        except Exception as e:
            logger.exception("Failing for taxId %r with %s" % (taxId, str(e)))
        return False

    def isArchaea(self, taxId):
        try:
            taxId = self.__mergeD[taxId] if taxId in self.__mergeD else taxId
            lineage = self.__getLineageD(taxId)
            return True if 2157 in lineage else False
        except Exception as e:
            logger.exception("Failing for taxId %r with %s" % (taxId, str(e)))
        return False

    def isOther(self, taxId):
        """ other/synthetic
        """
        try:
            taxId = self.__mergeD[taxId] if taxId in self.__mergeD else taxId
            lineage = self.__getLineageD(taxId)
            return True if 28384 in lineage else False
        except Exception as e:
            logger.exception("Failing for taxId %r with %s" % (taxId, str(e)))
        return False

    def isUnclassified(self, taxId):
        try:
            taxId = self.__mergeD[taxId] if taxId in self.__mergeD else taxId
            lineage = self.__getLineageD(taxId)
            return True if 12908 in lineage else False
        except Exception as e:
            logger.exception("Failing for taxId %r with %s" % (taxId, str(e)))
        return False

    #
    def __reload(self, urlTarget, taxDirPath, useCache=True):
        pyVersion = sys.version_info[0]
        taxNamePath = os.path.join(taxDirPath, "taxonomy_names-py%s.pic" % str(pyVersion))
        taxNodePath = os.path.join(taxDirPath, "taxonomy_nodes-py%s.pic" % str(pyVersion))
        taxMergedNodePath = os.path.join(taxDirPath, "taxonomy_nodes-merged-py%s.pic" % str(pyVersion))
        #
        if useCache and self.__mU.exists(taxNamePath) and self.__mU.exists(taxNamePath):
            tD = self.__mU.doImport(taxNamePath, format="pickle")
            nD = self.__mU.doImport(taxNodePath, format="pickle")
            mD = self.__mU.doImport(taxMergedNodePath, format="pickle")
            logger.debug("Taxonomy name length %d node length %d" % (len(tD), len(nD)))
        else:
            logger.info("Fetch taxonomy data from source %s" % urlTarget)
            nmL, ndL, mergeL = self.__fetchFromSource(urlTarget, taxDirPath)
            tD = self.__extractNames(nmL)
            ok = self.__mU.doExport(taxNamePath, tD, format="pickle")
            #
            nD = self.__extractNodes(ndL)
            ok = self.__mU.doExport(taxNodePath, nD, format="pickle") and ok
            #
            mD = self.__mergedTaxids(mergeL)
            ok = self.__mU.doExport(taxMergedNodePath, mD, format="pickle") and ok
        #
        return tD, nD, mD

    def __mergedTaxids(self, rowL):
        """ Extract taxonomy names and synonyms from NCBI taxonomy database dump file row list.
        """
        tD = {}
        try:
            for t in rowL:
                if len(t) < 2:
                    continue
                taxId = int(t[0])
                mergedTaxId = int(t[2])
                tD[taxId] = mergedTaxId

            logger.debug("Taxon merged dictionary length %d \n" % len(tD))
            #
        except Exception as e:
            logger.exception("Failing with %s" % str(e))
        return tD

    def __extractNames(self, rowL):
        """ Extract taxonomy names and synonyms from NCBI taxonomy database dump file row list.
        """
        tD = {}
        try:
            # csvL = []
            for t in rowL:
                if len(t) < 7:
                    continue
                taxId = int(t[0])
                name = t[2]
                nameType = t[6]
                # csvL.append({'t': taxId, 'name': name, 'type': nameType})
                #
                if nameType in ['scientific name', 'common name', 'synonym', 'genbank common name']:
                    if taxId not in tD:
                        tD[taxId] = {}
                    if nameType in ['scientific name']:
                        tD[taxId]['sn'] = name
                        continue
                    if 'cn' not in tD[taxId]:
                        tD[taxId]['cn'] = []
                    tD[taxId]['cn'].append(name)
                else:
                    pass
            logger.debug("Taxonomy dictionary length %d \n" % len(tD))
            #
        except Exception as e:
            logger.exception("Failing with %s" % str(e))
        return tD

    def __extractNodes(self, rowL):
        """ Extract taxonomy parent relationships from NCBI taxonomy database dump row list.
        """
        nD = {}
        try:
            for fields in rowL:
                taxid = int(fields[0].strip())
                parentTaxid = int(fields[2].strip())
                rank = str(fields[4]).strip()
                nD[taxid] = (parentTaxid, rank)
            logger.debug("Taxonomy parent dictionary length %d \n" % len(nD))
        except Exception as e:
            logger.exception("Failing with %s" % str(e))
        return nD

    def __fetchFromSource(self, urlTarget, taxDirPath):
        """  Fetch the ncbi taxonomy dump and extract name and node data.
        """
        pth, fn = os.path.split(urlTarget)
        #
        nmL = self.__mU.doImport(urlTarget, format='tdd', rowFormat='list', tarMember='names.dmp', uncomment=False)
        ndL = self.__mU.doImport(os.path.join(taxDirPath, fn), format='tdd', rowFormat='list', tarMember='nodes.dmp', uncomment=False)
        mergeL = self.__mU.doImport(os.path.join(taxDirPath, fn), format='tdd', rowFormat='list', tarMember='merged.dmp', uncomment=False)
        # deleteL = self.__mU.doImport(os.path.join(taxDirPath, fn), format='tdd', rowFormat='list', tarMember='delnodes.dmp')
        return nmL, ndL, mergeL
