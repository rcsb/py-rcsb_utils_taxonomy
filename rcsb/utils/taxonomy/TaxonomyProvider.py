##
# File: TaxonomyProvider.py
# Author:  J. Westbrook
# Date:    9-Mar-2019
# Version: 0.001
#
# Updates:
# 23-Mar-2019 jdw make cache file names python version specific
# 24-Mar-2019 jdw add leaf node to taxonomy lineage
# 25-Mar-2019 jdw add tests for merged taxons and method getMergedTaxId()
#  1-Apr-2019 jdw add method getChildren() for adjacent children.
# 24-Apr-2019 jdw Add filter option to node list generator and exclude synthetic root from exported tree.
# 14-May-2019 jdw cast all access to mergeD[]
# 23-Jul-2019 jdw adjustments to preserve ordering.
##

import collections
import logging
import os.path
import sys

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class TaxonomyProvider(object):
    def __init__(self, **kwargs):
        """"""
        self.__taxDirPath = os.path.abspath(kwargs.get("taxDirPath", "."))
        useCache = kwargs.get("useCache", True)
        #
        self.__urlTarget = kwargs.get("ncbiTaxonomyUrl", "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz")
        #
        self.__mU = MarshalUtil(workPath=self.__taxDirPath)
        #
        self.__nodeD = {}
        self.__nameD = {}
        self.__mergeD = {}
        self.__childD = {}
        self.__taxIdToNameD = {}
        #
        self.__nameD, self.__nodeD, self.__mergeD = self.__reload(self.__urlTarget, self.__taxDirPath, useCache=useCache)

    def testCache(self):
        # Lengths name 2133961 node 2133961 merge 54768
        logger.debug("Lengths name %d node %d merge %d", len(self.__nameD), len(self.__nodeD), len(self.__mergeD))
        if (len(self.__nameD) > 2100000) and (len(self.__nodeD) > 2100000) and (len(self.__mergeD) > 54000):
            return True
        return False

    def getTaxId(self, organismName):
        if not self.__taxIdToNameD:
            self.__buildTaxIdNameMap()
        try:
            ret = self.__taxIdToNameD[organismName.strip().upper()]
        except Exception:
            ret = None
        return ret

    def __buildTaxIdNameMap(self):
        tD = {}
        for taxId, nmD in self.__nameD.items():
            try:
                tD[nmD["sn"].strip().upper()] = taxId
                tD[nmD["alt"].strip().upper()] = taxId
                for nm in sorted(set(self.__nameD[int(taxId)]["cn"])):
                    tD[nm.strip().upper()] = taxId
            except Exception:
                pass
        self.__taxIdToNameD = tD

    def getMergedTaxId(self, taxId):
        try:
            taxId = self.__mergeD[int(taxId)] if isinstance(taxId, int) and int(taxId) in self.__mergeD else taxId
        except Exception:
            pass
        return taxId

    def getScientificName(self, taxId):
        try:
            taxId = self.__mergeD[int(taxId)] if isinstance(taxId, int) and int(taxId) in self.__mergeD else taxId
            return self.__nameD[int(taxId)]["sn"]
        except Exception:
            pass
        return None

    def getParentScientificName(self, taxId, depth=1):
        """Return the scientific name for the parent of the input taxId at the input lineage depth."""
        try:
            taxId = self.__mergeD[int(taxId)] if isinstance(taxId, int) and int(taxId) in self.__mergeD else taxId
            iL = self.getLineage(taxId)
            tId = iL[depth]
            return self.__nameD[int(tId)]["sn"]
        except Exception:
            pass
        return None

    def getAlternateName(self, taxId):
        """Approximately, the preferred common name."""
        try:
            taxId = self.__mergeD[int(taxId)] if isinstance(taxId, int) and int(taxId) in self.__mergeD else taxId
            return self.__nameD[int(taxId)]["alt"]
        except Exception:
            pass
        return None

    def getCommonNames(self, taxId):
        try:
            taxId = self.__mergeD[int(taxId)] if isinstance(taxId, int) and int(taxId) in self.__mergeD else taxId
            return sorted(set(self.__nameD[int(taxId)]["cn"]))
        except Exception:
            pass
        return None

    def getParentTaxid(self, taxId):
        try:
            taxId = self.__mergeD[int(taxId)] if isinstance(taxId, int) and int(taxId) in self.__mergeD else taxId
            return self.__nodeD[int(taxId)][0]
        except Exception:
            pass
        return None

    def getLineage(self, taxId):
        pList = []
        try:
            taxId = self.__mergeD[int(taxId)] if isinstance(taxId, int) and int(taxId) in self.__mergeD else taxId
            pList.append(taxId)
            pt = self.getParentTaxid(taxId)
            while (pt is not None) and (pt != 1):
                pList.append(pt)
                pt = self.getParentTaxid(pt)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        pList.reverse()
        return pList

    def getLineageWithNames(self, taxId):
        rL = []
        try:
            taxId = self.__mergeD[int(taxId)] if isinstance(taxId, int) and int(taxId) in self.__mergeD else taxId
            pTaxIdL = self.getLineage(taxId)
            for ii, pTaxId in enumerate(pTaxIdL, 1):
                nmL = [self.getScientificName(pTaxId)]
                cnL = self.getCommonNames(pTaxId)
                if cnL:
                    for cn in cnL:
                        if cn in nmL:
                            continue
                        nmL.append(cn)
                for nm in nmL:
                    rL.append((ii, pTaxId, nm))
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return rL

    def getParentList(self, taxId):
        """Return a list of tuples containing taxid & scientific name for
        parents of the input taxid.
        """
        taxId = self.__mergeD[int(taxId)] if isinstance(taxId, int) and int(taxId) in self.__mergeD else taxId
        pList = []
        pt = self.getParentTaxid(taxId)
        while (pt is not None) and (pt != "1"):
            pList.append((pt, self.getScientificName(pt)))
            pt = self.getParentTaxid(pt)
        return pList

    #
    ##

    def getLineageScientificNames(self, taxId):
        nmL = []
        try:
            taxId = self.__mergeD[int(taxId)] if isinstance(taxId, int) and int(taxId) in self.__mergeD else taxId
            lineage = self.getLineage(taxId)
            nmL = [self.getScientificName(taxid) for taxid in lineage]
        except Exception as e:
            logger.exception("Failing for taxId %r with %s", taxId, str(e))
        return nmL

    def __getLineageD(self, taxId):
        pD = {}
        try:
            taxId = self.__mergeD[int(taxId)] if isinstance(taxId, int) and int(taxId) in self.__mergeD else taxId
            pD[taxId] = True
            pt = self.getParentTaxid(taxId)
            while (pt is not None) and (pt != 1):
                pD[pt] = True
                pt = self.getParentTaxid(pt)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return pD

    def isBacteria(self, taxId):
        try:
            lineage = self.__getLineageD(taxId)
            return True if 2 in lineage else False
        except Exception as e:
            logger.exception("Failing for taxId %r with %s", taxId, str(e))
        return False

    def isEukaryota(self, taxId):
        try:
            taxId = self.__mergeD[int(taxId)] if isinstance(taxId, int) and int(taxId) in self.__mergeD else taxId
            lineage = self.__getLineageD(taxId)
            return True if 2759 in lineage else False
        except Exception as e:
            logger.exception("Failing for taxId %r with %s", taxId, str(e))
        return False

    def isVirus(self, taxId):
        try:
            taxId = self.__mergeD[int(taxId)] if isinstance(taxId, int) and int(taxId) in self.__mergeD else taxId
            lineage = self.__getLineageD(taxId)
            return True if 10239 in lineage else False
        except Exception as e:
            logger.exception("Failing for taxId %r with %s", taxId, str(e))
        return False

    def isArchaea(self, taxId):
        try:
            taxId = self.__mergeD[int(taxId)] if isinstance(taxId, int) and int(taxId) in self.__mergeD else taxId
            lineage = self.__getLineageD(taxId)
            return True if 2157 in lineage else False
        except Exception as e:
            logger.exception("Failing for taxId %r with %s", taxId, str(e))
        return False

    def isOther(self, taxId):
        """other/synthetic"""
        try:
            taxId = self.__mergeD[int(taxId)] if isinstance(taxId, int) and int(taxId) in self.__mergeD else taxId
            lineage = self.__getLineageD(taxId)
            return True if 28384 in lineage else False
        except Exception as e:
            logger.exception("Failing for taxId %r with %s", taxId, str(e))
        return False

    def isUnclassified(self, taxId):
        try:
            taxId = self.__mergeD[int(taxId)] if isinstance(taxId, int) and int(taxId) in self.__mergeD else taxId
            lineage = self.__getLineageD(taxId)
            return True if 12908 in lineage else False
        except Exception as e:
            logger.exception("Failing for taxId %r with %s", taxId, str(e))
        return False

    def exportNodeList(self, startTaxId=1, rootTaxId=1, filterD=None):
        """Test export taxonomy data in a particular node list data structure.

        Note: now excluding root node from node list.

        """
        try:
            dL = []
            taxIdList = self.getBfsTraverseList(startTaxId)
            logger.info("Full taxon list length %d", len(taxIdList))
            if filterD:
                taxIdList = [tId for tId in taxIdList if tId in filterD]
            logger.info("Filtered taxon list length %d", len(taxIdList))
            # rootTaxids = [131567, 10239, 28384, 12908]
            for taxId in taxIdList:
                sn = self.getScientificName(taxId)
                altName = self.getAlternateName(taxId)
                displayName = sn + " (" + altName + ")" if altName else sn
                #
                pTaxId = self.getParentTaxid(taxId)
                #
                if sn is None or not sn:
                    logger.info("Unexpected null taxon %r sn %r", taxId, sn)
                #
                # logger.debug("Scientific name (%d): %s (%r)" % (taxId, sn, altName))
                #
                if pTaxId == taxId:
                    lL = []
                else:
                    lL = self.getLineage(taxId)[1:]

                #
                # d = {'id': taxId, 'name': displayName, 'lineage': lL, 'parents': [pTaxId], 'depth': len(lL)}
                # d = {'id': str(taxId), 'name': displayName, 'lineage': [str(t) for t in lL], 'parents': [str(pTaxId)], 'depth': len(lL)}
                #
                if taxId == rootTaxId:
                    continue
                elif pTaxId == rootTaxId:
                    dD = {"id": str(taxId), "name": displayName, "depth": 0}
                else:
                    dD = {"id": str(taxId), "name": displayName, "parents": [str(pTaxId)], "depth": len(lL)}
                dL.append(dD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        taxTreePath = os.path.join(self.__taxDirPath, "taxonomy_node_tree.json")
        self.__mU.doExport(taxTreePath, dL, fmt="json", indent=3)
        return dL

    def getBfsTraverseList(self, startTaxId):
        """Traverse the taxonomy nodes in bfs order starting from startTaxId."""
        #
        tL = []
        visited = set([startTaxId])
        queue = collections.deque(visited)
        while queue:
            taxId = queue.popleft()
            tL.append(taxId)
            for childTaxId in self.getChildren(taxId):
                if childTaxId not in visited:
                    queue.append(childTaxId)
                    visited.add(childTaxId)
        return tL

    def getChildren(self, taxId):
        cL = []
        try:
            self.__childD = self.__getAdjacentDecendants(self.__nodeD) if not self.__childD else self.__childD
            cL = self.__childD[taxId]
        except Exception as e:
            logger.debug("For %r failing with %s", taxId, str(e))
        return cL
        #

    def __getAdjacentDecendants(self, nodeD):
        """Convert d[childTaxId] = (parentTaxid,rank) to d[parentTaxid] = [childTaxid, ... ]"""
        cD = {}
        logger.debug("Parent node lookup dictionary length %d", len(nodeD))
        for childTaxId, (parentTaxId, _) in nodeD.items():
            cD.setdefault(parentTaxId, []).append(childTaxId)
        #
        logger.debug("Adjacent child lookup dictionary length %d", len(cD))
        return cD

    #
    def __reload(self, urlTarget, taxDirPath, useCache=True):
        pyVersion = sys.version_info[0]
        taxNamePath = os.path.join(taxDirPath, "taxonomy_names-py%s.pic" % str(pyVersion))
        taxNodePath = os.path.join(taxDirPath, "taxonomy_nodes-py%s.pic" % str(pyVersion))
        taxMergedNodePath = os.path.join(taxDirPath, "taxonomy_nodes-merged-py%s.pic" % str(pyVersion))
        #
        logger.debug("Using taxonomy data path %s", taxDirPath)
        self.__mU.mkdir(taxDirPath)
        if not useCache:
            for fp in [taxNamePath, taxNodePath, taxMergedNodePath]:
                try:
                    os.remove(fp)
                except Exception:
                    pass
        #
        if useCache and self.__mU.exists(taxNamePath) and self.__mU.exists(taxNamePath):
            tD = self.__mU.doImport(taxNamePath, fmt="pickle")
            nD = self.__mU.doImport(taxNodePath, fmt="pickle")
            mD = self.__mU.doImport(taxMergedNodePath, fmt="pickle")
            logger.debug("Taxonomy name length %d node length %d", len(tD), len(nD))
        else:

            nmL, ndL, mergeL = self.__fetchFromSource(urlTarget, taxDirPath)
            tD = self.__extractNames(nmL)
            ok = self.__mU.doExport(taxNamePath, tD, fmt="pickle")
            #
            nD = self.__extractNodes(ndL)
            ok = self.__mU.doExport(taxNodePath, nD, fmt="pickle") and ok
            #
            mD = self.__mergedTaxids(mergeL)
            ok = self.__mU.doExport(taxMergedNodePath, mD, fmt="pickle") and ok
        #
        return tD, nD, mD

    def __mergedTaxids(self, rowL):
        """Extract taxonomy names and synonyms from NCBI taxonomy database dump file row list."""
        tD = {}
        try:
            for tV in rowL:
                if len(tV) < 2:
                    continue
                taxId = int(tV[0])
                mergedTaxId = int(tV[2])
                tD[taxId] = mergedTaxId

            logger.debug("Taxon merged dictionary length %d \n", len(tD))
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return tD

    def __extractNames(self, rowL):
        """Extract taxonomy names and synonyms from NCBI taxonomy database dump file row list."""
        # for matched enclosing quotes
        # reQu = re.compile(r'''(['\"])(.*?)\1$''')
        tD = {}
        try:
            # csvL = []
            for tV in rowL:
                if len(tV) < 7:
                    continue
                taxId = int(tV[0])
                name = tV[2].strip("'")
                #
                # if reQu.match(name):
                #    name = name.strip("'")
                #    logger.info("Matched quoted name: %r" % name)
                #
                nameType = tV[6]
                # csvL.append({'t': taxId, 'name': name, 'type': nameType})
                #
                if nameType in ["scientific name", "common name", "synonym", "genbank common name", "equivalent name", "acronym", "genbank acronym"]:
                    if taxId not in tD:
                        tD[taxId] = {}
                    if nameType in ["scientific name"]:
                        tD[taxId]["sn"] = name
                        continue
                    # take first common name
                    if nameType in ["common name"] and "alt" not in tD[taxId]:
                        tD[taxId]["alt"] = name
                    #
                    tD[taxId].setdefault("cn", []).append(name)
                else:
                    pass
            logger.debug("Taxonomy dictionary length %d \n", len(tD))
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return tD

    def __extractNodes(self, rowL):
        """Extract taxonomy parent relationships from NCBI taxonomy database dump row list."""
        nD = {}
        try:
            for fields in rowL:
                taxid = int(fields[0].strip())
                parentTaxid = int(fields[2].strip())
                rank = str(fields[4]).strip()
                nD[taxid] = (parentTaxid, rank)
            logger.debug("Taxonomy parent dictionary length %d \n", len(nD))
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return nD

    # def __fetchFromSource(self, urlTarget, taxDirPath):
    #     """  Fetch the ncbi taxonomy dump and extract name and node data.
    #     """
    #     _, fn = os.path.split(urlTarget)
    #     #
    #     nmL = self.__mU.doImport(urlTarget, fmt="tdd", rowFormat="list", tarMember="names.dmp", uncomment=False)
    #     ndL = self.__mU.doImport(os.path.join(taxDirPath, fn), fmt="tdd", rowFormat="list", tarMember="nodes.dmp", uncomment=False)
    #     mergeL = self.__mU.doImport(os.path.join(taxDirPath, fn), fmt="tdd", rowFormat="list", tarMember="merged.dmp", uncomment=False)
    #     # deleteL = self.__mU.doImport(os.path.join(taxDirPath, fn), fmt='tdd', rowFormat='list', tarMember='delnodes.dmp')
    #     return nmL, ndL, mergeL
    #
    def __fetchFromSource(self, urlTarget, taxDirPath):
        """Fetch the ncbi taxonomy dump and read name and node data (extract all members)"""
        logger.info("Fetch taxonomy data from source %s in %s", urlTarget, taxDirPath)
        #
        fileU = FileUtil()
        _, fn = os.path.split(urlTarget)
        #
        ok = False
        try:
            tarPath = os.path.join(taxDirPath, fn)
            ok1 = fileU.get(urlTarget, tarPath)
            ok2 = fileU.unbundleTarfile(tarPath, dirPath=taxDirPath)
            ok = ok1 & ok2
        except Exception as e:
            logger.exception("Failing taxonomy fetch from %r with %s", urlTarget, str(e))
        logger.info("%r fetch status is %r", urlTarget, ok)
        # ----  fallback ----
        if not ok:
            logger.info("Fetching taxonomy data from fallback to %s", taxDirPath)
            ok3 = True
            for fn in ["names", "nodes", "merged"]:
                filePath = os.path.join(taxDirPath, "%s.dmp.gz" % fn)
                urlFallback = "https://github.com/rcsb/py-rcsb_exdb_assets/raw/master/fall_back/NCBI/%s.dmp.gz" % fn
                ok1 = fileU.get(urlFallback, filePath)
                ok2 = fileU.uncompress(filePath, outputDir=taxDirPath)
                ok3 = ok1 and ok2 and ok3
            logger.info("Taxonomy fallback fetch status is %r", ok3)
        #
        nmL = self.__mU.doImport(os.path.join(taxDirPath, "names.dmp"), fmt="tdd", rowFormat="list", uncomment=False)
        ndL = self.__mU.doImport(os.path.join(taxDirPath, "nodes.dmp"), fmt="tdd", rowFormat="list", uncomment=False)
        mergeL = self.__mU.doImport(os.path.join(taxDirPath, "merged.dmp"), fmt="tdd", rowFormat="list", uncomment=False)
        # deleteL = self.__mU.doImport(os.path.join(taxDirPath, 'delnodes.dmp'), fmt='tdd', rowFormat='list')

        return nmL, ndL, mergeL
