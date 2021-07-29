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
# 21-Jul-2021 jdw  Make this provider a subclass of StashableBase
##

import collections
import logging
import os.path
from pickle import NONE
import sys

import networkx

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase

logger = logging.getLogger(__name__)


class TaxonomyProvider(StashableBase):
    def __init__(self, **kwargs):
        """ """
        dirName = "NCBI"
        if "cachePath" in kwargs:
            cachePath = os.path.abspath(kwargs.get("cachePath", None))
            self.__taxDirPath = os.path.join(cachePath, dirName)
        else:
            self.__taxDirPath = os.path.abspath(kwargs.get("taxDirPath", "."))
            cachePath, dirName = os.path.split(os.path.abspath(self.__taxDirPath))
        super(TaxonomyProvider, self).__init__(cachePath, [dirName])
        useCache = kwargs.get("useCache", True)
        self.__cleanup = kwargs.get("cleanup", True)
        #
        self.__urlTarget = kwargs.get("ncbiTaxonomyUrl", "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz")
        #
        self.__mU = MarshalUtil(workPath=self.__taxDirPath)
        #
        self.__nodeD = {}
        self.__nameD = {}
        self.__mergeD = {}
        self.__childD = {}
        self.__taxIdToNameD = {}
        self.__graph = None
        #
        self.__nameD, self.__nodeD, self.__mergeD = self.__reload(self.__urlTarget, self.__taxDirPath, useCache=useCache)

    def testCache(self):
        # Lengths name 2133961 node 2133961 merge 54768
        if (len(self.__nameD) > 2100000) and (len(self.__nodeD) > 2100000) and (len(self.__mergeD) > 54000):
            logger.info("Taxonomy lengths name %d node %d merge %d", len(self.__nameD), len(self.__nodeD), len(self.__mergeD))
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

    def getRank(self, taxId):
        try:
            rank = self.__nodeD[int(taxId)][1] if isinstance(taxId, int) and int(taxId) in self.__nodeD else None
        except Exception:
            pass
        return rank

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
        tD = nD = mD = {}
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
            logger.debug("Taxonomy names length %d nodes length %d", len(tD), len(nD))
        elif not useCache:

            nmL, ndL, mergeL = self.__fetchFromSource(urlTarget, taxDirPath)
            tD = self.__extractNames(nmL)
            ok = self.__mU.doExport(taxNamePath, tD, fmt="pickle")
            #
            nD = self.__extractNodes(ndL)
            ok = self.__mU.doExport(taxNodePath, nD, fmt="pickle") and ok
            #
            mD = self.__mergedTaxids(mergeL)
            ok = self.__mU.doExport(taxMergedNodePath, mD, fmt="pickle") and ok
            # Cleanup
            if self.__cleanup:
                fnList = [
                    "citations.dmp",
                    "delnodes.dmp",
                    "division.dmp",
                    "gc.prt",
                    "gencode.dmp",
                    "merged.dmp",
                    "names.dmp",
                    "nodes.dmp",
                    "readme.txt",
                    "taxdump.tar.gz",
                    "names.dmp.gz",
                    "nodes.dmp.gz",
                    "merged.dmp.gz",
                ]
                for fn in fnList:
                    try:
                        fp = os.path.join(taxDirPath, fn)
                        os.remove(fp)
                    except Exception:
                        pass
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

    def __fetchFromSource(self, urlTarget, taxDirPath):
        """Fetch the ncbi taxonomy dump and read name and node data (extract all members)"""
        logger.info("Fetch taxonomy data from source %s in %s", urlTarget, taxDirPath)
        #
        fileU = FileUtil()
        _, fn = os.path.split(urlTarget)
        #
        ok = False
        try:
            ok = True
            tarPath = os.path.join(taxDirPath, fn)
            ok1 = fileU.get(urlTarget, tarPath)
            ok = ok and ok1
            if ok:
                ok1 = fileU.unbundleTarfile(tarPath, dirPath=taxDirPath)
                ok = ok and ok1
        except Exception as e:
            logger.exception("Failing taxonomy fetch from %r with %s", urlTarget, str(e))
        #
        logger.info("Taxonomy primary fetch status (%r) using %r", ok, urlTarget)
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

    def getLowestCommonAncestorGen(self, taxId1, taxId2):
        """Return the lowest common ancestor for the input pair of taxonomy identifiers or None

        Args:
            taxId1 (int): first taxonomy identifier
            taxId2 (int): second taxonomy identifier

        Returns:
            (int): taxonomy identifier for the lowest common ancestor or None
        """
        try:
            if not self.__graph:
                self.__graph = self.__makeGraph()
            return networkx.lowest_common_ancestor(self.__graph, taxId1, taxId2) if self.__graph else None
        except Exception as e:
            logger.exception("Failing for %r %r with %s", taxId1, taxId2, str(e))
        return None

    def getLowestCommonAncestor(self, taxId1, taxId2):
        """Return the lowest common ancestor for the input pair of taxonomy identifiers or None.

        Args:
            taxId1 (int): first taxonomy identifier
            taxId2 (int): second taxonomy identifier

        Returns:
            (int): taxonomy identifier for the lowest common ancestor or None
        """
        if not self.__graph:
            self.__graph = self.__makeGraph()
        rD = dict(networkx.tree_all_pairs_lowest_common_ancestor(self.__graph, root=None, pairs=[(taxId1, taxId2)]))
        return rD[(taxId1, taxId2)] if (taxId1, taxId2) in rD else None

    def getLowestCommonAncestors(self, taxIdPairList):
        """Return the lowest common ancestors for the input pair list of taxonomy identifiers.

        Args:
            taxIdPairList ([type]): [description]

        Returns:
            (dict): {(taxid1,taxid2): lca), ... }
        """
        try:
            if not self.__graph:
                self.__graph = self.__makeGraph()
            return dict(networkx.tree_all_pairs_lowest_common_ancestor(self.__graph, root=None, pairs=taxIdPairList))
        except Exception as e:
            logger.input("Failing for %r with %s", taxIdPairList, str(e))
        return {}

    def __makeGraph(self):
        """Create the networkx graph object..."""
        try:
            graph = networkx.DiGraph()
            graph.name = "Taxonomy"
            for taxId, (parentTaxId, rank) in self.__nodeD.items():
                if taxId == parentTaxId:
                    continue
                graph.add_node(taxId, rank=rank)
                graph.add_edge(parentTaxId, taxId)
            return graph
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return None

    def compareTaxons(self, queryTaxId, refTaxId):
        """Return a summary status, lowest common ancestor and lca rank for input
        query and reference taxonomy identifiers.

                Args:
                    queryTaxId (int): query taxonomy identifier
                    refTaxId (int): reference taxonomy identifier

                Returns:
                    (tuple): status, lcaTaxId, lcaRank (lca: lowerest common ancestor)

                    status: matched, query is ancestor, query is descendant, alternate subspecies
                            alternate variants, alternate strain/serotype/isolate/genotype,
                            alternate strain, orthologous match (by lca), lowest common ancestor.

        """
        try:
            lcaTaxId = None
            status = None
            lcaRank = NONE
            refRank = self.getRank(refTaxId)
            queryRank = self.getRank(queryTaxId)
            if queryTaxId == refTaxId:
                status = "matched"
                lcaTaxId = queryTaxId
                lcaRank = self.getRank(lcaTaxId)
                return status, lcaTaxId, lcaRank
            #
            refTaxIdL = self.getLineage(refTaxId)
            if queryTaxId in refTaxIdL:
                status = "query is ancestor"
                lcaTaxId = queryTaxId
                lcaRank = self.getRank(lcaTaxId)
                return (
                    status,
                    lcaTaxId,
                    lcaRank,
                )
            #
            queryTaxIdL = self.getLineage(queryTaxId)
            if refTaxId in queryTaxIdL:
                status = "query is descendant"
                lcaTaxId = refTaxId
                lcaRank = self.getRank(lcaTaxId)
                return status, lcaTaxId, lcaRank
            #
            lcaTaxId = None
            for taxId in reversed(queryTaxIdL):
                if taxId in refTaxIdL:
                    lcaTaxId = taxId
                    status = "lowest common ancestor"
                    break
            lcaRank = self.getRank(lcaTaxId)
            if lcaRank and lcaRank in ["serotype", "serogroup"]:
                status = "alternate variants"
            elif lcaRank and lcaRank in [
                "clade",
                "class",
                "family",
                "genotype",
                "genus",
                "infraorder",
                "isolate",
                "kingdom",
                "order",
                "parvorder",
                "phylum",
                "subfamily",
                "subgenus",
                "superkingdom",
                "superorder",
                "tribe",
                "species group",
            ]:
                status = "orthologous match (by lca)"
            elif lcaRank and "species" in lcaRank and refRank and "subspecies" in refRank and queryRank and "subspecies" in queryRank:
                status = "alternate subspecies"
            elif (
                lcaRank
                and "species" in lcaRank
                and refRank
                and refRank in ["strain", "serotype", "serovar", "serogroup", "biotype", "isolate", "no rank", "genotype"]
                and queryRank
                and queryRank in ["strain", "serotype", "serovar", "serogroup", "biotype", "isolate", "no rank", "genotype"]
            ):
                status = "alternate strain/serotype/isolate/genotype"
            elif lcaRank and "no rank" in lcaRank and refRank and refRank in ["species", "strain"] and queryRank and queryRank in ["species", "strain"]:
                status = "orthologous match (by lca)"
            elif lcaRank and lcaRank in ["strain"] and refRank and refRank in ["no rank", "strain"] and queryRank and queryRank in ["no rank", "strain"]:
                status = "alternate strain"

        except Exception as e:
            logger.exception("Failing for %r and %r with %s", queryTaxId, refTaxId, str(e))
        return status, lcaTaxId, lcaRank
