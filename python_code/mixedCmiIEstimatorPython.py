import numpy as np
import itertools
from scipy import special
from joblib import Parallel, delayed, cpu_count
# from CMIEstimator import getEpsDistance
# from scipy import spatial
from scipy.spatial import cKDTree

# The input data is purely continuous and well tailored, the output is the relative distance array
def getDistArray(data):
    """
        calculate the distance from each element in data to another element
    :param data: list of list , each list inside the list is a column
    :return:list of list of list (list of matrices) containing the distance between one point and another point
    """
    N = data.shape[0]
    nDim = data.shape[1]
    # inds = sorted(list(data.columns))
    inds = list(data.columns)
    disArray = []
    for m in inds:
        dataDim = data[m]

        # calculate the distance of Manhattan, return list
        # dataDim est un dataframe, convert to list
        dataDim = list(map(float, list(dataDim)))
        listForDistance = list(np.abs(np.subtract.outer(list(dataDim), list(dataDim))))
        disArray.append([list(element) for element in listForDistance])
    return disArray

def getEpsilonDistance(k, disArray):
    """
        Get epsilon used in Knn for each point
    :param k: the parameter k for k nearest neighbors
    :param disArray: the distance btw points
    :return: the maximum distance for each index of matrices comparing the different dimensions inside disArray
    """
    # takes the most time
    epsilonDisArray = []
    N = len(disArray)
    lengthDisArray0 = len(disArray[0])
    lengthDisArray00 = len(disArray[0][0])
    for elementIndex in range(lengthDisArray0):
        listTemp = []
        for element2Index in range(lengthDisArray00):
            listTemp.append(max([disArray[m][elementIndex][element2Index] for m in range(N)]))
        # sort the list to take the k nearest neighbor
        listTempSorted = sorted(listTemp)
        # take the k éme point
        epsilonDisArray.append(2*sorted(listTemp)[k])

    return epsilonDisArray
def parralelGetEpsilon(k, lenRaws, lenColumns, elementIndex, data):
    listTemp = []
    for element2Index in range(lenRaws):
        maxVal = -float('inf')
        for numColumn in range(lenColumns):
            val = np.abs(data.iloc[elementIndex, numColumn] - data.iloc[element2Index, numColumn])
            if maxVal < val:
                maxVal = val
        listTemp.append(maxVal)
    return 2*sorted(listTemp)[k]
def getEpsDistOptimizedParallel(k, data):
    lenColumns = data.shape[1]
    lenRaws = data.shape[0]
    resultParallel = Parallel(n_jobs=cpu_count()-1)(delayed(parralelGetEpsilon)(
        k, lenRaws, lenColumns, elementIndex, data) for elementIndex in range(lenRaws))
    return resultParallel
def getEpsilonDistanceOptimized(k, data):
    lenColumns = data.shape[1]
    lenRaws = data.shape[0]
    epsilonDisArray = []
    for elementIndex in range(lenRaws):
        listTemp = []
        for element2Index in range(lenRaws):
            maxVal = -float('inf')
            for numColumn in range(lenColumns):
                val = np.abs(data.iloc[elementIndex, numColumn] - data.iloc[element2Index, numColumn])
                if maxVal < val:
                    maxVal = val
            listTemp.append(maxVal)
        epsilonDisArray.append(2*sorted(listTemp)[k])
    return epsilonDisArray
def getEpsilonDistanceFast(k, data):
    tree_xyz = cKDTree(data)
    epsarray = tree_xyz.query(data, k=[k + 1], p=np.inf,
                              eps=0., workers=-1)[0][:, 0].astype(np.float64)
    # distArray = getDistArray(data2)
    # epsilonDis = getEpsilonDistance(k, distArray)
    # print("epsarray", 2*epsarray)
    # print("epsilonDis", epsilonDis)

    return 2*epsarray
def condEntropyEstimator(data, k, dN):
    """
        calculate the conditional entropy estimator
    :param data: the continuous data: data frame
    :param k: the parameter k for k nearest neighbors
    :param dN: number of continuous dimensions
    :return:
    """
    N = data.shape[0]
    if N == 1:
        return 0
    # distArray = getDistArrayParallel(data)
    # distArray2 = getDistArray2(data)

    # distArray = getDistArray(data)
    # epsilonDis = getEpsilonDistance(k, distArray)

    # epsilonDis = getEpsilonDistanceOptimized(k, data)
    # epsilonDis = getEpsDistOptimizedParallel(k, data)
    # cython
    # distArray = getEpsDistance.getDistArray(data)
    # epsilonDis = getEpsDistance.getEpsilonDistance(k, distArray)
    # epsilonDis = getEpsDistance.getEpsilonDistanceOptimized(k, data.to_numpy())
    # print("epsilonDis",  epsilonDis)
    epsilonDis = getEpsilonDistanceFast(k, data.to_numpy())
    if 0 in epsilonDis:
        # delete all null values
        epsilonDis = list(filter(lambda value: value != 0, epsilonDis))
        N = len(epsilonDis)
        if N == 0:
            return 0
    # calculate the entropy using the famous equation
    # if list(epsilonDis)[0] == 0:
    # epsilonDis = list(map(float, epsilonDis))
    entropy = -special.digamma(k) + special.digamma(N) + (dN*sum(np.log(epsilonDis)))/N
    return entropy

def calcDfComb(data, dimDis, allCombinations):
    """
        Calculate the data frames for all combination in combinations
    :param data: the data frame source
    :param dimDis: The index for discrete columns
    :param allCombinations: All possible unique combinations
    :return: list of all combinations
    """
    result = []
    for element in allCombinations:
        dataComb = data
        indiceElement = 0
        for i in dimDis:
            # select all element in the data frame having a specific value
            dataComb = dataComb.loc[dataComb[i] == element[indiceElement]]
            indiceElement += 1

        # To investigate if len(dataComb.index) == 0:
        result.append(dataComb)

    return result

def mixedEntroEstimator(data, dimCon, dimDis):
    """
        Calculate the information entropy for mixed data: the principal function
    :param data: the data frame which contains all the information
    :param dimCon: the indexes in the data frame - continuous ones
    :param dimDis: the indexes in the data frame - discrete ones
    :return: the value of the entropy
    """
    # Input data should be as matrix
    dN = len(dimCon)
    if data is not None:
        N = data.shape[0]
    # dataDis: data discrete taken from data
    # dataCon: data continuous taken from data
    dataDis = []
    dataCon = []
    estimatorCont = 0
    estimatorDisc = 0
    if len(dimCon) != 0:
        # Select only continuous dimensions
        dataCon = data[dimCon]
    if len(dimDis) != 0:
        # Select only discrete dimensions
        dataDis = data[dimDis]

    # if len(dimCon) == 0 and len(dimDis) == 0:
    #     print("The data is Null")

    # If the data is purely continuous!
    if len(dimDis) == 0 and len(dimCon) != 0:
        estimatorCont = condEntropyEstimator(dataCon, max(1, int(0.2*N)), dN)

    # This list takes all unique combinations in the list
    classByDimList = []
    # data frame takes all the combinations of points
    listDfComb = []
    # Calculate the probability of different bins
    probBins = []

    if len(dimDis) != 0:
        for i in dimDis:
            if len(dataDis) > 0:
                classByDimList.append(list(np.unique(dataDis[i])))

        # search all combination of the list
        allCombinations = list(itertools.product(*classByDimList))
        listDfComb = calcDfComb(data, dimDis, allCombinations)

        for element in listDfComb:
            probBins.append(len(element.index)/N)

        for proba in probBins:
            # To delete the possibility having log(0)!
            if proba != 0:
                estimatorDisc -= proba*np.log(proba)

    if len(dimDis) != 0 and len(dimCon) != 0:
        # if the data is mixed
        for i in range(len(probBins)):
            proba = probBins[i]
            if proba != 0 and listDfComb != 0:
                # define k to delete problems
                k = max(1, int(0.2*len(listDfComb[i].index)))
                estimatorCont += proba*condEntropyEstimator(dataCon.iloc[list(listDfComb[i].index), :], k, dN)

    finalResult = estimatorDisc + estimatorCont
    return finalResult

def mixedEstimator(data, xind, yind, zind, isCat):
    """
        This method estimate the conditional mutual information
        between X, Y / Z. Here we calculate all the elements
    :param data: is a dataframe which contain many columns including those of X,Y,Z
    :param xind: Are the indexes of X, list
    :param yind: Are the indexes of Y, list
    :param zinds: Are the indexes of Z, list
    :param isCat: Are the indexes of discrete variables, list
    :return:
    """
    # to delete the possibility to have double values, we use the set the transform it to a list
    xDimCon = list(set(xind).difference(isCat))
    xDimDis = list(set(xind).difference(xDimCon))

    yDimCon = list(set(yind).difference(isCat))
    yDimDis = list(set(yind).difference(yDimCon))

    zDimCon = list(set(zind).difference(isCat))
    zDimDis = list(set(zind).difference(zDimCon))

    conXYZ = list(set(xDimCon + yDimCon + zDimCon))
    disXYZ = list(set(xDimDis + yDimDis + zDimDis))
    # print("here", conXYZ, disXYZ)
    hXYZ = mixedEntroEstimator(data, conXYZ, disXYZ)

    conXZ = list(set(xDimCon + zDimCon))
    disXZ = list(set(xDimDis + zDimDis))
    hXZ = mixedEntroEstimator(data, conXZ, disXZ)

    conYZ = list(set(yDimCon + zDimCon))
    disYZ = list(set(yDimDis + zDimDis))
    hYZ = mixedEntroEstimator(data, conYZ, disYZ)

    conZ = list(set(zDimCon))
    disZ = list(set(zDimDis))
    # calculate hZ, in case of Z null this will give 0 value
    hZ = mixedEntroEstimator(data, conZ, disZ)

    # to parallelize the jobs
    # listeParam = [(conXYZ, disXYZ), (conXZ, disXZ), (conYZ, disYZ), (conZ, disZ)]
    # resultParallel = Parallel(n_jobs=5)(delayed(mixedEntroEstimator)(element[0], element[1], element[2]) for element in listeParam)
    # resultParallel = Parallel(n_jobs=4, prefer="threads")(delayed(mixedEntroEstimator)(data, element[0], element[1]) for element in listeParam)

    # cmi = resultParallel[1] + resultParallel[2] - resultParallel[0] - resultParallel[3]
    cmi = hXZ + hYZ - hXYZ - hZ
    return cmi
