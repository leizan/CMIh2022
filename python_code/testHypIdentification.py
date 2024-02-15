# Importing Pandas as pd
import numpy as np
import pandas as pd
import random
from joblib import Parallel, delayed, cpu_count
# from CMIEstimator.mixedCmiIEstimatorPython import mixedEstimator
from src.lib.causal_discovery_from_mixed_data.mixedCmiIEstimatorPython import mixedEstimator
import multiprocessing as mp

class TestHypIdentification:
    def __init__(self, data, xind, yind, zinds, isCat, kin, B):
        self.xind, self.yind, self.zinds = xind, yind, zinds
        self.nodes = self.xind + self.yind + self.zinds
        self.data = data
        self.isCat = isCat
        self.kin, self.B = kin, B
        self.zDis = list(set(self.zinds).intersection(self.isCat))
        # self.zDis = [1, 2]
        self.processDataIdentification()
        # self.searchK(self.kin)
        # pvalueCalc = self.calcPvalue()
        # self.swapData(3)
        # self.searchKNN('c', 5)
    def processDataIdentification(self):
        """
            This function will process and normalize data
        :return:
        """
        # get only columns with continuous dimension
        nodes = set(self.xind + self.yind + self.zinds)
        continuousNodes = list(nodes.difference(self.isCat))
        # NB rank takes only continous columns ! so it will be a problem if a column is
        # discrete and tou take it as continuous
        self.data[continuousNodes] = self.data[continuousNodes].rank(method="first")
        # print(self.data)
    def pvaluefunc4Parallel(self, Ioriginal):
        dataswapped = self.swapData(self.kin)
        Iswapped = mixedEstimator(data=dataswapped, xind=self.xind, yind=self.yind,
                                  zind=self.zinds, isCat=self.isCat)
        # print("Iswapped", Iswapped)
        pvalueCalc = 1 if Iswapped >= Ioriginal else 0
        return pvalueCalc
    def parallelCalcPvalue(self):
        """
            This function parallelize the code for pvalue calculus
        :return:
        """
        Ioriginal = mixedEstimator(data=self.data, xind=self.xind, yind=self.yind,
                                    zind=self.zinds, isCat=self.isCat)
        resultParallel= Parallel(n_jobs=cpu_count()-1)(delayed(self.pvaluefunc4Parallel)(Ioriginal) for i in range(1, self.B+1))
        pvalueCalc = sum(resultParallel)/self.B
        return pvalueCalc
    def calcPvalue(self):
        """
            This func will generate B locally permuted samples and calculate the p_value
        :return:
        """
        Ioriginal = mixedEstimator(data=self.data, xind=self.xind, yind=self.yind,
                                    zind=self.zinds, isCat=self.isCat)
        pvalueCalc = 0
        for i in range(1, self.B+1):
            dataswapped = self.swapData(self.kin)
            Iswapped = mixedEstimator(data=dataswapped, xind=self.xind, yind=self.yind,
                                        zind=self.zinds, isCat=self.isCat)
            pvalueCalc += 1 if Iswapped >= Ioriginal else 0
        pvalueCalc /= self.B
        return pvalueCalc
    def searchKNN(self, indexPointDf, kin):
        """
            This function calculates the knn for a given index point in the dataframe
            # NB, if many are equal, we take the all possible points even if have more k points
            since the points are the "same" if their distance is equal
            # NB, we should have order operator between discrete values
        :param indexPointDf: the index of the point in the dataframe
        :param k: the k to use for this point
        :return: the list of indexes in df of knn
        """
        # calculate the distance to all other points, in discrete and continuous dimensions
        # the distance is max(distance of point to other point in 1 dimension)
        # in case of discrete one it's 0 or inf, continuous is |a-b|
        listNodeDistance = {}
        # we will use this var to know the cardinal of point that has the same values in the discrete column
        cardSameDiscPoint = 0
        for index in self.data.index:
            boolSameValues = True
            maxDistance = -np.inf
            for node in self.zinds:
                if node in self.zDis:
                    newBool = True if self.data.loc[index, node] == self.data.loc[indexPointDf, node] else False
                    boolSameValues = boolSameValues and newBool
                    distance = 0 if self.data.loc[index, node] == self.data.loc[indexPointDf, node] else np.Inf
                else:
                    distance = np.abs(self.data.loc[index, node] - self.data.loc[indexPointDf, node])
                if distance > maxDistance:
                    maxDistance = distance

            listNodeDistance[index] = maxDistance
            # if we have the same discrete values and we have discrete values, we increase the card
            if boolSameValues is True and len(self.zDis) > 0:
                cardSameDiscPoint += 1
        listResultSorted = sorted(listNodeDistance, key=listNodeDistance.get)
        k = self.calcK(cardSameDiscPoint, kin)
        if k < 1:
            # just in case of invalid k
            k = 1
        if k > len(listResultSorted):
            k = len(listResultSorted)
        valueKnn = listNodeDistance[listResultSorted[k-1]]
        listKnnPoint = []
        # we will go through all values and see if the values is small or equal than valueKnn
        for key in listResultSorted:
            if listNodeDistance[key] <= valueKnn:
                listKnnPoint.append(key)
        return listKnnPoint
    def calcK(self, cardSameDiscPoint, kin):
        """
            This function calculate k = kin if z continuous
            or min(kin, cardSameDiscPoint) else
        :param cardSameDiscPoint: the cardinal of point discrete
        :param kin: the k initial to use
        :return: the final k to use in the method
        """
        if cardSameDiscPoint == 0:
            # z is only continuous
            return kin
        else:
            return min(kin, cardSameDiscPoint)
    # def searchK(self, kin):
        #code pb
        # if len(self.zDis) == 0:
        #     return kin
        # # self.data[self.data.columns[self.zDis]]
        # self.data["discreteZ"] = self.data['First_Name'] + self.data['Last_Name']
        # self.data["discreteZ"] = self.data[self.zDis].apply(
        #     lambda x: ','.join(x.dropna().astype(str)),
        #     axis=1
        # )
        # print(self.data)

        # k = self.calcK(cardSameDiscPoint, kin)
        # if k < 1:
        #     # just in case of invalid k
        #     k = 1
        # if k > len(listResultSorted):
        #     k = len(listResultSorted)
    def swapData(self, kin):
        """
            This function swap the data in the x dim, using nearest neighbors
            without changing other values
        :param kin: the k initial
        :return:
        """
        self.dataSwapped = self.data.copy()
        listIndexesUsed = set()
        for index in self.data.index:
            listKnnPoint = self.searchKNN(index, kin)
            listNodeNotUsed = list(set(listKnnPoint).difference(listIndexesUsed))
            index2UseSwap = random.choice(listKnnPoint)
            if len(listNodeNotUsed) > 0:
                index2UseSwap = random.choice(listNodeNotUsed)
            self.dataSwapped.loc[index, self.xind] = self.data.loc[index2UseSwap, self.xind]
            listIndexesUsed.add(index2UseSwap)
        return self.dataSwapped

# if __name__ == '__main__':
#     testHypIdentification = testHypIdentification()
