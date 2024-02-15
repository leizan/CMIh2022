# Author: Lei ZAN <lzan@easyvista.com>
#
# License: GNU General Public License v3.0

from __future__ import print_function

import warnings

import pandas as pd
from numba import jit
from scipy import spatial
from causallearn.utils.cit import CIT
import numpy as np
from src.lib.causal_discovery_from_mixed_data.independence_tests_base import CondIndTest

from src.lib.causal_discovery_from_mixed_data.mixedCmiIEstimatorPython import mixedEstimator
from src.lib.causal_discovery_from_mixed_data.testHypIdentification import TestHypIdentification


class CMIh(CondIndTest):
    r"""Conditional mutual information test based on nearest-neighbor estimator.

    Parameters
    ----------
    knn : int or float, optional (default: 0.2)
        Number of nearest-neighbors which determines the size of hyper-cubes
        around each (high-dimensional) sample point. If smaller than 1, this is
        computed as a fraction of T, hence knn=knn*T. For knn larger or equal to
        1, this is the absolute number.

    shuffle_neighbors : int, optional (default: 5)
        Number of nearest-neighbors within Z for the shuffle surrogates which
        determines the size of hyper-cubes around each (high-dimensional) sample
        point.

    transform : {'ranks', 'standardize',  'uniform', False}, optional
        (default: 'ranks')
        Whether to transform the array beforehand by standardizing
        or transforming to uniform marginals.

    workers : int (optional, default = -1)
        Number of workers to use for parallel processing. If -1 is given
        all processors are used. Default: -1.

    significance : str, optional (default: 'shuffle_test')
        Type of significance test to use. For CMIknn only 'fixed_thres' and
        'shuffle_test' are available.

    **kwargs :
        Arguments passed on to parent class CondIndTest.
    """
    @property
    def measure(self):
        """
        Concrete property to return the measure of the independence test
        """
        return self._measure

    def __init__(self,
                 knn=0.2,
                 shuffle_neighbors=5,
                 significance='shuffle_test',
                 transform='ranks',
                 workers=-1,
                 **kwargs):
        # Set the member variables
        self.knn = knn
        self.shuffle_neighbors = shuffle_neighbors
        self.transform = transform
        self._measure = 'cmi_knn'
        self.two_sided = False
        self.residual_based = False
        self.recycle_residuals = False
        self.workers = workers
        # Call the parent constructor
        CondIndTest.__init__(self, significance=significance, **kwargs)
        # Print some information about construction
        if self.verbosity > 0:
            if self.knn < 1:
                print("knn/T = %s" % self.knn)
            else:
                print("knn = %s" % self.knn)
            print("shuffle_neighbors = %d\n" % self.shuffle_neighbors)

    @jit(forceobj=True)
    def _get_nearest_neighbors(self, array, xyz, knn):
        """

        Parameters
        ----------
        array : array-like
            data array with X, Y, Z in rows and observations in columns

        xyz : array of ints
            XYZ identifier array of shape (dim,).

        knn : int or float
            Number of nearest-neighbors which determines the size of hyper-cubes
            around each (high-dimensional) sample point. If smaller than 1, this
            is computed as a fraction of T, hence knn=knn*T. For knn larger or
            equal to 1, this is the absolute number.

        Returns
        -------
        k_xz, k_yz, k_z : tuple of arrays of shape (T,)
            Nearest neighbors in subspaces.
        """

        array = array.astype(np.float64)
        xyz = xyz.astype(np.int32)

        dim, T = array.shape

        # Add noise to destroy ties...
        array += (1E-6 * array.std(axis=1).reshape(dim, 1)
                  * self.random_state.random((array.shape[0], array.shape[1])))

        if self.transform == 'standardize':
            # Standardize
            array = array.astype(np.float64)
            array -= array.mean(axis=1).reshape(dim, 1)
            std = array.std(axis=1)
            for i in range(dim):
                if std[i] != 0.:
                    array[i] /= std[i]
            # array /= array.std(axis=1).reshape(dim, 1)
            # FIXME: If the time series is constant, return nan rather than
            # raising Exception
            if np.any(std == 0.):
                warnings.warn("Possibly constant array!")
                # raise ValueError("nans after standardizing, "
                #                  "possibly constant array!")
        elif self.transform == 'uniform':
            array = self._trafo2uniform(array)
        elif self.transform == 'ranks':
            array = array.argsort(axis=1).argsort(axis=1).astype(np.float64)

        array = array.T
        tree_xyz = spatial.cKDTree(array)
        epsarray = tree_xyz.query(array, k=[knn+1], p=np.inf,
                                  eps=0., workers=self.workers)[0][:, 0].astype(np.float64)

        # To search neighbors < eps
        epsarray = np.multiply(epsarray, 0.99999)

        # Subsample indices
        x_indices = np.where(xyz == 0)[0]
        y_indices = np.where(xyz == 1)[0]
        z_indices = np.where(xyz == 2)[0]

        # Find nearest neighbors in subspaces
        xz = array[:, np.concatenate((x_indices, z_indices))]
        tree_xz = spatial.cKDTree(xz)
        k_xz = tree_xz.query_ball_point(xz, r=epsarray, eps=0., p=np.inf, workers=self.workers, return_length=True)

        yz = array[:, np.concatenate((y_indices, z_indices))]
        tree_yz = spatial.cKDTree(yz)
        k_yz = tree_yz.query_ball_point(yz, r=epsarray, eps=0., p=np.inf, workers=self.workers, return_length=True)

        if len(z_indices) > 0:
            z = array[:, z_indices]
            tree_z = spatial.cKDTree(z)
            k_z = tree_z.query_ball_point(z, r=epsarray, eps=0., p=np.inf, workers=self.workers, return_length=True)
        else:
            # Number of neighbors is T when z is empty.
            k_z = np.full(T, T, dtype=np.float64)

        return k_xz, k_yz, k_z

    def get_dependence_measure(self, array, xyz):
        """

        Parameters
        ----------
        array : array-like
            data array with X, Y, Z in rows and observations in columns

        xyz : array of ints
            XYZ identifier array of shape (dim,).

        Returns
        -------
        val : float
            Conditional mutual information estimate.
        """
        # print('In function get_dependence_measure:')
        # print('array:')
        # print(array)
        # print('xyz')
        # print(xyz)

        dim, T = array.shape

        if self.knn < 1:
            knn_here = max(1, int(self.knn*T))
        else:
            knn_here = max(1, int(self.knn))

        xind = [0]
        yind = [1]
        zind = []
        isCat = []

        index_data = list(xyz)
        if len(index_data) > 2:
            zind = [i+2 for i in range(len(index_data)-2)]
        column_index = 0
        for ele in array:
            if column_index == 0:
                if isinstance(ele[0], str):
                    isCat.append(column_index)
                    data = np.array([eval(i) for i in ele]).reshape(-1, 1)
                else:
                    data = ele.argsort().argsort().reshape(-1, 1) + 1
            else:
                if isinstance(ele[0], str):
                    isCat.append(column_index)
                    data = np.concatenate((data, np.array([eval(i) for i in ele]).reshape(-1, 1)), axis=1)
                else:
                    data = np.concatenate((data, ele.argsort().argsort().reshape(-1, 1) + 1), axis=1)
            column_index += 1
        data = pd.DataFrame(data)
        val = mixedEstimator(data=data, xind=xind, yind=yind, zind=zind, isCat=isCat)
        return val



    def get_shuffle_significance(self, array, xyz, value,
                                 return_null_dist=False):
        """Returns p-value for nearest-neighbor shuffle significance test.

        For non-empty Z, overwrites get_shuffle_significance from the parent
        class  which is a block shuffle test, which does not preserve
        dependencies of X and Y with Z. Here the parameter shuffle_neighbors is
        used to permute only those values :math:`x_i` and :math:`x_j` for which
        :math:`z_j` is among the nearest niehgbors of :math:`z_i`. If Z is
        empty, the block-shuffle test is used.

        Parameters
        ----------
        array : array-like
            data array with X, Y, Z in rows and observations in columns

        xyz : array of ints
            XYZ identifier array of shape (dim,).

        value : number
            Value of test statistic for unshuffled estimate.

        Returns
        -------
        pval : float
            p-value
        """
        dim, T = array.shape

        # Skip shuffle test if value is above threshold
        # if value > self.minimum threshold:
        #     if return_null_dist:
        #         return 0., None
        #     else:
        #         return 0.

        # max_neighbors = max(1, int(max_neighbor_ratio*T))
        xind = [0]
        yind = [1]
        zind = []
        isCat = []

        index_data = list(xyz)
        if len(index_data) > 2:
            zind = [i+2 for i in range(len(index_data)-2)]
        column_index = 0
        for ele in array:
            if column_index == 0:
                if isinstance(ele[0], str):
                    isCat.append(column_index)
                    data = np.array([eval(i) for i in ele]).reshape(-1, 1)
                else:
                    data = ele.argsort().argsort().reshape(-1, 1) + 1
            else:
                if isinstance(ele[0], str):
                    isCat.append(column_index)
                    data = np.concatenate((data, np.array([eval(i) for i in ele]).reshape(-1, 1)), axis=1)
                else:
                    data = np.concatenate((data, ele.argsort().argsort().reshape(-1, 1) + 1), axis=1)
            column_index += 1

        if len(isCat) == 0:
            fisherz_obj = CIT(data, "fisherz") # construct a CIT instance with data and method name
            pval = fisherz_obj(xind[0], yind[0], zind)
        else:
            data = pd.DataFrame(data)
            testHypIdentification = TestHypIdentification(data=data, xind=xind, yind=yind, zinds=zind,
                                                      isCat=isCat, B=self.sig_samples, kin=self.shuffle_neighbors)
            pval = testHypIdentification.parallelCalcPvalue()

        return pval

    # @jit(forceobj=True)
    # def get_restricted_permutation(self, T, shuffle_neighbors, neighbors, order):
    #
    #     restricted_permutation = np.zeros(T, dtype=np.int32)
    #     used = np.array([], dtype=np.int32)
    #
    #     for sample_index in order:
    #         m = 0
    #         use = neighbors[sample_index, m]
    #
    #         while ((use in used) and (m < shuffle_neighbors - 1)):
    #             m += 1
    #             use = neighbors[sample_index, m]
    #
    #         restricted_permutation[sample_index] = use
    #         used = np.append(used, use)
    #
    #     return restricted_permutation


if __name__ == '__main__':

    import numpy as np

    random_state = np.random.default_rng(seed=42)
    cmi = CMIh(mask_type=None,
                   significance='shuffle_test',
                   fixed_thres=None,
                   sig_samples=1000,
                   sig_blocklength=1,
                   transform='none',
                   knn=0.1,
                   verbosity=0)

    T = 1000
    dimz = 1

    # Continuous data
    z = random_state.standard_normal((T, dimz))
    x = (0.8*z[:,0] + random_state.standard_normal(T)).reshape(T, 1)
    y = (0.8*z[:,0] + random_state.standard_normal(T)).reshape(T, 1)

    print('X _|_ Y')
    print(cmi.run_test_raw(x, y, z=None))
    print('X _|_ Y | Z')
    print(cmi.run_test_raw(x, y, z=z))
