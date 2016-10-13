import os
import pickle
from multiprocessing import Process
from multiprocessing.pool import Pool

import Orange
import itertools

import sys

import multiprocessing
from Orange.classification import NaiveBayesLearner, LogisticRegressionLearner, RandomForestLearner, KNNLearner, \
    SVMLearner
from Orange.data import Table
from Orange.distance import Euclidean
from Orange.regression import LinearRegressionLearner, RandomForestRegressionLearner, KNNRegressionLearner, SVRLearner
from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVC

from cp.classification import InductiveClassifier
from cp.evaluation import run, CrossSampler, RandomSampler
from cp.nonconformity import InverseProbability, SVMDistance, KNNDistance, KNNFraction, LOOClassNC, ProbabilityMargin, \
    AbsError, AbsErrorRF, ErrorModelNC, AbsErrorNormalized, LOORegrNC, AbsErrorKNN, AvgErrorKNN, ClassNC
from cp.regression import InductiveRegressor


dataset_class = ['AmesHansen6511StdzdRDkitDescForModelling.tab', 'MitoToxStdzdRDkitDescForModelling.tab']
dataset_class = ['dataSets/'+dataset for dataset in dataset_class]

dataset_regr = ['hERGCuratedTropsha5000StdzdRDkitDescForModelling.tab', 'SolubilityMolSARStdzdDupRmDescForModelling.tab']
dataset_regr = ['dataSets/'+dataset for dataset in dataset_regr]

classifiers = ['NaiveBayesLearner()', 'LogisticRegressionLearner()', 'SVMLearner(probability=True)',
               'RandomForestLearner(n_estimators=30)']
nc_class_str = [
    "InverseProbability(%s)" % m for m in classifiers] + [
    "LOOClassNC(LogisticRegressionLearner(), Euclidean, 10, relative=False, include=True, neighbourhood='variable')"]

regressors = ['LinearRegressionLearner()', 'SVRLearner()', 'RandomForestRegressionLearner(n_estimators=30)']
nc_regr_str = [
    "AbsError(%s)" % m for m in regressors] + [
    "AbsErrorRF(%s, RandomForestRegressor())" % m for m in regressors] + [
    "LOORegrNC(LinearRegressionLearner(), Euclidean, 10, relative=True, include=True, neighbourhood='variable')"]

err_rates = [0.05, 0.1, 0.2]

def evaluate_nc_dataset_eps(nc_str, dataset, eps, id):
    nc = eval(nc_str)
    tab = Table(dataset)
    res = None
    for rep in range(100):
        if isinstance(nc, ClassNC):
            r = run(InductiveClassifier(nc), eps, RandomSampler(tab, 2, 1))
        else:
            r = run(InductiveRegressor(nc), eps, RandomSampler(tab, 2, 1))
        if res is None:
            res = r
        else:
            res.concatenate(r)
        print(rep+1, nc_str, dataset, eps)
        with open('results/qsar/%d.p' % id, 'wb') as f:
            pickle.dump((res, rep+1), f)

def evaluate(nc_strs, datasets, epss, ps, args):
    for nc_str in nc_strs:
        for dataset in datasets:
            for eps in epss:
                print(nc_str, dataset, eps)
                p = multiprocessing.Process(target=evaluate_nc_dataset_eps, args=[nc_str, dataset, eps, len(args)])
                p.start()
                ps.append(p)
                args.append([nc_str, dataset, eps])

if __name__ == '__main__':
    os.makedirs('results/qsar', exist_ok=True)

    ps = []
    args = []
    id = 0
    evaluate(nc_class_str, dataset_class, err_rates, ps, args)
    evaluate(nc_regr_str, dataset_regr, err_rates, ps, args)
    print('SPAWNED')
    with open('results/qsar/args.p', 'wb') as f:
        pickle.dump(args, f)
    for p in ps:
        p.join()
    print([p.exitcode for p in ps])
    print('DONE')
