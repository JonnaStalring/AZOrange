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
from Orange.distance import Euclidean, MahalanobisDistance
from Orange.regression import LinearRegressionLearner, RandomForestRegressionLearner, KNNRegressionLearner, SVRLearner
from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVC

from cp.classification import InductiveClassifier
from cp.evaluation import run, CrossSampler
from cp.nonconformity import InverseProbability, SVMDistance, KNNDistance, KNNFraction, LOOClassNC, ProbabilityMargin, \
    AbsError, AbsErrorRF, ErrorModelNC, AbsErrorNormalized, LOORegrNC, AbsErrorKNN, AvgErrorKNN
from cp.regression import InductiveRegressor


dataset_class = ['breast-w.tab', 'credit-a.tab', 'hepatitis.tab', 'iono.tab', 'vote.tab']  # 'spambase.tab'
dataset_class = ['benchmark/classification/'+dataset for dataset in dataset_class]

dataset_regr = ['concrete.tab', 'friedman.tab', 'housing.tab', 'laser.tab', 'mortgage.tab']  # 'abalone.tab'
dataset_regr = ['benchmark/regression/'+dataset for dataset in dataset_regr]

classifiers = ['NaiveBayesLearner()', 'LogisticRegressionLearner()', 'SVMLearner(probability=True)',
               'KNNLearner(n_neighbors=10)', 'KNNLearner(n_neighbors=30)',
               'RandomForestLearner(n_estimators=10)', 'RandomForestLearner(n_estimators=30)']
nc_class_str = [
    "InverseProbability(%s)" % m for m in classifiers] + [
    "ProbabilityMargin(%s)" % m for m in classifiers] + [
    "SVMDistance(SVC())"] + [
    "KNNDistance(Euclidean, %d)" % k for k in [1, 10, 30]] + [
    "KNNFraction(Euclidean, %d, weighted=%s)" % (k, w) for k in [10, 30] for w in [False, True]] + [
    "LOOClassNC(NaiveBayesLearner(), Euclidean, %d, relative=%s, include=%s, neighbourhood='%s')" % (k, r, i, n)
    for k in [10, 30] for r in [False, True] for i in [False, True] for n in ['fixed', 'variable']] + [
    "LOOClassNC(LogisticRegressionLearner(), Euclidean, %d, relative=%s, include=%s, neighbourhood='%s')" % (k, r, i, n)
    for k in [10, 30] for r in [False, True] for i in [False, True] for n in ['fixed', 'variable']] + [
    "LOOClassNC(RandomForestLearner(n_estimators=10), Euclidean, %d, relative=%s, include=%s, neighbourhood='%s')" % (k, r, i, n)
    for k in [10, 30] for r in [False, True] for i in [False, True] for n in ['fixed', 'variable']]

nc_class_str2 = [
    "LOOClassNC(LogisticRegressionLearner(), MahalanobisDistance(), %d, relative=%s, include=%s, neighbourhood='%s')" % (k, r, i, n)
    for k in [10, 30] for r in [False, True] for i in [False, True] for n in ['fixed', 'variable']]

regressors = ['LinearRegressionLearner()', 'SVRLearner()',
              'KNNRegressionLearner(n_neighbors=10)', 'KNNRegressionLearner(n_neighbors=30)',
              'RandomForestRegressionLearner(n_estimators=10)', 'RandomForestRegressionLearner(n_estimators=30)']
nc_regr_str = [
    "AbsError(%s)" % m for m in regressors] + [
    "AbsErrorRF(%s, RandomForestRegressor())" % m for m in regressors] + [
    "ErrorModelNC(LinearRegressionLearner(), %s, loo=%s)" % (m, l) for m in regressors for l in [False, True]] + [
    "ErrorModelNC(RandomForestRegressionLearner(n_estimators=10), %s, loo=%s)" % (m, l) for m in regressors for l in [False, True]] + [
    "AvgErrorKNN(Euclidean, %d)" % k for k in [10, 30]] + [
    "AbsErrorKNN(Euclidean, %d, average=%s, variance=%s)" % (k, a, v)
    for k in [10, 30] for a in [False, True] for v in [False, True]] + [
    "AbsErrorNormalized(LinearRegressionLearner(), Euclidean, %d, exp=%s, rf=%s)" % (k, e, r)
    for k in [10, 30] for e in [False, True] for r in ['None', 'RandomForestRegressor()']] + [
    "AbsErrorNormalized(RandomForestRegressionLearner(n_estimators=10), Euclidean, %d, exp=%s, rf=%s)" % (k, e, r)
    for k in [10, 30] for e in [False, True] for r in ['None', 'RandomForestRegressor()']] + [
    "LOORegrNC(LinearRegressionLearner(), Euclidean, %d, relative=%s, include=%s, neighbourhood='%s')" % (k, r, i, n)
    for k in [10, 30] for r in [False, True] for i in [False, True] for n in ['fixed', 'variable']] + [
    "LOORegrNC(RandomForestRegressionLearner(n_estimators=10), Euclidean, %d, relative=%s, include=%s, neighbourhood='%s')" % (k, r, i, n)
    for k in [10, 30] for r in [False, True] for i in [False, True] for n in ['fixed', 'variable']]

nc_regr_str2 = [
    "LOORegrNC(RandomForestRegressionLearner(n_estimators=10), MahalanobisDistance(), %d, relative=%s, include=%s, neighbourhood='%s')" % (k, r, i, n)
    for k in [10, 30] for r in [False, True] for i in [False, True] for n in ['fixed', 'variable']]

def evaluate_ncs(tab, cp, nc_strs, id):
    res = {}
    for nc_str in nc_strs:
        nc = eval(nc_str)
        r = run(cp(nc), 0.1, CrossSampler(tab, 5), rep=5)
        res[nc_str] = r
        print(id, nc_str)
    with open('results/%s.p' % id, 'wb') as f:
        pickle.dump(res, f)

def evaluate_datasets(datasets, cp, nc_str, id):
    res = {}
    for dataset in datasets:
        dataset_id = dataset.split('/')[-1].split('.')[0]
        imp = Orange.preprocess.Impute()
        rc = Orange.preprocess.preprocess.RemoveConstant()
        tab = rc(imp(Table(dataset)))
        nc = eval(nc_str)
        r = run(cp(nc), 0.1, CrossSampler(tab, 5), rep=5)
        res[dataset_id] = r
        print(nc_str, dataset_id)
    print(nc_str.upper())
    with open('results/nc/%d.p' % id, 'wb') as f:
        pickle.dump(res, f)

def evaluate(datasets, nc_strs, cp, ps, args):
    for nc_str in nc_strs:
        print(nc_str)
        p = multiprocessing.Process(target=evaluate_datasets, args=[datasets, cp, nc_str, len(args)])
        p.start()
        ps.append(p)
        args.append(nc_str)

if __name__ == '__main__':
    os.makedirs('results/nc', exist_ok=True)
    ps = []
    args = []
    id = 0
    if len(sys.argv) == 1 or 'class' in sys.argv:
        #evaluate(dataset_class, nc_class_str, InductiveClassifier, ps, args)
        evaluate(dataset_class, nc_class_str2, InductiveClassifier, ps, args)
    if len(sys.argv) == 1 or 'regr' in sys.argv:
        #evaluate(dataset_regr, nc_regr_str, InductiveRegressor, ps, args)
        evaluate(dataset_regr, nc_regr_str2, InductiveRegressor, ps, args)
    print('SPAWNED')
    with open('results/nc/args.p', 'wb') as f:
        pickle.dump(args, f)
    for p in ps:
        p.join()
    print([p.exitcode for p in ps])
    print('DONE')
