import os

from Orange.classification import LogisticRegressionLearner
from Orange.data import Table

from cp.classification import TransductiveClassifier, CrossClassifier
from cp.evaluation import RandomSampler, LOOSampler, ResultsClass, run_train_test
from cp.nonconformity import InverseProbability

tab = Table(os.path.join(os.path.dirname(__file__), './dataSets/MitoToxStdzdRDkitDescForModelling.tab'))
trains, tests = [], []
lo, hi = 10, 40
for rep in range(30):
    train, test = next(RandomSampler(tab, 100, len(tab)-100))
    trains.append(train)
    tests.append(test)
    for a, b in LOOSampler(train[:lo]):
        assert(len(set(inst.get_class() for inst in a)) > 1)
for n in range(lo, hi, 2):
    print(n)
    tcp = TransductiveClassifier(InverseProbability(LogisticRegressionLearner()))
    ccp = CrossClassifier(InverseProbability(LogisticRegressionLearner()), n)
    tr, cr = ResultsClass(), ResultsClass()
    for train, test in zip(trains, tests):
        tr.concatenate(run_train_test(tcp, 0.1, train[:n], test))
        cr.concatenate(run_train_test(ccp, 0.1, train[:n], test))
    print(tr.accuracy(), tr.multiple_criterion(), tr.time())
    print(cr.accuracy(), cr.multiple_criterion(), cr.time())
