import Orange
import numpy as np

from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import LogisticRegression

from nonconformist.icp import IcpClassifier
from nonconformist.nc import ProbEstClassifierNc, margin

from myicp import ICP


def split_data(data, n_train, n_test):
    n_train = n_train*len(data)//(n_train+n_test)
    n_test = len(data)-n_train
    ind = np.random.permutation(len(data))
    return data[ind[:n_train]], data[ind[n_train:n_train+n_test]]

#data = Orange.data.Table("../data/usps.tab")
data = Orange.data.Table("iris")

for sig in np.linspace(0.0, 0.4, 11):
    errs, szs = [], []
    for rep in range(10):
        #train, test = split_data(data, 7200, 2098)
        train, test = split_data(data, 2, 1)
        train, calib = split_data(train, 2, 1)

        #icp = IcpClassifier(ProbEstClassifierNc(DecisionTreeClassifier(), margin))
        icp = IcpClassifier(ProbEstClassifierNc(LogisticRegression(), margin))
        #icp = ICP()
        icp.fit(train.X, train.Y)
        icp.calibrate(calib.X, calib.Y)
        pred = icp.predict(test.X, significance=sig)

        acc = sum(p[y] for p, y in zip(pred, test.Y))/len(pred)
        err = 1-acc
        sz = sum(sum(p) for p in pred)/len(pred)
        errs.append(err)
        szs.append(sz)
    print(sig, np.mean(errs), np.mean(szs))
