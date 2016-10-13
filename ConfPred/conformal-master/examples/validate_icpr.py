import Orange
from Orange.preprocess import Impute
import numpy as np

from sklearn.tree import DecisionTreeRegressor

from nonconformist.icp import IcpRegressor
from nonconformist.nc import RegressorNc, abs_error, abs_error_inv


def split_data(data, n_train, n_test):
    n_train = n_train*len(data)//(n_train+n_test)
    n_test = len(data)-n_train
    ind = np.random.permutation(len(data))
    return data[ind[:n_train]], data[ind[n_train:n_train+n_test]]

data = Orange.data.Table("auto-mpg")
imp = Impute()
data = imp(data)

for sig in np.linspace(0.01, 0.1, 10):
    errs, szs = [], []
    for rep in range(10):
        train, test = split_data(data, 2, 1)
        train, calib = split_data(train, 2, 1)

        icp = IcpRegressor(RegressorNc(DecisionTreeRegressor(), abs_error, abs_error_inv))
        icp.fit(train.X, train.Y)
        icp.calibrate(calib.X, calib.Y)
        pred = icp.predict(test.X, significance=sig)

        acc = sum(p[0] <= y <= p[1] for p, y in zip(pred, test.Y))/len(pred)
        err = 1-acc
        sz = sum(p[1]-p[0] for p in pred)/len(pred)
        errs.append(err)
        szs.append(sz)
    print(sig, np.mean(errs), np.mean(szs))
