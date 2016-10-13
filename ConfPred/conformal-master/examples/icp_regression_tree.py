#!/usr/bin/env python3

"""
Example: inductive conformal regression using DecisionTreeRegressor
"""

import numpy as np

from sklearn.tree import DecisionTreeRegressor
import Orange

from nonconformist.acp import AggregatedCp, CrossSampler
from nonconformist.icp import IcpRegressor
from nonconformist.nc import RegressorNc, abs_error, abs_error_inv

# -----------------------------------------------------------------------------
# Setup training, calibration and test indices
# -----------------------------------------------------------------------------
data = Orange.data.Table('iris')
X, y = data.X[:,:3], data.X[:, 3]

idx = np.random.permutation(y.size)
train = idx[:idx.size // 3]
calibrate = idx[idx.size // 3:2 * idx.size // 3]
test = idx[2 * idx.size // 3:]

# -----------------------------------------------------------------------------
# Train and calibrate
# -----------------------------------------------------------------------------
icp = IcpRegressor(RegressorNc(DecisionTreeRegressor(), abs_error, abs_error_inv))
icp.fit(X[train, :], y[train])
icp.calibrate(X[calibrate, :], y[calibrate])

acp = AggregatedCp(IcpRegressor(RegressorNc(DecisionTreeRegressor(), abs_error, abs_error_inv)), sampler=CrossSampler())
acp.fit(X[train, :], y[train])

# -----------------------------------------------------------------------------
# Predict
# -----------------------------------------------------------------------------
print('# Inductive')
prediction = icp.predict(X[test, :], significance=0.1)
for pred, actual in zip(prediction[:5], y[test]):
    print(pred, actual)

print('\n# Cross')
prediction = acp.predict(X[test, :], significance=0.1)
for pred, actual in zip(prediction[:5], y[test]):
    print(pred, actual)
