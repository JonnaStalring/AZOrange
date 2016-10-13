#!/usr/bin/env python3

"""
Example: inductive conformal classification using DecisionTreeClassifier
"""

import numpy as np

from sklearn.tree import DecisionTreeClassifier
import Orange

from nonconformist.acp import CrossConformalClassifier, AggregatedCp, CrossSampler
from nonconformist.icp import IcpClassifier
from nonconformist.nc import ProbEstClassifierNc, inverse_probability

# -----------------------------------------------------------------------------
# Setup training, calibration and test indices
# -----------------------------------------------------------------------------
data = Orange.data.Table('iris')
X, y = data.X, data.Y

idx = np.random.permutation(y.size)
train = idx[:idx.size // 3]
calibrate = idx[idx.size // 3:2 * idx.size // 3]
test = idx[2 * idx.size // 3:]

# -----------------------------------------------------------------------------
# Train and calibrate
# -----------------------------------------------------------------------------
icp = IcpClassifier(ProbEstClassifierNc(DecisionTreeClassifier(), inverse_probability))
icp.fit(X[train, :], y[train])
icp.calibrate(X[calibrate, :], y[calibrate])

ccp = CrossConformalClassifier(IcpClassifier(ProbEstClassifierNc(DecisionTreeClassifier(), inverse_probability)))
ccp.fit(X[train, :], y[train])

acp = AggregatedCp(IcpClassifier(ProbEstClassifierNc(DecisionTreeClassifier(), inverse_probability)), CrossSampler())
acp.fit(X[train, :], y[train])

# -----------------------------------------------------------------------------
# Predict
# -----------------------------------------------------------------------------
print('# Inductive')
prediction = icp.predict(X[test, :], significance=0.1)
for pred, actual in zip(prediction[:5], y[test]):
    print(pred, actual)

print('\n# Cross')
prediction = ccp.predict(X[test, :], significance=0.1)
for pred, actual in zip(prediction[:5], y[test]):
    print(pred, actual)

print('\n# Aggre')
prediction = acp.predict(X[test, :], significance=0.1)
for pred, actual in zip(prediction[:5], y[test]):
    print(pred, actual)