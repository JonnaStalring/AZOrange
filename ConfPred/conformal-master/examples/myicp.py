import numpy as np
from sklearn.linear_model import LogisticRegression


class ICP:
    def __init__(self):
        self.lr = LogisticRegression()

    def fit(self, X, y):
        self.model = self.lr.fit(X, y)

    def calibrate(self, X, y):
        self.classes = np.unique(y)
        ncs = 1 - self.model.predict_proba(X)
        self.cal_scores = np.array(sorted(p[yi] for p, yi in zip(ncs, y)))

    def predict(self, X, significance=0.1):
        p = np.zeros((X.shape[0], self.classes.size))
        for j in range(X.shape[0]):
            ncs = 1 - self.model.predict_proba(X[j:j+1])[0]
            for i, c in enumerate(self.classes):
                nc = ncs[i]
                cal_scores = self.cal_scores
                n_cal = cal_scores.size
                idx_left = np.searchsorted(cal_scores, nc, 'left')
                p[j, i] = (n_cal - idx_left + 1) / (n_cal + 1)
        return p > significance


if __name__ == '__main__':
    import Orange
    iris = Orange.data.Table('iris')
    ind = np.random.permutation(len(iris))
    X, y = iris.X[ind], iris.Y[ind]
    icp = ICP()
    icp.fit(X[:50], y[:50])
    icp.calibrate(X[50:100], y[50:100])
    pred = icp.predict(X[100:], 0.3)
