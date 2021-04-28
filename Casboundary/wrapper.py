import numpy as np

class ClassifierWrapper(object):
    
    def __init__(self, clf, alpha=3.0):
        self.clf = clf
        self.alpha = alpha
    
    def fit(self, X, y):
        self.clf.fit(X, y)
    
    def fit_proba(self, X, y):
        probs = self.clf.predict_proba(X)
        y_pred = self.clf.classes_[np.argmax(probs, axis=1)]
        self.thresholds_ = np.zeros(len(self.clf.classes_))

        for i, c in enumerate(self.clf.classes_):
            idx = np.where(y == c)[0]
            assert np.all(y[idx] == c)

            probs_c = probs[idx, i]
            probs_c = np.concatenate((probs_c, 1.0 + (1.0 - probs_c)))
            std_c = np.std(probs_c, ddof=1)
            self.thresholds_[i] = max(0.5, 1.0 - self.alpha * std_c)
    
    def predict(self, X):
        probs = self.clf.predict_proba(X)
        y_pred = np.argmax(probs, axis=1)
        y_pred_proba = np.max(probs, axis=1)
        thresholds = self.thresholds_[y_pred]
        y_pred = np.where(y_pred_proba >= thresholds, y_pred, -1)
        return np.array([self.clf.classes_[yp] if yp >= 0 else 'unknown' for yp in y_pred])
    
    def set_params(self, **params):
        self.clf.set_params(**params)
