from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils.validation import check_X_y, check_array, check_is_fitted
from sklearn.utils.multiclass import unique_labels
from sklearn.utils.extmath import softmax
from sklearn.linear_model import RidgeClassifier, RidgeClassifierCV
from sklearn.neighbors import NearestCentroid
from sklearn.metrics.pairwise import pairwise_distances
import numpy as np
from mlxtend.feature_selection import ColumnSelector

class DummyClassifierReturnOriginal(BaseEstimator, ClassifierMixin):

    def __init__(self):
        return

    def fit(self, X, y):

        # Check that X and y have correct shape
        X, y = check_X_y(X, y)
        # Store the classes seen during fit
        self.classes_ = unique_labels(y)

        self.X_ = X
        self.y_ = y
        self._estimator_type = 'classifier'
        # Return the classifier
        return self

    def predict(self, X):

        # Check is fit had been called
        check_is_fitted(self)

        # Input validation
        X = check_array(X)


        return np.around(np.array(X)[:, 0])
    
    def predict_proba(self, X):
        # Check is fit had been called
        check_is_fitted(self)

        # Input validation
        X = check_array(X)
        
        return np.concatenate([1-X, X], axis=1)

class RidgeClassifierWithProb(RidgeClassifier):
    def predict_proba(self, X):
        d = self.decision_function(X)
        d_2d = np.c_[-d, d]
        return softmax(d_2d)

class RidgeClassifierCVWithProb(RidgeClassifierCV):
    def predict_proba(self, X):
        d = self.decision_function(X)
        d_2d = np.c_[-d, d]
        return softmax(d_2d)

class RidgeClassifierMultiClassWithProb(RidgeClassifier):
    def predict_proba(self, X):
        return softmax(self.decision_function(X))


class NearestCentroidWithProb(NearestCentroid):
    def predict_proba(self, X):
        distances = pairwise_distances(X, self.centroids_, metric=self.metric)
        return softmax(distances)

class ColumnSelectorWPred(ColumnSelector):
    def fit(self, X, y):

        # Check that X and y have correct shape
        X, y = check_X_y(X, y)
        # Store the classes seen during fit
        self.classes_ = unique_labels(y)
        
        super().fit(X, y)

#         self.X_ = X
#         self.y_ = y
#         self._estimator_type = 'classifier'
        # Return the classifier
        
        return self
    def predict(self, X):
        return self.fit_transform(X)
    
    def predict_proba(self, X):
        return self.fit_transform(X)
    