## Meetings & reports

**July 18, 2016**

Nonconformities:
* Bugfix in LOO nonconformities (caching nearest neighbours)

Evaluation
* Add accuracy by class
* Add benchmark datasets
* Run methods also measure time
* Fix range widths for empty regression predictions
* Evaluate classification and regression nonconformities on benchmark datasets
* Evaluate best classification and regression nonconformities on QSAR datasets

Documentation:
* Fix LOORegrNC formula
* Fix AbsErorrRF example
* Add PDF documentation

**July 1, 2016**

Nonconformity
* Experimental NC: |y - \hat{y}\frac{i}{sigma}|
* RF variance in AbsErrorNormalized
* Error model normalization with LOO - ErrorModelNC
* Variance normalization based on predictions of trees in a random forest - AbsErrorRF
* Leave-one-out regression nonconformity - LOORegrNC
* Leave-one-out classification nonconformity - LOOClassNC
* Distance to SVM boundary for binary classification - SVMDistance

Classification
* Mondrian option for all types
* LOOClassifier

Regression
* Fix CrossRegressor (random folds)
* LOORegressor

Other
* Speed-up nearest neighbours
* Mahalanobis distance (Orange development version, not released yet)

Documentation
* Minor clarifications

**June 17, 2016**

Nonconformity
* Average absolute error (eq. 5) - AvgErrorKNN
* Normalized absolute error (eq. 10 / Papadopoulos eq. 31, 32) - AbsErrorNormalized
  * see test_abs_error_normalized for evaluation

Evaluation
* Regression Results
  * Interdecile mean
* Classification Results
  * Fixed confidence measure
  * singleton_correct, average confidence and credibility
  * confusion - count predictions with a given actual and predicted class (e.g. for confusion matrix)
* Run evaluation with specific datasets - run_train_test(cp, eps, train, test, calibrate)

Documentation
* Tutorial
* Module-level documentation
* Clarification of sampling used in run() with an example.
* Accuracy clarification

Other
* Datasets for calibration plots in the same folder

**June 3, 2016**

Documentation
* Source code documentation (/doc/_build/html/index.html)
* Examples of use for individual classes and methods

Evaluation
* Calibration plots (/plots)
* Other datasets

TODO
* Average difference in response to KNN (eq. 5)
* Transductive regression
* User documentation

**May 27, 2016**

Non-conformity measures:
* Classification
  * Probability margin
  * Sum of distances to KNN of same vs. KNN of other class
  * (weighted) Fraction of other class within KNN
* Regression
  * Difference to average response of KNN (normalized by average and/or variance)

Evaluation:
* Run testing with a sampler (random, cross, leave-one-out)
* Accuracy over a test set (validity)
* Classification
  * confidence and credibility of a single prediction
  * various test set criterions (multiple, singleton, empty)
* Regression
  * median width, mean width, std. dev., interdecile range

TODO:
* Average difference in response to KNN (eq. 5, tricky inverse)
* Transductive regression
* Calibration plots on standard datasets
* Documentation


**May 13, 2016**

Implemented and validated methods:
* Classification
  * Transductive
  * Inductive 
  * Cross
* Regression
  * ~~Transductive~~ (difficult, useful?)
  * Inductive
  * Cross

Implemented non-conformity scores:
* Classification
  * Inverse probability with arbitrary classifier (probability of not predicting the correct class)
  * Nearest neighbours (distance ratio between nearest neighbours of same and other class)
* Regression
  * Absolute error


**April 15, 2016**

* first tests of nonconformists
  * Classification. Data is iris, usps. Tested with sklearn classification. For non-conformity scores we have used inverse probabilities.
  * Regression. Data is iris (prediction of one of the leaf measurements).
  * Calibration curves (significance vs error rate) works.
  * All the code is in the github.
* comparison conformal prediction uses smoothing, by default. We have not notice any formal description of this particular approach.
* cross-conformite uses stratification with weighted non-conformity scores.
* some implementations in nonconformist do not stringly replicate the methods from the Vovk et al papers. Example is cross-conformity and aggregated conformity (with cross sampler), where one would expect the same results. The results were different, but some changes in code (say, changing from > to >=) can help to get the same results. 
* Bugs. Cross-conformal predictor ignores parameter that defines number of folds. 
* LOO CP cannot be implemented using cross-conformal predictions. (average accross folds and not accross data sets).
* NCS: classification strictly uses either probability of the classifier.

Conclusion: nonconformists has a nice architecture, but does not implement everything that we would need, and some implementation choices are to us unknown. The code, besides the parameters of functions, is not well documented. We are leaning towards our own implementation that would reuse some code from nonconformists but where all decision choices would be clear.


## Experiments in conformal predictions

Currently, we are using (nonconformist)[https://github.com/donlnz/nonconformist] and checking its functional completness and utility.

### Methods

* Classification
  * Transductive (not implemented)
  * Inductive [Single|Cross]
* Regression
  * Transductive (not implemented)
  * Inductive [Single|Cross]

### Nonconformity scores

* Classification
  * margin
  * inverse_probability
* Regression
  * abs_error
  * sign_error
