from AZutilities import dataUtilities
import unittest
import os
import time
from trainingMethods import AZorngConsensus
from trainingMethods import AZorngCvSVM
from trainingMethods import AZorngCvANN
from trainingMethods import AZorngRF
from AZutilities import dataUtilities

import orange
from AZutilities import evalUtilities
from AZutilities import miscUtilities
import AZOrangeConfig as AZOC
import AZorngTestUtil
import orngImpute


def float_compare(a, b):
    return abs(a-b)<0.00001

class ConsensusClassifierTest(AZorngTestUtil.AZorngTestUtil):
    def setUp(self):
        """Creates the training and testing data set attributes. """
        testDataRegPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Reg_No_metas_Train.tab")
        irisPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/iris2.tab")
        # Read in the data
        self.DataReg = dataUtilities.DataTable(testDataRegPath)
        self.irisData = dataUtilities.DataTable(irisPath)

    def test_saveloadClass(self):
        """Test the save/load for a classification model - Using Majority"""
        learners = [AZorngRF.RFLearner(), AZorngCvSVM.CvSVMLearner(), AZorngCvANN.CvANNLearner()]

        learner = AZorngConsensus.ConsensusLearner(learners = learners)
        classifier = learner(self.irisData)

        self.assertEqual(len(classifier.domain),len(self.irisData.domain))
        self.assertEqual(len(classifier.imputeData) , len(classifier.domain))
        self.assertEqual(len(classifier.basicStat), len(classifier.domain))
        self.assertEqual(classifier.NTrainEx, len(self.irisData))

        predictions = []
        for ex in self.irisData:
            predictions.append(classifier(ex))
        
        scratchdir = miscUtilities.createScratchDir(desc="ConsensusSaveLoadTest")
        classifier.write(os.path.join(scratchdir,"./CM.model"))

        predictionsL = []
        Loaded = AZorngConsensus.Consensusread(os.path.join(scratchdir,"./CM.model"))
        self.assertEqual(len(Loaded.domain),len(self.irisData.domain))
        self.assertEqual(len(Loaded.imputeData) , len(Loaded.domain))
        self.assertEqual(len(Loaded.basicStat), len(Loaded.domain)) 
        self.assertEqual(Loaded.NTrainEx, 0.8*len(self.irisData))
        for ex in self.irisData:
            predictionsL.append(Loaded(ex))

        self.assertEqual(predictions,predictionsL)

        miscUtilities.removeDir(scratchdir) 


    def test_saveloadClass2(self):
        """Test the save/load for a classification model - Using probabilities average"""
        learners = [AZorngRF.RFLearner(), AZorngCvANN.CvANNLearner()]

        learner = AZorngConsensus.ConsensusLearner(learners = learners)
        classifier = learner(self.irisData)
        predictions = []
        for ex in self.irisData:
            predictions.append(classifier(ex))

        scratchdir = miscUtilities.createScratchDir(desc="ConsensusSaveLoadTest")
        classifier.write(os.path.join(scratchdir,"./CM.model"))

        predictionsL = []
        Loaded = AZorngConsensus.Consensusread(os.path.join(scratchdir,"./CM.model"))
        for ex in self.irisData:
            predictionsL.append(Loaded(ex))

        self.assertEqual(predictions,predictionsL)
        self.assertEqual(len(Loaded.domain),len(self.irisData.domain))
        self.assertEqual(len(Loaded.imputeData) , len(Loaded.domain))
        self.assertEqual(len(Loaded.basicStat), len(Loaded.domain))
        self.assertEqual(Loaded.NTrainEx, len(self.irisData) - int(0.2 * len(self.irisData)))
        #self.assertEqual(Loaded.NTrainEx, len(self.irisData))

        miscUtilities.removeDir(scratchdir)

       
    def test_saveloadReg(self):
        """Test the save/load for a regression model - Using average of N classifiers"""
        learners = [AZorngRF.RFLearner(), AZorngCvSVM.CvSVMLearner(), AZorngCvANN.CvANNLearner()]

        learner = AZorngConsensus.ConsensusLearner(learners = learners)
        classifier = learner(self.DataReg)
        predictions = []
        for ex in self.DataReg:
            predictions.append(classifier(ex))

        scratchdir = miscUtilities.createScratchDir(desc="ConsensusSaveLoadTest")
        classifier.write(os.path.join(scratchdir,"./CM.model"))

        predictionsL = []
        Loaded = AZorngConsensus.Consensusread(os.path.join(scratchdir,"./CM.model"))
        for ex in self.DataReg:
            predictionsL.append(Loaded(ex))

        self.assertEqual([round(pred.value,4) for pred in predictions],
                         [round(pred.value,4) for pred in predictionsL],
                         "Loaded model predictions differ: Pred. 1 (saved/loaded):"+str(predictions[0])+" / "+str(predictionsL[0]))

        self.assertEqual(len(Loaded.domain),len(self.DataReg.domain))
        self.assertEqual(len(Loaded.imputeData) , len(Loaded.domain))
        self.assertEqual(len(Loaded.basicStat), len(Loaded.domain))
        self.assertEqual(Loaded.NTrainEx, len(self.DataReg)*0.8)

        miscUtilities.removeDir(scratchdir)
 

    def test_FeedLearnersReg(self):
        """Test the creation of Consensus feeding Learners for regression"""
        # The Learners can be individualy costumized before passing them to the Consensus
        learners = [AZorngCvSVM.CvSVMLearner(), AZorngCvANN.CvANNLearner(), AZorngRF.RFLearner()]

        # Passing now the learnersObj instead
        learner = AZorngConsensus.ConsensusLearner(learners = learners)
        classifier = learner(self.DataReg)
        predictions = []
        for ex in self.DataReg:
            predictions.append(classifier(ex))

        scratchdir = miscUtilities.createScratchDir(desc="ConsensusSaveLoadTest")
        classifier.write(os.path.join(scratchdir,"./CM.model"))

        predictionsL = []
        Loaded = AZorngConsensus.Consensusread(os.path.join(scratchdir,"./CM.model"))
        for ex in self.DataReg:
            predictionsL.append(Loaded(ex))

        self.assertEqual([round(pred.value,4) for pred in predictions],[round(pred.value,4) for pred in predictionsL],"Loaded model predictions differ: Pred. 1 (saved/loaded):"+str(predictions[0])+" / "+str(predictionsL[0]))

        self.assertEqual(len(Loaded.domain),len(self.DataReg.domain))
        self.assertEqual(len(Loaded.imputeData) , len(Loaded.domain))
        self.assertEqual(len(Loaded.basicStat), len(Loaded.domain))
        self.assertEqual(Loaded.NTrainEx, len(self.DataReg))

        miscUtilities.removeDir(scratchdir)


    def test_FeedClassifiersClass(self):
        """Test the creation of Consensus feeding Classifiers"""

        learners = [AZorngCvSVM.CvSVMLearner(), AZorngCvANN.CvANNLearner(), AZorngRF.RFLearner()]
        classifiers = [l(self.irisData) for l in learners]

        classifier = AZorngConsensus.ConsensusClassifier(classifiers = classifiers)

        predictions = []
        for ex in self.irisData:
            predictions.append(classifier(ex))

        scratchdir = miscUtilities.createScratchDir(desc="ConsensusSaveLoadTest")
        classifier.write(os.path.join(scratchdir,"./CM.model"))

        predictionsL = []
        Loaded = AZorngConsensus.Consensusread(os.path.join(scratchdir,"./CM.model"))
        for ex in self.irisData:
            predictionsL.append(Loaded(ex))

        self.assertEqual(predictions,predictionsL)

        self.assertEqual(len(Loaded.domain),len(self.irisData.domain))
        self.assertEqual(len(Loaded.imputeData) , len(Loaded.domain))
        self.assertEqual(len(Loaded.basicStat), len(Loaded.domain))
        self.assertEqual(Loaded.NTrainEx, len(self.irisData))

        miscUtilities.removeDir(scratchdir)


    def test_FeedClassifiersReg(self):
        """Test the feeding of regression classifiers """
        #DataSet = dataUtilities.DataTable("/home/palmeida/dev/OpenAZOTesteInstall/tests/source/data/linearTrain.tab")
        DataSet = self.DataReg
        learners = [AZorngCvSVM.CvSVMLearner(), AZorngCvANN.CvANNLearner(), AZorngRF.RFLearner()]
        classifiers = [l(DataSet) for l in learners]

        classifier = AZorngConsensus.ConsensusClassifier(classifiers = classifiers)
        predictions = []
        for ex in DataSet:
            predictions.append(classifier(ex))

        scratchdir = miscUtilities.createScratchDir(desc="ConsensusSaveLoadTest")
        classifier.write(os.path.join(scratchdir,"./CM.model"))

        predictionsL = []
        Loaded = AZorngConsensus.Consensusread(os.path.join(scratchdir,"./CM.model"))
        for ex in DataSet:
            predictionsL.append(Loaded(ex))
        self.assertEqual([round(pred.value,4) for pred in predictions],
                         [round(pred.value,4) for pred in predictionsL],
                         "Loaded model predictions differ: Pred. 1 (saved/loaded):"+str(predictions[0])+" / "+str(predictionsL[0]))

        self.assertEqual(len(Loaded.domain),len(DataSet.domain))
        self.assertEqual(len(Loaded.imputeData) , len(Loaded.domain))
        self.assertEqual(len(Loaded.basicStat), len(Loaded.domain))
        self.assertEqual(Loaded.NTrainEx, len(DataSet))

        miscUtilities.removeDir(scratchdir)


    def test_FeedLearnersClass(self):
        """Test the creation of Consensus feeding Learners for classification"""
        #The Learners can be individualy costumized before passing them to the Consensus
        learners = [AZorngCvSVM.CvSVMLearner(), AZorngCvANN.CvANNLearner(), AZorngRF.RFLearner()]

        #Passing now the learnersObj instead
        learner = AZorngConsensus.ConsensusLearner(learners = learners)
        classifier = learner(self.irisData)
        predictions = []
        for ex in self.irisData:
            predictions.append(classifier(ex))

        scratchdir = miscUtilities.createScratchDir(desc="ConsensusSaveLoadTest")
        classifier.write(os.path.join(scratchdir,"./CM.model"))

        predictionsL = []
        Loaded = AZorngConsensus.Consensusread(os.path.join(scratchdir,"./CM.model"))
        for ex in self.irisData:
            predictionsL.append(Loaded(ex))

        self.assertEqual(predictions,predictionsL)

        self.assertEqual(len(Loaded.domain),len(self.irisData.domain))
        self.assertEqual(len(Loaded.imputeData) , len(Loaded.domain))
        self.assertEqual(len(Loaded.basicStat), len(Loaded.domain))
        self.assertEqual(Loaded.NTrainEx, len(self.irisData))

        miscUtilities.removeDir(scratchdir)

        
    def test_CreateLearnerWithObjectMapping(self):
        """ Test the creation of learners with an object map """
        # Arrange
        learners = {'firstLearner':AZorngCvSVM.CvSVMLearner(),
                    'secondLearner':AZorngCvANN.CvANNLearner(),
                    'thirdLearner':AZorngRF.RFLearner()}

        # Act
        learner = AZorngConsensus.ConsensusLearner(learners = learners)

        # Assert
        self.assertEqual(len(learner.learners), len(learners))

    def test_CreateLearnerWithObjectMappingWithoutExpression(self):
        """ Test with name variable mapping defined but not expression given """
        # Arrange
        learners = {'firstLearner':AZorngCvSVM.CvSVMLearner(),
                    'secondLearner':AZorngCvANN.CvANNLearner(),
                    'thirdLearner':AZorngRF.RFLearner()}
        
        learner = AZorngConsensus.ConsensusLearner(learners = learners)

        # Act
        classifier = learner(self.DataReg)
        
        # Assert
        self.assertEqual(classifier, None)


    def test_CanCreateClassifierUsingObjMapping(self):
        """ Test with name variable mapping defined but not expression given """
        # Arrange
        learners = {'firstLearner':AZorngCvSVM.CvSVMLearner(),
                    'secondLearner':AZorngCvANN.CvANNLearner(),
                    'thirdLearner':AZorngRF.RFLearner()}
        
        discreteExpression = ""
        regressionExpression = "(firstLearner + secondLearner + thirdLearner) / 2"
        
        learner = AZorngConsensus.ConsensusLearner(learners = learners, expression = regressionExpression)

        # Act
        classifier = learner(self.DataReg)
        
        # Assert
        self.assertNotEqual(classifier, None)
        self.assertEqual(len(classifier.classifiers), 3)
        self.assertEqual(classifier.expression, regressionExpression)


    def test_AverageNRegressionExpressionUsingObjMap(self):
        """ Test regular expression using average N regression with object map """
        # Arrange
        learners = {'firstLearner':AZorngCvSVM.CvSVMLearner(),
                    'secondLearner':AZorngCvANN.CvANNLearner(),
                    'thirdLearner':AZorngRF.RFLearner()}
        
        # Construct expression learner/classifier
        regressionExpression = "(firstLearner + secondLearner + thirdLearner) / 3"
        expressionLearner = AZorngConsensus.ConsensusLearner(learners = learners, expression = regressionExpression)
        expressionClassifier = expressionLearner(self.DataReg)

        # Construct default learner/classifier
        defaultLearners = [AZorngRF.RFLearner(), AZorngCvANN.CvANNLearner(), AZorngCvSVM.CvSVMLearner()]
        defaultLearner = AZorngConsensus.ConsensusLearner(learners = defaultLearners)
        defaultClassifier = defaultLearner(self.DataReg)
        
        # Act
        expressionPredictions = []
        for ex in self.DataReg:
            expressionPredictions.append(expressionClassifier(ex))

        defaultPredictions = []
        for ex in self.DataReg:
            defaultPredictions.append(defaultClassifier(ex))

        # Assert
        for index in range(len(expressionPredictions)):
            self.assertEqual(True, float_compare(expressionPredictions[index], defaultPredictions[index]))

    def test_CreateLogicalExpressionConsensusLearner(self):
        """ Test creation of logical expression consensus learner """
        # Arrange
        
        # Construct expression learner/classifier
        learners = {'firstLearner':AZorngCvSVM.CvSVMLearner(),
                    'secondLearner':AZorngCvANN.CvANNLearner(),
                    'thirdLearner':AZorngRF.RFLearner()}
        discreteExpression = ["firstLearner == Iris-setosa -> Iris-setosa", "-> Iris-virginica"]
        discreteLearner = AZorngConsensus.ConsensusLearner(learners = learners, expression = discreteExpression)
        discreteClassifier = discreteLearner(self.irisData)

        verifiedLearner = AZorngCvSVM.CvSVMLearner()
        verifiedClassifier = verifiedLearner(self.irisData)
        
        # Act
        result = []
        verifiedResult = []
        for ex in self.irisData:
            result.append(discreteClassifier(ex))
            verifiedResult.append(verifiedClassifier(ex))

        # Assert
        for index, item in enumerate(result):
            if not result[index].value == verifiedResult[index].value:
                print "Not equal on index: ", index
            self.assertEqual(result[index].value, verifiedResult[index].value)

    def test_SaveLoadCustomLogicalExpression(self):
        """ Test save/load functionality with a custom logical expression """
        # Arrange
        
        # Construct expression learner/classifier
        learners = {'firstLearner':AZorngCvSVM.CvSVMLearner(),
                    'secondLearner':AZorngCvANN.CvANNLearner(),
                    'thirdLearner':AZorngRF.RFLearner()}
        discreteExpression = ["firstLearner == Iris-setosa -> Iris-setosa", "-> Iris-virginica"]
        discreteLearner = AZorngConsensus.ConsensusLearner(learners = learners, expression = discreteExpression)
        discreteClassifier = discreteLearner(self.irisData)

        result = []
        for ex in self.irisData:
            result.append(discreteClassifier(ex))
        
        # Act
        scratchdir = miscUtilities.createScratchDir(desc="ConsensusSaveLoadTest")
        discreteClassifier.write(os.path.join(scratchdir,"./CM.model"))

        resultLoaded = []
        loaded = AZorngConsensus.Consensusread(os.path.join(scratchdir,"./CM.model"))
        self.assertNotEqual(loaded, None)
        for ex in self.irisData:
            resultLoaded.append(loaded(ex))

        # Assert
        for index, item in enumerate(result):
            if not result[index].value == resultLoaded[index].value:
                print "Not equal on index: ", index
            self.assertEqual(result[index].value, resultLoaded[index].value)

        self.assertEqual(len(loaded.domain),len(self.irisData.domain))
        self.assertEqual(len(loaded.imputeData) , len(loaded.domain))
        self.assertEqual(len(loaded.basicStat), len(loaded.domain))
        self.assertEqual(loaded.NTrainEx, len(self.irisData))

        miscUtilities.removeDir(scratchdir)

    def test_SaveLoadCustomRegressionExpression(self):
        """ Test save/load custom expression using average N regression with object map """
        # Arrange
        learners = {'firstLearner':AZorngCvSVM.CvSVMLearner(),
                    'secondLearner':AZorngCvANN.CvANNLearner(),
                    'thirdLearner':AZorngRF.RFLearner()}
        
        # Construct expression learner/classifier
        regressionExpression = "(firstLearner + secondLearner + thirdLearner) / 3"
        expressionLearner = AZorngConsensus.ConsensusLearner(learners = learners, expression = regressionExpression)
        expressionClassifier = expressionLearner(self.DataReg)

        # Construct default learner/classifier
        result = []
        for ex in self.DataReg:
            result.append(expressionClassifier(ex))

        # Act
        scratchdir = miscUtilities.createScratchDir(desc="ConsensusSaveLoadTest")
        expressionClassifier.write(os.path.join(scratchdir,"./CM.model"))

        resultLoaded = []
        loaded = AZorngConsensus.Consensusread(os.path.join(scratchdir,"./CM.model"))
        self.assertNotEqual(loaded, None)
        for ex in self.DataReg:
            resultLoaded.append(loaded(ex))

        # Assert
        for index, item in enumerate(result):
            if not float_compare(result[index].value, resultLoaded[index].value):
                print "Not equal on index: ", index
            self.assertEqual(float_compare(result[index].value, resultLoaded[index].value), True)

        self.assertEqual(len(loaded.domain),len(self.DataReg.domain))
        self.assertEqual(len(loaded.imputeData) , len(loaded.domain))
        self.assertEqual(len(loaded.basicStat), len(loaded.domain))
        self.assertEqual(loaded.NTrainEx, len(self.DataReg))

        miscUtilities.removeDir(scratchdir)


if __name__ == "__main__":
    #unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(ConsensusClassifierTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

