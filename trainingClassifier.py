from parsingfasta import *
import numpy as np
from intervaltree import IntervalTree
from sklearn import svm
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import train_test_split 
import itertools
from pyfasta import Fasta

print("extracting sequence")
seqs = Fasta('data/hg19.fa')
print("extracting TrainingExons")
positiveTrainingExons = read_bed('data/positiveData.bed')
negativeTrainingExons = read_bed('data/negativeData.bed')

#################################
# Build kmer features
k = 3
key_kmers = []
for i in range(1, k+1):
	key_kmers += [ "".join(c) for c in itertools.product('ACGT', repeat=i) ]
print(key_kmers)

################################
#function for extracting features for sequences of given exons
def extractFeatures(exons,seqs):
    kmer_features = []
    for exon in exons:
        exon_data = get_exon_data(exon, seqs)
        kmer_features.append([ exon_data.count(kmer) / (exon.end - exon.start) 
             for kmer in key_kmers ])
    
    return np.array(kmer_features)

################################

def createLabels(nPos, nNeg):
    p=np.ones((nPos,), dtype=np.float)
    n=np.zeros((nNeg,), dtype=np.float)
    return np.concatenate((p, n))

################################
#putting it all together 
print("Building features")
positive_training_features=extractFeatures(positiveTrainingExons,seqs)
negative_training_features=extractFeatures(negativeTrainingExons,seqs)
total_training_data=np.concatenate((positive_training_features,negative_training_features))

trainingLabels=createLabels(len(positiveTrainingExons,),len(negativeTrainingExons))

print("total training data shape : ", total_training_data.shape)
print("training data labels : ", trainingLabels.shape)
#print("total test data shape : ", total_test_data.shape)
#print("test data labels : ", testLabels.shape)

#################################
#training and testing 

#splitting data into training and testing samples, I use big test size as I will not
#Be actually using Test but only training with cross validation
#And becasue my computer can't hnadle too much data I will use only 20% of available data 
X_train, X_test, y_train, y_test = train_test_split(
     total_training_data, trainingLabels, test_size=0.9, random_state=1)

print ( "train data shape :", X_train.shape, "train label shape :", y_train.shape )
print ( "test data shape :", X_test.shape, "test label shape :", y_test.shape )

print("Training")
clf = svm.SVC()
#scores = cross_val_score(clf, X_train[0:10000], y_train[0:10000], cv=5)

#just because it would take too long to train on all the data I am limiting the training
#to only 10,000 points. 

recall=cross_val_score(clf, X_train, y_train, cv=10, scoring='recall')
precision=cross_val_score(clf,X_train, y_train, cv=10,scoring='precision')

f1Score=cross_val_score(clf,X_train, y_train,cv=10, scoring='f1')
aucScore=cross_val_score(clf,X_train, y_train, cv=10,scoring='roc_auc')


print("f1 socre score :", f1Score, "AUC score:" , aucScore)
print("AUPRC:" , precision/recall)

#clf = svm.SVC()
#clf.fit(total_training_data, trainingLabels)
#score = clf.score(total_test_data, testLabels)     
#print("score on test data: " , score)
#sklearnmetric 
