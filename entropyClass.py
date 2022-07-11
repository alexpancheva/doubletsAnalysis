import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import LatentDirichletAllocation as LDA
from sklearn.preprocessing import FunctionTransformer
import seaborn as sb
from scipy.stats import entropy
import random


class EntropyAnalysis(object): 
    def __init__(self, fileName):
        self.counts = pd.read_csv(fileName,header=0,index_col=0)
        self.counts = self.counts.T
        self.transformed = None
        self.weights = None
        self.features = None
        self.model = None
        
    def preprocessing(self):
        scalar = StandardScaler(with_mean=False)
        transformer = FunctionTransformer(np.log1p, validate=True)
        logged = transformer.transform(self.counts)
        transformed = scalar.fit_transform(logged)
        self.transformed = transformed
        
    
    def fit_LDA(self,n_topics,doc_topic_prior=None,topic_word_prior=None):
        lda = LDA(n_components=n_topics,random_state=0,doc_topic_prior=doc_topic_prior, topic_word_prior=topic_word_prior)
        self.model = lda
        output = lda.fit_transform(self.transformed)
        self.weights = pd.DataFrame(output, index=self.counts.index)
        self.features = pd.DataFrame(lda.components_ / lda.components_.sum(axis=1)[:, np.newaxis], columns=self.counts.columns)
        
    
    def calculateEntropy(self):
        return entropy (self.weights.T)
        
        