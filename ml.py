from astropy.table import Table
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier

class features:
    def __init__(self):
        self.feature_names = list(self.create().columns)[:-1]
        self.features = self.create()[self.feature_names].values
        self.targets = self.create()['class'].values
        self.target_names = ['galaxy', 'star']
    
    def create(self, tab='sdss.fits'):
        cat = Table.read(tab,format='fits')
        cat = pd.DataFrame(np.array(cat))
        cat['class'] = (cat['class']/3).astype(int).astype(str)
        cat.drop_duplicates(['objID'], inplace=True)
        cat.drop(columns=['objID'], inplace=True)
        cat.drop(['ra','dec'],axis=1,inplace=True)
#        cat = cat[cat.MAG_APER14 > 80]
        return cat

def rfc(n_estimators=100,max_features=np.sqrt(len(features().feature_names)).astype("int"),random_state=0):
    feat = features()
    X_train, X_test, y_train, y_test = train_test_split(
        feat.features, feat.targets, random_state=random_state)
    forest = RandomForestClassifier(n_estimators=n_estimators, max_features=max_features, random_state=random_state)
    forest.fit(X_train, y_train)
    print("Accuracy on training set: {:.3f}".format(forest.score(X_train, y_train))) 
    print("Accuracy on test set: {:.3f}".format(forest.score(X_test, y_test)))
    num_g = np.sum(feat.create()['class'] == '1')
    num_s = np.sum(feat.create()['class'] == '2')
    per = num_g/(num_g+num_s)
    print('number of galaxies: {}'.format(num_g))
    print('number of stars: {}'.format(num_s))
    print('percentage of galaxies over the sample: {}'.format(per))
