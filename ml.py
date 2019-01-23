from astropy.table import Table
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
from python_codes.ml_supports import *

class features:
    def __init__(self):
        self.feature_names = list(self.create().columns)[:-1]
        self.features = self.create()[self.feature_names].values
        self.targets = self.create()['class'].values
        self.target_names = ['1', '2'] # 1=galaxy, 2=star
    
    def create(self, tab='sdss.fits'):
        cat = Table.read(tab,format='fits')
        cat = pd.DataFrame(np.array(cat))
        cat['class'] = (cat['class']/3).astype(int).astype(str)
        cat.drop_duplicates(['objID'], inplace=True)
        cat.drop(columns=['objID'], inplace=True)
        cat.drop(['ra','dec'],axis=1,inplace=True)
#        cat = cat[cat.MAG_APER14 > 80]
        return cat

def rfc(n_estimators=100,max_features=np.sqrt(len(features().feature_names)).astype("int"), k=10, random_state=1027):
    table = features()
    forest = RandomForestClassifier(n_estimators=n_estimators, max_features=max_features, random_state=random_state)

    pred = kfold_predictions(forest, table.features, table.targets, k)
    
    score = calculate_accuracy(pred, table.targets)
    print("Accuracy score: {}".format(score))

    model_cm = confusion_matrix(y_true=table.targets, y_pred=pred, labels=table.target_names)

    # Plot the confusion matrix using the provided functions.
    plt.figure()
    plot_confusion_matrix(model_cm, classes=table.target_names, normalize=False)
    plt.show()
