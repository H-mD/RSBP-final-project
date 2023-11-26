import joblib
import pandas as pd
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score

class Predict:

    def __init__(self, filepath):
        self.filename = filepath
        data = pd.read_table(self.filename)
        cv = joblib.load('./model/cv.joblib')
        classifier = joblib.load('./model/classifier.joblib')

        data['words'] = data.apply(lambda x: self.Kmers_funct(x['sequence']), axis=1)
        data = data.drop('sequence', axis=1)

        texts = list(data['words'])
        for item in range(len(texts)):
            texts[item] = ' '.join(texts[item])

        self.y_data = data.iloc[:, 0].values
        self.x_data = cv.transform(texts)
        self.y_pred = classifier.predict(self.x_data)

    def Kmers_funct(self, seq, size=6):
        return [seq[x:x+size].lower() for x in range(len(seq) - size + 1)]

    def get_metrics(self):
        accuracy = accuracy_score(self.y_data, self.y_pred)
        precision = precision_score(self.y_data, self.y_pred, average='weighted')
        recall = recall_score(self.y_data, self.y_pred, average='weighted')
        f1 = f1_score(self.y_data, self.y_pred, average='weighted')
        return accuracy, precision, recall, f1

    def get_result(self):
        result = []
        matrics = pd.crosstab(pd.Series(self.y_data, name='Actual'), pd.Series(self.y_pred, name='Predicted'))
        accuracy, precision, recall, f1 = self.get_metrics()
        result.append("accuracy = %.3f" % (accuracy))
        result.append("precision = %.3f" % (precision))
        result.append("recall = %.3f" % (recall))
        result.append("f1 = %.3f" % (f1))
        
        return matrics, result