import joblib
import pandas as pd
import matplotlib.pyplot as plt
from io import BytesIO
import base64
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score
import matplotlib
matplotlib.use('Agg')

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

    def get_matrix(self):
        matrics = pd.crosstab(pd.Series(self.y_data, name='Actual'), pd.Series(self.y_pred, name='Predicted'))
        plt.figure(figsize=(8, 6))
        plt.imshow(matrics, cmap='Blues', alpha=0.3, interpolation='nearest')
        for i in range(len(matrics.index)):
            for j in range(len(matrics.columns)):
                plt.text(j, i, matrics.iloc[i, j], ha='center', va='center', color='Black')
        plt.tick_params(axis='x', which='both', bottom=False, top=True, labelbottom=False, labeltop=True)
        plt.gca().xaxis.set_label_position('top')
        plt.xlabel('Predicted')
        plt.ylabel('Actual')
        plt.title('Confusion Matrix')
        plt.xticks(range(len(matrics.columns)), matrics.columns)
        plt.yticks(range(len(matrics.index)), matrics.index)
        plt.tight_layout()
        img = BytesIO()
        plt.savefig(img, format='png')
        img.seek(0)
        plt.close()

        img = base64.b64encode(img.getvalue()).decode()

        return img

    def get_metrics(self):
        accuracy = accuracy_score(self.y_data, self.y_pred)
        precision = precision_score(self.y_data, self.y_pred, average='weighted')
        recall = recall_score(self.y_data, self.y_pred, average='weighted')
        f1 = f1_score(self.y_data, self.y_pred, average='weighted')
        return accuracy, precision, recall, f1

    def get_result(self):
        result = []
        accuracy, precision, recall, f1 = self.get_metrics()
        result.append("accuracy = %.3f" % (accuracy))
        result.append("precision = %.3f" % (precision))
        result.append("recall = %.3f" % (recall))
        result.append("f1 = %.3f" % (f1))
        
        return result, self.get_matrix()