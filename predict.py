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
        true_pred = sum(1 for actual, pred in zip(self.y_data, self.y_pred) if actual == pred)
        return accuracy, precision, recall, f1, true_pred
    
    def get_conclusion(self, true_pred, accuracy):
        if accuracy > 0.95:
            kekerabatan = 'dekat'
        else:
            kekerabatan = 'jauh'
        conclusion = f"Dari hasil confusion matrix yang diperoleh, prediksi benar sebanyak {true_pred} sequence-class yang ditunjukkan oleh cell berwarna pada confusion matrix, sehingga didapatkan akurasi {accuracy*100}% yang dapat diartikan bahwa input data-test memiliki kekerabatan yang {kekerabatan} dengan data-training yang digunakan."
        return conclusion

    def get_result(self):
        result = {}
        accuracy, precision, recall, f1, true_pred = self.get_metrics()
        result['accuracy'] = round(accuracy, 4)
        result['precision'] = round(precision, 4)
        result['recall'] = round(recall, 4)
        result['f1'] = round(f1, 4)
        
        return result, self.get_matrix(), self.get_conclusion(true_pred, result['accuracy'])