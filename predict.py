import joblib
import pandas as pd
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score

filename = "./test-file/dog.txt"

def Kmers_funct(seq, size=6):
    return [seq[x:x+size].lower() for x in range(len(seq) - size + 1)]

def get_metrics(y_test, y_predicted):
    accuracy = accuracy_score(y_test, y_predicted)
    precision = precision_score(y_test, y_predicted, average='weighted')
    recall = recall_score(y_test, y_predicted, average='weighted')
    f1 = f1_score(y_test, y_predicted, average='weighted')
    return accuracy, precision, recall, f1

cv = joblib.load('./model/cv.joblib')
classifier = joblib.load('./model/classifier.joblib') 

data = pd.read_table(filename)

data['words'] = data.apply(lambda x: Kmers_funct(x['sequence']), axis=1)
data = data.drop('sequence', axis=1)

texts = list(data['words'])
for item in range(len(texts)):
    texts[item] = ' '.join(texts[item])

y_data = data.iloc[:, 0].values

x_data = cv.transform(texts)

y_pred = classifier.predict(x_data)

print("Confusion matrix\n")
print(pd.crosstab(pd.Series(y_data, name='Actual'), pd.Series(y_pred, name='Predicted')))
accuracy, precision, recall, f1 = get_metrics(y_data, y_pred)
print("accuracy = %.3f \nprecision = %.3f \nrecall = %.3f \nf1 = %.3f" % (accuracy, precision, recall, f1))