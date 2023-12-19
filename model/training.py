import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.naive_bayes import MultinomialNB
import joblib

def Kmers_funct(seq, size=6):
    return [seq[x:x+size].lower() for x in range(len(seq) - size + 1)]

human_dna = pd.read_table('./model/human.txt')

human_dna['words'] = human_dna.apply(lambda x: Kmers_funct(x['sequence']), axis=1)
human_dna = human_dna.drop('sequence', axis=1)

human_texts = list(human_dna['words'])
for item in range(len(human_texts)):
    human_texts[item] = ' '.join(human_texts[item])

y_human = human_dna.iloc[:, 0].values

cv = CountVectorizer(ngram_range=(4,4))
X = cv.fit_transform(human_texts)

X_train, X_test, y_train, y_test = train_test_split(X, y_human, test_size = 0.20, random_state=42)

classifier = MultinomialNB(alpha=0.19)
classifier.fit(X_train, y_train)

joblib.dump(cv, './model/cv.joblib')
joblib.dump(classifier, './model/classifier.joblib')