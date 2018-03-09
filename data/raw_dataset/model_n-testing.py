import numpy as np
from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score
from sklearn.externals import joblib


import parser_n_vector as par

train_filename = 'jpred2_train.fasta'
test_filename = 'jpred2_test.fasta'

# Training the model

X_train, y_train = par.input_vecctors(par.parse_fasta_to_dict(train_filename), 5)

clf = SVC(kernel = 'linear', gamma = 0.01, C = 3, cache_size = 1000)
model = clf.fit(X_train, y_train)

score = cross_val_score(clf, X_train, y_train, cv = 10, verbose = True)  ####get scores for the train set
print(score)

joblib.dump(model, 'simple.pkl')

#Testing the model

clf1 = joblib.load('simple.pkl')
X_test, y_test = par.input_vecctors(par.parse_fasta_to_dict(test_filename), 5)

result=clf1.predict(X_test)

count=0
for i in range (len(y_test)):
    if y_test[i] == result[i]:
        count+=1
correct_percentage=count/len(y_test)*100
#print (y_test)
print (result)
print (correct_percentage)

