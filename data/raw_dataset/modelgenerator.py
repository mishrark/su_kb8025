import numpy as np
from sklearn import svm

##############parsing raw fasta file into dictionaries'''

filename = 'jpred2_train.fasta'
data_dict = {}

with open(filename, 'r') as infile:
	for i, line in enumerate(infile):
		line = line.strip()
		if i % 3 == 0:
			id = line[1:]
		elif i % 3 == 1:
			seq = line
		elif i % 3 == 2:
			topo = line
			data_dict[id] = [seq,topo]

#print(data)

# ###########separating sequences and secondary structures into two lists

seq_list = []
ss_list = []

seq_ss = (data_dict.values())

for i in seq_ss:
	seq_list.append(i[0])
	ss_list.append(i[1])

#print(seq_list)
#print(ss_list)				

############Specifying Window size

window_size = 5
pad_size = int((window_size/ 2 )-0.5)
pad = 'X' * pad_size



##################prepares input vector for sklearn'''

aa_dic = {  'A':[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],
                'C':[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,], 
                'D':[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],
                'E':[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],
                'F':[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],
                'G':[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,], 
                'H':[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,], 
                'I':[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,],
                'K':[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,],
                'L':[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,],
                'M':[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,],
                'N':[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,],
                'P':[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,],
                'Q':[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,],
                'R':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,],
                'S':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,],
                'T':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,],
                'V':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,],
                'W':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,],
                'Y':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,],
                'X':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,] }

windowed_seq = []

for seq in seq_list:
	seq = pad + seq +pad

	for i in range(len(seq)):
		window = seq[i:i + window_size]
		if len(window) == window_size:
			windowed_seq.append(window)
#print(windowed_seq)

X_values = []

for i in windowed_seq:
	tmp = []
	for j in i:
		tmp.extend(aa_dic[j])
	X_values.append(tmp)
#print(X_values)

#############Preparing labels

ss = {'H':1, 'E': 0, 'C':-1}

Y_values = []

for item in ss_list:
	for sec_str in item:
		Y_values.append(ss[sec_str])
#print(Y_values)


#Training and Model

X_train = np.array(X_values)
#print(X_train.shape)
y_train = np.array(Y_values)
#print(Y_train.shape)

clf = svm.SVC(gamma = 0.01, C =1)
clf.fit(X_train, y_train)

# Testing the model

#clf.predict(X_test)
#tried on a test sequence and the score is only .46.

