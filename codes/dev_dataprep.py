#Code to prepare the test data
#A running window for every 15 words is used as input. 
#For example 1-15 words, 2-16 words and so on.
import htkmfc as htk
import numpy as np
from gensim.models import Word2Vec
import random

path = '/home/neerajs/work/neeks/spkrChange/dBase/libriSpeech/test-listening/M_5_spkrs_8_pairs_v0/books/change_trans.txt'

data = open(path).read().split('\n')
model = Word2Vec.load('training_model.model')
index = set(model.wv.index2word)

def running_window():
	label = 0
	for i in range(len(data)):
	
		line = data[i].split()
		changeword = line.index('$$') #'$$' denotes the end of sentence 1 

		line.remove('$$')
	
		spk1 = line[0].split('_')[1]
		spk2 = line[0].split('_')[6]

		if spk1!=spk2:
			change = 1
		else:
			change = 0
	
		print changeword,line[changeword]
		print spk1,spk2,change
		print line
		vector = np.empty((1,1920))
		vlen = 0

		for j in range(1,len(line)-15):
			word = np.empty((128,))
			word = np.array(model[line[j].lower()])
			for k in range(j+1,j+15):
				word = np.concatenate((word,model[line[k].lower()]))
			if change==1 and k>=changeword+2 and j<changeword:
				label = 1
			elif change==1 and j<changeword and k<changeword+2:
				label = 0
			elif change==1 and j>=changeword:
				label = 0
			elif change==0:
				label=0
			word = np.reshape(word,((1,word.shape[0])))
			
			if vlen==0:
				vector = word
				vlen = 1
				l = np.array([label])
				l = np.reshape(l,(1,1))
				vector = np.hstack((vector,l))	
				print label
			else:
				print vector.shape,word.shape
	
				l = np.array([label])
				l = np.reshape(l,(1,1))
				word = np.hstack((word,l))
				vector = np.vstack((vector,word))
				print label
	
		w = htk.open('data_for_nn/dev2/'+line[0],mode='w',veclen=vector.shape[1])
		w.writeall(vector)

		w.close()

def make_data():
 for i in range(len(data)):
 
 	line = data[i].split()
	changeword = line.index('$$')

	line.remove('$$')

	spk1 = line[0].split('_')[1]
	spk2 = line[0].split('_')[6]

	if spk1!=spk2:
		change = 1
	else:
		change = 0


 	vector = np.empty((128,))
	word = np.array(model[line[1].lower()])
	
	end = random.randint(3,5)

	for j in range(2,changeword+end):
		word = np.concatenate((word,model[line[j].lower()]))
	word = np.reshape(word,((1,word.shape[0])))

	l = np.array([change])
	l = np.reshape(l,(1,1))

	word = np.hstack((word,l))
	
	w = htk.open('data_for_nn/dev/'+line[0],mode='w',veclen=word.shape[1])
	print word.shape
	w.writeall(word)
	w.close()

running_window()
