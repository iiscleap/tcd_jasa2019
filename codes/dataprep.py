#Code to prepare sentences for RNN training and validation
from gensim.models import Word2Vec
import htkmfc as htk
import numpy as np
import pickle
import random

path = 'train960/speakers/' #Speakers used for training
dev_list = ['374','2843','5456','7447','7505'] #Speakers used for testing

train_list = open('train960/speakers.list').read().split('\n')

for i in range(len(dev_list)):
	train_list.remove(dev_list[i])

train_list.remove('')

model = Word2Vec.load('training_model.model')
index = set(model.wv.index2word)
change = 0

#Each sentence is a fixed length of 15 words 
w = htk.open('data_for_nn/train_data_15.htk',mode='w',veclen=1921)

arrays = {}

def make_data(n1,n2):

	s1 = train_list[n1]
	s2 = train_list[n2]

	if s1!=s2: #For sentences with a change
		f1 = open(path+s1+'.txt').readlines()
		f2 = open(path+s2+'.txt').readlines()
		#Pick two random sentences
		line1 = f1[random.randint(0,len(f1)-1)].split()
		line2 = f2[random.randint(0,len(f2)-1)].split()

	else: #For sentences with no change
		f1 = open(path+s1+'.txt').readlines()
		ind = random.randint(0,len(f1)-2)
		line1 = f1[ind].split()
		line2 = f1[ind+1].split()

	if s1!=s2:
		change = 1

	else:
		change = 0

	r = random.randint(0,1)

	if s1==s2:
		if r==0:
			if len(line1)>15:
				sentence = line1[1:16]
			elif len(line2) > 15:
				sentence = line2[1:16]
			else:
				sentence = line1[1:11] + line2[1:6]
		else:
			sentence = line1[1:11] + line2[1:6]

	else:
		sentence = line1[1:11] + line2[1:6]

	word = np.empty((128,))
	word = np.array(model[sentence[0].lower()])

	for j in range(1,len(sentence)):
		word = np.concatenate((word,model[sentence[j].lower()]))
	word = np.reshape(word,((1,word.shape[0])))
			
	l = np.array([change])
	l = np.reshape(l,(1,1))

	word = np.hstack((word,l))	

	print word.shape
	if word.shape[1]==1921:
		w.writeall(word)


def get_data():
	for i in range(700000):
		print i
		n1 = random.randint(0,len(train_list)-1)
		n2 = random.randint(0,len(train_list)-1)
		if n1!=n2:
			make_data(n1,n2)
		else:
			if n1 == len(train_list)-1:
				n2 = 0
			else:
				n2 = n1+1
			make_data(n1,n2)

	for i in range(700000):
		print i
		n1 = random.randint(0,len(train_list)-1)
		make_data(n1,n1)	


get_data()
w.close()
