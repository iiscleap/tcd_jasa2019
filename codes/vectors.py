#This code generates the vectorized features of every word in the libriSpeech training data

import glob
from gensim.models import Word2Vec,Doc2Vec
from preprocess import clean
import os
import re

#Function to clean every sentence of the text
def prune(data):
        ind1 = 0
        ind2 = len(data)
        text = data.split('\n')
        for i in range(len(text)):
                if 'START OF THIS PROJECT GUTENBERG EBOOK' in text[i]:
                        ind1 = i+1
                if 'END OF THIS PROJECT GUTENBERG EBOOK' in text[i]:
                        ind2 = i

        text = text[ind1:ind2]
        text = ' '.join(text)

        text = re.split('[.?]',text)

        return text

def prune1(data):
        p = re.compile('\s\*\*\*\r\n(.*)\*\*\*\sEND',flags=re.DOTALL)
        string = ' '.join(re.findall(p,data))

        if string=='' or string==' ' :
                p=re.compile('<pre>(.*)</pre>',flags=re.DOTALL)
                string = ' '.join(re.findall(p,data))

	string = re.split('[.?]',string)

        return string


def GetData(filepath):
        words = []
	ct = 1
        for r,d,files in os.walk(filepath):
                for f in files:
                        if f.endswith('.txt'):
                                filepath = os.path.join(r,f)
                                print filepath
				
                                with open(filepath, 'r') as fp:
                                        data = fp.read()
                                        data = prune1(data)
                                        for j in range(len(data)):
                                                text = clean(data[j])
                                                text = [text]
                                                for x in text:
                                                        if len(x)!=0:
                                                                words.append(x.split())
				print 'done: ',ct,'files'
				ct=ct+1
        return words
                                                                       


if __name__ == '__main__':
        sentences = GetData('ascii')
	print '\nbuilding model\n'
	#Build word2vec model. The vectors for each word is calculated with respect to the previous and next word
        model = Word2Vec(sentences,window=3,iter=100,size=128)
        model.save('new_word2vec_model.model')
