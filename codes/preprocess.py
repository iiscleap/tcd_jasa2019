#Preprocessing file
from nltk.corpus import stopwords
import re

stop_words = ['project','gutenberg','ebook','ebooks','wwwgutenbergorg','online','transcribed','copyright','united','states','email','newsletter','']

def clean(text):
	text = text.lower()
	text = re.sub('[#$%|~\-&\'"`*+=!?;()\{\}^\[\]<>:;_@]', '', text) 
	text = re.sub('/', '', text)
	text = re.sub('http[a-z]*\s','',text)
	text = re.sub('character set encoding\s[a-z]*','',text)
	text = re.sub('\r\n',' ',text)
	text = re.sub('\n',' ',text)
	text = re.sub(' +',' ',text)
	text = re.sub('\d+','',text)
	text = re.sub('chapter\s[a-z]*\s','',text)
	words = [word for word in text.split() if word not in stop_words]
	return ' '.join(words)
