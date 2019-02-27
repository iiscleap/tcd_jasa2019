#Code to build an LSTM network
#1 LSTM layer with 3 fully connected layers
import keras
import numpy as np
import htkmfc as htk
import scipy.io as sio
import time
import keras
import keras.backend as K
import theano
import h5py
from keras.models import Sequential,Model
from keras.layers import Dense,Convolution2D,Dropout,MaxPooling2D,Input,Flatten,Activation,Merge,Dropout,LSTM
from keras.layers.normalization import BatchNormalization
from keras.layers.wrappers import Bidirectional,TimeDistributed
from keras.utils import np_utils
from keras.optimizers import Adam,SGD
from keras.callbacks import ModelCheckpoint,CSVLogger
from keras.models import load_model
import sys

np.random.seed(2308)
perc=.5

#Function to extract data
def Data_Getter():
	print('Getting and prepping data')
        train = htk.open('data_for_nn/train_data_15.htk')
        train_data=train.getall()
        np.random.shuffle(train_data)
        print('train_done')
        print(train_data.shape)

        val = htk.open('data_for_nn/val_data1.htk')
        val_data=val.getall()
        np.random.shuffle(val_data)
        print('val_done')
        print(val_data.shape)

        Y_train=train_data[:,-1]
        X_train=train_data[:,:-1]
        print(X_train.shape)
        del train_data
        time.sleep(5)
        Y_train=Y_train.reshape(Y_train.shape[0],1)
        Y_train=Y_train.astype(np.int8)
        Y_train=np_utils.to_categorical(Y_train,2)

	
	    Y_val=val_data[:,-1]
        X_val=val_data[:,:-1]
        print(X_val.shape)
        del val_data
        time.sleep(5)
        Y_val=Y_val.reshape(Y_val.shape[0],1)
        Y_val=Y_val.astype(np.int8)
        Y_val=np_utils.to_categorical(Y_val,2)

        print(X_train.shape,X_val.shape,Y_train.shape,Y_val.shape)
        return (X_train,X_val,Y_train,Y_val)



class Master_Controller(keras.callbacks.Callback):
        def on_train_begin(self,logs={}):
                self.val_acc=[]
                self.best_val=0
        def on_epoch_end(self,epoch,logs={}):
                print('\n')
                print(' lr',model.optimizer.lr.get_value())
                model.optimizer.momentum.set_value(model.optimizer.momentum.get_value()+np.float32(.025))
                print('momentum value',model.optimizer.momentum.get_value())
                self.val_acc.append(logs.get('val_acc'))

		print('best_val:',self.best_val)

        #Condition to save the best weights after each epoch and change learning rate
		if self.val_acc[epoch]<=self.best_val:
			file=h5py.File('best_model_text_data_15_pass2.h5','r')
			weight = []
			for i in range(len(file.keys())):
	   			weight.append(file['weight'+str(i)][:])
			model.set_weights(weight)
		
			print('the best model uptil now loaded again : no weight change :(')
			model.optimizer.lr.set_value(K.cast_to_floatx(model.optimizer.lr.get_value()/2))	
		else:
			self.best_val=self.val_acc[epoch]
			file = h5py.File('best_model_text_data_15_pass2.h5','w')
			weight = model.get_weights()
			for i in range(len(weight)):
    				file.create_dataset('weight'+str(i),data=weight[i])
			file.close()
			print('validation accuracy increased and hence model saved')






def mean_var_finder(data):

        mean=np.mean(data,axis=0)
        var=np.var(data,axis=0)

        return (mean,var)

def cnn_reshaper(Data):

        Data=np.reshape(Data,(Data.shape[0],15,128))
        return Data

def normalizer(data,mean,var):
        data=(data-mean)/np.sqrt(var)
        return data

csv_logger=CSVLogger('text_data_15_pass2.log')
lr_red=Master_Controller()
Context,Context_test,Y_train,Y_test=Data_Getter()
print(Context.shape,Y_train.shape,Context_test.shape,Y_test.shape)
Context=cnn_reshaper(Context)
Context_test=cnn_reshaper(Context_test)


mean,var=mean_var_finder(Context)
Context=normalizer(Context,mean,var)

Context_test=normalizer(Context_test,mean,var)
print('saving mean and var')
sio.savemat('mean_var_text_data_15_pass2.mat',{'mean':mean,'var':var})

model=Sequential()
model.add(LSTM(512,input_shape=(15,128)))
model.add(Dense(128,activation='relu'))
model.add(Dense(64,activation='relu'))
model.add(Dense(2,activation='softmax'))

model.summary()

sgd=SGD(lr=.04,momentum=0.5,nesterov=True)
model.compile(loss='categorical_crossentropy',optimizer=sgd,metrics=['accuracy'])
model.fit(Context,Y_train,nb_epoch=20,batch_size=16,validation_data=(Context_test,Y_test),callbacks=[lr_red])
scores=model.predict(Context_test,batch_size=16)

model.save('model_text_data_15_pass2.h5')
print('mat file name is text_data_15_pass2.mat')



