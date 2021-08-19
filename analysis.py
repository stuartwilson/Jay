import h5py
import pylab as pl

f=h5py.File('logs/data.h5')
x = f['x'][:]
y = f['y'][:]
X = f['X'][:]
f.close()

F = pl.figure(figsize=(5,5))
f = F.add_subplot(111)
f.scatter(x,y,c=X,s=3)
pl.show()
