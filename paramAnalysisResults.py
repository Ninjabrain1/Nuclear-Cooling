import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

def plot(data):
	vn = []
	vdi = []
	vD = []
	vval = []
	verr = []
	for d in data:
		vn.append(d[0])
		vdi.append(d[1])
		vD.append(d[2])
		vval.append(d[3])
		verr.append(d[4])
	fig = plt.figure(figsize=(6, 6))
	ax3 = fig.add_subplot(111, projection='3d')
	plt.tight_layout()
	col = ['g' if b else 'r' for b in vval]
	size = [500 if b else 50 for b in vval]
	shape = ['.' if e==None else '$'+e+'$' for e in verr]
	for n_, di_, D_, col_, size_, shape_ in zip(vn, vdi, vD, col, size, shape):
		ax3.scatter([n_], [di_], [D_], c=col_, s=size_, alpha=0.5, marker=shape_)
	plt.xlabel('n')
	plt.ylabel('di')
	ax3.set_zlabel('D')
	plt.show()

def plot2D(data):
	vn = []
	vdi = []
	vval = []
	verr = []
	for d in data:
		vn.append(d[0])
		vdi.append(d[1])
		vval.append(d[3])
		verr.append(d[4])
	fig = plt.figure(figsize=(6, 6))
	col = ['g' if b else 'r' for b in vval]
	size = [500 if b else 50 for b in vval]
	shape = ['.' if e==None else '$'+e+'$' for e in verr]
	for n_, di_, col_, size_, shape_ in zip(vn, vdi, col, size, shape):
		plt.scatter([n_], [di_], c=col_, s=size_, alpha=0.5, marker=shape_)
	plt.xlabel(r'$n$')
	plt.ylabel(r'$d_i$')
	plt.yticks([6*0.001 + j*0.001 for j in range(dim)])
	plt.xticks([20 + i*5 for i in range(dim)])
	plt.tight_layout()
	plt.show()

def readFile(fileName):
	fileFolder = os.path.dirname(os.path.realpath(__file__))
	fileDir = os.path.join(fileFolder, fileName)
	filehandle = open(fileDir, "r")
	s = filehandle.read()
	filehandle.close()
	return s



data = eval(readFile("sparse.txt"))
plot(data)