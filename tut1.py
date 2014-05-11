import numpy
from matplotlib import pyplot as plt
class advect:
    def __init__(self,npart=300,u=1.0,dx=1.0):
        x=numpy.zeros(npart)
        x[npart/3:2*npart/3]=1.0;
        self.x=x
        self.u=u
        self.dx=dx
    def get_bc_periodic(self):
        #the first and last cells are ghost cells
        #so the first cell should be equal to the 
        #next-to-last cell, and the last cell should
        #be equal to the second cell
        self.x[0]=self.x[-2]
        self.x[-1]=self.x[1]
    def update(self,dt=1.0):
        self.get_bc_periodic()
        delt=self.x[1:]-self.x[0:-1]
        self.x[1:-1]+=self.u*dt/self.dx*delt[1:]

if __name__=='__main__':
    stuff=advect()
    plt.ion()
    plt.plot(stuff.x)
    plt.show()
    for i in range(0,300):
        stuff.update()
        plt.clf()
        plt.plot(stuff.x)
        plt.draw()
        
