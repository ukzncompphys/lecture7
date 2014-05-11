import numpy
from matplotlib import pyplot as plt

class Particles:
    def __init__(self,npart=300,xmax=1.0,u=1.0):
        #distribute particles evenly between 0 and xmax
        self.x=numpy.arange(npart)/(0.0+npart)*xmax
        self.u=u*(xmax-self.x)/xmax
    def update(self,dt=0.01):
        self.x+=self.u*dt
    def get_density(self,dx=0.01):
        #Since the particles don't naturally live on a grid, we need to create one
        #in order to calculate the density.  let's let the grid just run from the minimum to the maximum
        #particle positions.
        xmin=numpy.min(self.x)
        xmax=numpy.max(self.x)
        nbin=numpy.round(1+(xmax-xmin)/dx)
        myind=numpy.round( (self.x-xmin)/dx)
        rho=numpy.zeros(nbin)

        assert(myind.max()<nbin) #we screwed up if this fails
        for i in numpy.arange(0,myind.size):
            rho[myind[i]]+=1.0
        #let's also create the X vector so we can plot the density
        xvec=numpy.arange(0,nbin)*dx+xmin
        return rho,xvec
if __name__=='__main__':
    part=Particles(npart=30000)
    plt.ion()
    plt.plot(part.x)
    plt.show()

    plt.clf()
    for ii in range(0,200):
        part.update(dt=0.01)
        rho,x=part.get_density()
        plt.plot(x,rho)
        plt.draw()
                
