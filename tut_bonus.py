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

    def get_state_fft(self,t,dt=1.0):
        #add this routine to jump to an arbitrary time using an FFT
        #from lecture, tutorial 3 etc., we have f(dt,x)=exp(ikx/dx)[1+Cexp(-ik)-C]
        #fft definition is exp(2pi*ikx), so we'll have to take that into account
        #basic idea is that every time we take a step, we pick up another factor of
        #[1+Cexp(-ik)-C], so final state is [1+Cexp(-ik)-C]^nstep
        nstep=t/dt
        xft=numpy.fft.fft(self.x)
        kvec=numpy.arange(0,self.x.size)
        kvec=kvec/(self.x.size+0.0)*2*numpy.pi
        C=dt*self.u/self.dx
        i=numpy.complex(0,1.0)
        fac=1+C*numpy.exp(-i*kvec)-C
        newft=xft*(fac**nstep)
        return numpy.real(numpy.fft.ifft(newft))
    def get_half_k(self,t,dt=1.0):
        nstep=t/dt
        C=dt*self.u/self.dx
        alpha=2.0**(-0.5/nstep)
        cosk=(alpha-1+2*C-2*C*C)/(2*C*(1-C))
        kcrit=numpy.arccos(cosk)
        return kcrit

        

if __name__=='__main__':
    stuff=advect()
    tmax=stuff.x.size
    newx=stuff.get_state_fft(tmax)
    plt.clf()
    #these guys are the same
    plt.plot(stuff.x)
    plt.plot(newx)
    plt.draw()
    

    #now let's call the analytic expression to see the k at which the amplitude is down by a factor of 2
    #and let's make a plot of it.
    tvec=numpy.array([10,20,50,100,200,300,600,900,1200])
    dt=0.5
    kcrit=stuff.get_half_k(tvec,dt)
    

    plt.clf()
    plt.plot(tvec,kcrit)    
    h=plt.gca()
    h.loglog()
    h.set_xlabel('Elapsed time')
    h.set_ylabel('Half-power K')
    h.set_title('Half power wavenumber as a function of time for dt=' + repr(dt))
    plt.draw()
    plt.savefig('half_power_k.png')
    

    #now let's compare some analytic k_1/2's with those observed from the evolution equation 
    tt=55
    kk=stuff.get_half_k(tt,dt)
    stuff.x[:]=0
    stuff.x[150]=1
    newx=stuff.get_state_fft(tt,dt)
    newft=numpy.abs(numpy.fft.fft(newx))
    oldft=numpy.abs(numpy.fft.fft(stuff.x))
    rat=newft/oldft
    i=0
    while rat[i]>0.5:
        i=i+1
    #let's print out the two values.  Due to pi's in the definition of FFT's, we'll need to include them here.
    print 'k observed is ' + repr(i*numpy.pi/stuff.x.size) + ' and expected is ' + repr(kk)

