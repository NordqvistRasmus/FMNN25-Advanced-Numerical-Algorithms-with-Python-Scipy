import numpy as np
import matplotlib.pyplot
import scipy as sp

class CubicSpline(object):

    def __init__(self,deBoorPoints,KnotSeq):

        self.deBoorPoints = deBoorPoints
        self.KnotSeq = KnotSeq

    def __call__(self,u):

        self.u = u        
        self.s_u = self.deBoorsAlgo(u)
        
        return self.s_u
  
    def deBoorsAlgo(self,u):

        
        xVal = np.zeros(len(u))

        yVal = np.zeros(len(u))
 
        
        for i in range(0,len(u)):

            index = self.FirstIndex(u[i])
        
            xVal[i] = self.Blossom(u[i],3,index,0)

            yVal[i] = self.Blossom(u[i],3,index,1)
          
        return np.array([xVal,yVal])

    def Blossom(self,x, rank, i,xy):
        if rank == 0:
            return self.deBoorPoints[i-2,xy]
        if self.KnotSeq[i+1] == self.KnotSeq[i+rank-3]:
            dleft = 0
        else:
            alfa=(self.KnotSeq[i+1] - x)/(self.KnotSeq[i+1] - self.KnotSeq[i+rank-3])
            dleft = alfa * self.Blossom(x, rank-1, i,xy)
        if self.KnotSeq[i+1] == self.KnotSeq[i+rank-3]:
              dright = 0
        else:
            alfa=(self.KnotSeq[i+1] - x)/(self.KnotSeq[i+1] - self.KnotSeq[i+rank-3])
            dright = (1-alfa)* self.Blossom(x, rank-1, i+1,xy)
        return dleft + dright
        

    def FirstIndex(self,u):

        return (self.KnotSeq > u).argmax()-1
 

    def plot(self):
        
        
        xValues = self.s_u[0]
        
        yValues = self.s_u[1]
        #matplotlib.pyplot.figure()
        
        matplotlib.pyplot.plot(xValues,yValues,'k')
        
        #matplotlib.pyplot.plot(self.deBoorPoints[:,0],self.deBoorPoints[:,1],'r--')
        #matplotlib.pyplot.show()
        

    #def makeBasis(knots,j):

     #   controlBase = np.zeros([len(knots)-2,2])
      #  controlBase[j]=[1,1]
     #   knots = np.linspace(0,1,len(controlBase)-2)
      #  knots = np.hstack(([0,0],knots,[1,1]))
     #   return CubicSpline(controlBase,knots)
    
    
    
    def interpol(points):
        d = points #Our interpolationpoints
        knots = np.linspace(0.,1.,len(d)-2)
        knots = np.concatenate([[0,0],knots,[1,1]])

        grevabs=np.zeros([len(knots)-2])
        for i in range(len(knots)-2):
            grevabs[i]=(knots[i]+knots[i+1]+knots[i+2])/3 #Our Greville abscissae
        
        grevabs[-1] =.99 #why is this needed?
    
        VanderMatrix = np.zeros((len(grevabs),len(grevabs))) #Empty matrix for our A in Ax=b

        for i in range(len(grevabs)):
            for j in range(len(grevabs)):
                controlBase = np.zeros([len(knots)-2,2]) #For creating our b-spline basis
                controlBase[j]=[1,1]
                bspline=CubicSpline(controlBase,knots) #Creates the b-spline base
                u=(np.array([grevabs[i]])) 
                VanderMatrix[i,j] = bspline(u)[0] #Evalues the point u using our bspline

                
        x = sp.linalg.solve(VanderMatrix,d[:,0]) #Solves for our vector of x coords
        
        y = sp.linalg.solve(VanderMatrix,d[:,1]) ##Solves for our vector of y coords
    
        controlpoints = np.zeros((len(x),2))
    
        for i in range(len(x)): #Gatheres the coordinates for returning as controlpoints
            controlpoints[i,0]=x[i]
            controlpoints[i,1]=y[i]
        return controlpoints #returns controlpoints
        
        

def main():
   
    deBoorPoints = np.array([[-12.73564, 9.03455],
                                [-26.77725, 15.89208],
                                [-42.12487, 20.57261],
                                [-15.34799, 4.57169],
                                [-31.72987, 6.85753],
                                [-49.14568, 6.85754],
                                [-38.09753, -1e-05],
                                [-67.92234, -11.10268],
                                [-89.47453, -33.30804],
                                [-21.44344, -22.31416],
                                [-32.16513, -53.33632],
                                [-32.16511, -93.06657],
                                [-2e-05, -39.83887],
                                [10.72167, -70.86103],
                                [32.16511, -93.06658],
                                [21.55219, -22.31397],
                                [51.377, -33.47106],
                                [89.47453, -33.47131],
                                [15.89191, 0.00025],
                                [30.9676, 1.95954],
                                [45.22709, 5.87789],
                                [14.36797, 3.91883],
                                [27.59321, 9.68786],
                                [39.67575, 17.30712]])

    KnotSeq = np.linspace(0, 1, 26)
    KnotSeq[ 1] = KnotSeq[ 2] = KnotSeq[ 0]
    KnotSeq[-3] = KnotSeq[-2] = KnotSeq[-1]
    u = np.linspace(0,0.999,1000)

    Values = CubicSpline(deBoorPoints,KnotSeq)

    Values(u)

    Values.plot()


def interpolate():
    #points=np.array([[-10,5],[ -9,-3],[ -8,5],[ -7,-2],[ 6,8],[5,-5],[4,1],[3,-7],[0,7],[1,-3]])
    points=np.array([[0,3],[3,6],[5.5,6.5],[4,0],[2,-3],[0,-6],
                     [-2,-3],[-4,0],[-5.5,6.5],[-3,5.5],[0,3],[0,3]]) 
    
    deBoorPoints=CubicSpline.interpol(points)
    
    knots = np.linspace(0.,1.,len(deBoorPoints)-2)
    KnotSeq = np.concatenate([[0,0],knots,[1,1]])
    
    Values = CubicSpline(deBoorPoints,KnotSeq)
    u = np.linspace(0,0.999,100)
    Values(u)
    Values.plot()

interpolate()
