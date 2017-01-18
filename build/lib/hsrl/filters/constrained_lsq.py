import scipy.optimize as opt
import numpy as np

class opt:

    def __init__(self,x,data,fx0=None,fxn=None,dfx0=None,dfxn=None):
        """Compute constrained polynomial least-squares fit
           x[i] = independent variable vector
           data = meaurements at points specified in x
           optional constraints:
              fx0  = constrained value of fit at x[0]
              fxn  = constrained value of fit at x[-1]
              dfx0 = constrained value of fitted df/dx at x[0]
              dfxn = constrained value of fitted df/dx at x[-1]
        """      
        self.x = x
        self.data = data
        self.fx0  = fx0
        self.fxn =fxn
        self.dfx0 = dfx0
        self.dfxn = dfxn



    def func(self,coef,sign=1.0):
        f = sign* np.sum((self.data-np.polyval(coef,self.x))**2)
        return f

    def func_deriv(self,coef, sign=1.0):
       """ Derivative of objective function """
       order = len(coef)-1
       df = np.zeros(order+1)
       for i in range(order+1):
            df[i] = sign*2.0*np.sum((self.data-np.polyval(coef,self.x)) * (-self.x**(order-i)))
       return np.array(df)

    def f0_constraint(self,coef):
         order = len(coef)-1
         constraint = 0.0
         for i in range(order+1):
             constraint = constraint + coef[i] * self.x[0]**(order-i)
         constraint = constraint - self.fx0
         return constraint

    def fn_constraint(self,coef):
             #value at end point constraint
             order=len(coef)-1
             constraint =0.0
             for i in range(order+1):
                 constraint = constraint + coef[i] * self.x[-1]**(order-i)
             constraint = constraint - self.fxn
             return constraint
    def dfx0_constraint(self,coef):
             #slope at start point constraint
             n=len(coef)
             constraint = 0.0
             for i in range(n-1):
                 constraint = constraint + (n-i-1)*coef[i]*self.x[0]**(n-i-2)
             constraint = constraint - self.dfx0
             return np.array(constraint)
    def dfxn_constraint(self,coef):
             #slope at end point constraint
             n=len(coef)
             constraint = 0.0
             for i in range(n-1):
                 constraint = constraint + (n-i-1)*coef[i]*self.x[-1]**(n-i-2)
             constraint = constraint - self.dfxn
             return np.array(constraint)
             #return np.array(2*coef[0]*self.x[-1] +coef[1]-self.dfxn)
    def f0_jacobian(self,coef):
         order = len(coef)-1
         jacobian = np.zeros(order+1)
         for i in range(order+1):
             jacobian[i] = self.x[0]**(order-i)
         return jacobian
    def fn_jacobian(self,coef):  
             order = len(coef)-1
             jacobian = np.zeros(order+1)
             for i in range(order+1):
                 jacobian[i] = self.x[-1]**(order-i)
             return jacobian
    def dfx0_jacobian(self,coef):
             order =len(coef)-1
             jacobian = np.zeros(order+1)
             for i in range(order):
                 jacobian[i] = (order-i)*self.x[0]**(order-i-1)
             return jacobian 
    def dfxn_jacobian(self,coef):
             order =len(coef)-1
             jacobian = np.zeros(order+1)
             for i in range(order):
                 jacobian[i] = (order-i)*self.x[-1]**(order-i-1)
             return jacobian
         #return np.array([2*self.x[-1], 1.0,0.0])
     
    def fit(self,coef):
       """Compute constrained polynomial least-squares fit
          coef = first guess polynomial coeficients
               = p[0]*x**n +p[1]*x**(n-1)...... +p[n]*x + p[n+1]
          coef = [ p[0],p[1]........p[n],p[n+1]     
          order of fitted polynomial is determined by length of the guess
       """     
       cons =[]
       if not self.fxn == None:
           cons.append({'type': 'eq',
               'fun' : self.fn_constraint,
               'jac' : self.fn_jacobian})
           
       if not self.fx0 == None:
            cons.append({'type':'eq',
               'fun': self.f0_constraint,
               'jac': self.f0_jacobian})
            
       if not self.dfx0 == None:
            cons.append({'type':'eq',
               'fun': self.dfx0_constraint,
               'jac': self.dfx0_jacobian})
            
       if not self.dfxn == None:
            cons.append({'type':'eq',
               'fun': self.dfxn_constraint,
               'jac': self.dfxn_jacobian})
            
       if self.fx0 == None and self.fxn == None \
                 and self.dfx0 == None and self.dfxn == None:
          #unconstrained fit 
          res = opt.minimize(self.func,coef, method='SLSQP',
                  options={'xtol': 1e-8, 'disp': True})
  
       else:
           res = opt.minimize(self.func, coef, jac=self.func_deriv
                ,constraints=cons, method='SLSQP', options={'disp': True})
       return res
    

       """
       cons = ({'type': 'eq',
           'fun' : lambda x: np.array(coef[0]*self.x0**2 +coef[1].self.x0 + coef[0]-self.dfx0),
           'jac' : lambda x: np.array([self.x0**2, self.x0,1.0])})
    
       res = minimize(func, [-1.0,1.0], args=(-1.0,), jac=func_deriv,
                constraints=cons, method='SLSQP', options={'disp': True})



if __name__ == '__main__':
    from numpy.random import *
    x = np.arange(100)
    y = x**2
    ydata = y + 0.1*y*(random_sample(len(y))-0.5)
    #ydata = y
    # opt(x_vector,y_data,fx0=None,fxn=None,dfxn=None])
    c = opt(x,ydata,dfxn=100.0) 
   #c=opt(x,ydata,0.00,7000.0,0.0)
    #length of initial guess sets order of fit
    coef0 =[0.0,0.0,0.0,0.95,0.0,0.0]
    res=c.fit(coef0)
    print res
    print 'res.x'
    print res.x
    import matplotlib.pylab as plt

    
    #y0 =   coef0[0]*x**2 + coef0[1]*x + coef0[2]
    #yfit = res.x[0]*x**2 + res.x[1]*x + res.x[2]
    y0 = np.polyval(coef0,x)
    yfit = np.polyval(res.x,x)
    print 'slope= ',yfit[-1]-yfit[-2]
    plt.figure(1)
    plt.plot(x,y0,'c',x,ydata,'.k',x,yfit,'r')

    plt.show()
"""
