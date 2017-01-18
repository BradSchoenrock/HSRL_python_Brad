import numpy as np


def cubic_fit(bin0,bin1,v0,v1,s0,s1):
    """cubic_fit(bin_vec,bin0,bin1,v0,v1,s0,s1)
       returns coeficents of a cubic fit specified at bin_vec locations
       where the fit has value and derivative v0, s0 at bin0 and v1, s1 at bin1.
       a + b*bin + c*bin^2 + d*bin^3 are returned as a vector [d,c,v,a]"""

   
   
    in_array = np.array([[bin0**3,bin0**2,bin0,1.0]
                        ,[bin1**3,bin1**2,bin1,1.0]
                        ,[3.0*bin0**2,2.0*bin0,1.0,0.0]
                        ,[3.0*bin1**2,2.0*bin1,1.0,0.0]])
    
    in_vector = np.array([v0,v1,s0,s1])
    z = np.linalg.solve(in_array,in_vector)
    
    return z

   
    
if __name__ == '__main__':
        #cubic_fit(bin_vec,bin0,bin1,v0,v1,s0,s1)
      
        print
        print
        bin_vec = np.arange(11)
        print bin_vec
        import matplotlib.pylab as plt
        
        a = 0.0
        b = 1.0
        c = 0.1
        d = 10**-3
        print 'a=',a,' b=',b,' c=',c,' d=',d
        import matplotlib.pylab as plt
        plt.figure(1)
        plt.plot(bin_vec,.001+a + b*bin_vec + c*bin_vec**2 + d*bin_vec**3,'b')


        #inputs computed from polynomial
        bin0 = 0
        bin1 = 10
        x0= a + b*bin0 + c*bin0**2 + d*bin0**3
        x1= a +b*bin1 + c*bin1**2 + d*bin1**3
        s0= b + 2.0*c*bin0 + 3.0*d*bin0**2
        s1= b + 2.0*c*bin1 + 3.0*d*bin1**2
           
        print 'value at  0= ' ,x0
        print 'value at 10= ' ,x1
        print 'slope at  0= ' ,s0
        print 'slope at 10= ' ,s1
        print       
        a,b,c,d = cubic_fit(bin_vec,bin0,bin1, x0,x1,s0,s1)
        print 'a= ',a
        print 'b= ',b
        print 'c= ',c
        print 'd= ',d
     
       
        plt.plot(bin_vec,a + b*bin_vec + c*bin_vec**2 + d*bin_vec**3,'r')         
        plt.show()
