import numpy as np


def ninth_order_fit(bin0,bin1,bin2,bin3,bin4,v0,v1,v2,v3,v4,s0,s1,s2,s3,s4):
    """seventh_order_fit(bin0,bin1,bin2,bin3,v0,v1,v2,v3,s0,s1,s2,s3)
       returns coeficents of a 9-th order polynomial
       where the fit has values v0, v1, v2, v3 and v4 at locations bin0, bin1, bin2, bin3 and bin4
       and derivatives s0, s1, s2, s3 and s4 at bin0,bin1, bin2, bin3 and bin4.
       p9*bin^9 + ........ +p2*bin^2 + p1*bin + p0 are returned as a vector [p9....p2,p1,p0]"""

    in_array =np.array( \
        [[bin0**9,bin0**8,bin0**7,bin0**6,bin0**5,bin0**4,bin0**3,bin0**2,bin0,1.0]
        ,[bin1**9,bin1**8,bin1**7,bin1**6,bin1**5,bin1**4,bin1**3,bin1**2,bin1,1.0]
        ,[bin2**9,bin2**8,bin2**7,bin2**6,bin2**5,bin2**4,bin2**3,bin2**2,bin2,1.0]
        ,[bin3**9,bin3**8,bin3**7,bin3**6,bin3**5,bin3**4,bin3**3,bin3**2,bin3,1.0]
        ,[bin4**9,bin4**8,bin4**7,bin4**6,bin4**5,bin4**4,bin4**3,bin4**2,bin4,1.0]
        ,[9.0*bin0**8,8.0*bin0**7,7.0*bin0**6,6.0*bin0**5,5.0*bin0**4,4.0*bin0**3,3.0*bin0**2,2.0*bin0,1.0,0.0]
        ,[9.0*bin1**8,8.0*bin1**7,7.0*bin1**6,6.0*bin1**5,5.0*bin1**4,4.0*bin1**3,3.0*bin1**2,2.0*bin1,1.0,0.0] 
        ,[9.0*bin2**8,8.0*bin2**7,7.0*bin2**6,6.0*bin2**5,5.0*bin2**4,4.0*bin2**3,3.0*bin2**2,2.0*bin2,1.0,0.0]
        ,[9.0*bin3**8,8.0*bin3**7,7.0*bin3**6,6.0*bin3**5,5.0*bin3**4,4.0*bin3**3,3.0*bin3**2,2.0*bin3,1.0,0.0]
        ,[9.0*bin4**8,8.0*bin4**7,7.0*bin4**6,6.0*bin4**5,5.0*bin4**4,4.0*bin4**3,3.0*bin4**2,2.0*bin4,1.0,0.0]])
    in_vector = np.array([v0,v1,v2,v3,v4,s0,s1,s2,s3,s4])
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
