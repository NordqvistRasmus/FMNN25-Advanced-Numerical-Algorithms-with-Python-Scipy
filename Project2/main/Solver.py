
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 21:12:56 2019
@author: Mattias Lundström1
"""
from  scipy import *
from scipy.linalg import inv
from  pylab import *
from OptimizationProblem import OptimizationProblem
    
class Solver:
    
    def __init__(self, problem):
        self.problem = problem
        self.function = problem.objective_function
        self.n = problem.x_0.shape[0]
        self.delta_grad = 1e-8
        self.delta_values_grad = diag([self.delta_grad for i in range(self.n)])
        self.delta_hess = 1e-4
        self.delta_mat_hess = diag([self.delta_hess for i in range(self.n)])
        
    
    def __call__(self):
        pass

    def newton(self, mode='default', tol = 1e-12, maxIteration = 99):
        iterations = 0
        x_k = self.problem.x_0
        x_next = self._newton_step(x_k,mode)
        
        while(norm(x_k-x_next)>tol and iterations < maxIteration):
            x_k = x_next
            x_next = self._newton_step(x_k, mode)
            iterations +=1
        if norm(x_k-x_next) > tol:
            pass
            #raise exception didn't converge 
        
        return x_next
    
    def _newton_step(self, x_k, mode):
        
        if self.problem.gradient is not None:
            gradient = self.problem.gradient(x_k)
        else:
            gradient = self._gradient(x_k)

        #Detta kommer behövas göras om, Hessian kan inte ligga här då den behöver 
        #alpha.
        hessian = self._hessian(x_k)
        
        if mode == 'default':
            alpha = 1
        
        elif mode == 'exact':
            alpha = self.exact_line_search(x_k, gradient, hessian)
            
        elif mode == 'inexact':
            alpha = self.inexact_line_search(x_k)
        
        x_next = solve(hessian, hessian@x_k - alpha*gradient)
        
        return x_next
    
    def _gradient(self, x_k):
        gradient = array([(self.function(x_k + self.delta_values_grad[:,i])
                            -self.function(x_k - self.delta_values_grad[:,i]))/(2*self.delta_grad)
                            for i in range(self.n)])
        print(norm(gradient))
        return gradient
    
    def _search_dir(self, x_k):
        pass    
    
    def exact_line_search(self, x_k, gradient, hessian, tol=1e-8, alpha_0=20):
        
        delta_grad = self.delta_grad
        delta_hess = self.delta_hess
        hessian = inv(hessian)
        s = -hessian@gradient
        f_alpha = lambda alpha: self.function(x_k + alpha*s)
        deriv = lambda alpha: (f_alpha(alpha+delta_grad)-f_alpha(alpha-delta_grad))/(delta_grad*2)
        sec_deriv = lambda alpha: (f_alpha(alpha+delta_hess) - 2*f_alpha(alpha) +
                                   f_alpha(alpha-delta_hess))/delta_hess**2
        alpha_k = alpha_0
        alpha_next = alpha_next = alpha_k - deriv(alpha_k)/sec_deriv(alpha_k)
        while abs(alpha_next-alpha_k) > tol:
            alpha_k = alpha_next
            alpha_next = alpha_k - deriv(alpha_k)/sec_deriv(alpha_k)
        return alpha_next
    
    #Some what of a skeleton    
    def inexact_line_search(self, x_k, gradient, hessian, tol = 1e-8, alpha_0 = 20, 
                            rho = 0.1, sigma = 0.7, tau = 0.1, chi = 9):

        #Skulle kanske kunna vara hjälpfunktioner
        if (condition):
            lc = True
        else:
            lc = False
            
        if (condition):
            rc = True
        else:
            rc = False
        
        while not (lc and rc):
            if not lc:
                #block1
                #extrapolation
            else:
                #block2
                #interpolation
            lc = ...
            rc = ...
        return alpha_0
    
    def _hessian(self, x_k):
        hessian = array([[self._second_part_div(x_k, i, j) for j in range(self.n)]
                          for i in range(self.n)])
        hessian = 0.5*(hessian+hessian.T)
        #print(hessian)
        return hessian
        
    def _second_part_div(self, x_k, i, j):
        div = (self.function(x_k+self.delta_mat_hess[:,i] + self.delta_mat_hess[:,j]) -  
               self.function(x_k+self.delta_mat_hess[:,i])-
               self.function(x_k+self.delta_mat_hess[:,j])+
               self.function(x_k))/self.delta_hess**2            # consider changing delta to something bigger
        return div
        
    
class GoodBroydenSolver(QuasiNewton):
    
    def _hessian(self, x_k, alpha, hessian, gradient):
        delta =  alpha*(hessian@gradient)
        gamma = gradient(x_k) - gradient(x_k-delta)
        u = delta - hessian@gamma
        a = 1 /u@gamma
        return hessian + a*outer(u, u)
     
class BadBroydenSolver(QuasiNewton):
    
    def _hessian(self, x_k, alpha, hessian, gradient):
        delta =  alpha*(hessian@gradient)
        gamma = gradient(x_k) - gradient(x_k-delta)
        u = gamma - hessian@delta 
        a = 1 / gamma@gamma
        return hessian + a*outer(u,gamma)
        
class DFP2Solver(QuasiNewton):
    
    def _hessian(self, x_k, alpha, hessian, gradient):
        delta =  alpha*(hessian@gradient)
        gamma = gradient(x_k) - gradient(x_k-delta)
        hg = hessian@gamma
        u1 = outer(delta,delta)
        a1 = delta@gamma
        u2 = outer(hg,gamma)@hessian
        a2 = gamma@hg
        return hessian + a1*u1 - a2*u2
    
class BFGS2Solver(QuasiNewton):
    
    def _hessian(self, x_k, alpha, hessian, gradient):
        delta =  alpha*(hessian@gradient)
        gamma = gradient(x_k) - gradient(x_k-delta)
        hg = hessian@gamma
        dg = delta@gamma
        u1 = gamma@hg
        a1 = a2 = a3 = dg
        u2 = outer(delta,delta)
        u3 = outer(dg,hessian) + outer(dg,delta).T #Transponat för motsat ordning 
        return hessian+(1+a1*u1)*(a2*u2)-a3*u3
    
    
if __name__ == '__main__':
    #function = lambda x: (x[0]-1)**2 + x[1]**2
    function = lambda x: 100*((x[1]-x[0]**2)**2)+(1-x[0])**2
    op = OptimizationProblem(function, array([2,2]))
    s = Solver(op)
    zero = s.newton(mode='exact')
    print(zero)
    
        
    
    
