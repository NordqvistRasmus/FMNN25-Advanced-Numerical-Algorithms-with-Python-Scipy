
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

    def newton(self, mode='default', tol = 1e-11, maxIteration = 1000):
        iterations = 0
        x_k = self.problem.x_0
        x_next = self._newton_step(x_k,mode)
        
        while(norm(x_k-x_next)>tol
              and norm(self._gradient(x_k)) > tol
              and iterations < maxIteration):
            x_k = x_next
            x_next = self._newton_step(x_k, mode)
            iterations +=1
        if norm(x_k-x_next) > tol:
            pass
            #raise exception didn't converge 
        print("I converged in " ,iterations, " iterations")
        return x_next
    
    def _newton_step(self, x_k, mode):
        
        if self.problem.gradient is not None:
            gradient = self.problem.gradient(x_k)
        else:
            gradient = self._gradient(x_k)
        
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
    
    def exact_line_search(self, x_k, gradient, hessian, tol=1e-6, alpha_0=1):
        
        delta_grad = self.delta_grad
        delta_hess = self.delta_hess
        hessian = inv(hessian)
        s = -hessian@gradient
        f_alpha = lambda alpha: self.function(x_k + alpha*s)
        deriv = lambda alpha: (f_alpha(alpha+delta_grad)-f_alpha(alpha-delta_grad))/(delta_grad*2)
        sec_deriv = lambda alpha: (f_alpha(alpha+delta_hess) - 2*f_alpha(alpha) +
                                   f_alpha(alpha-delta_hess))/delta_hess**2
        if deriv(alpha_0) < tol:
            return alpha_0
        alpha_k = alpha_0
        alpha_next = alpha_k - deriv(alpha_k)/sec_deriv(alpha_k)
        
        while abs(alpha_next-alpha_k) > tol and abs(deriv(alpha_k)) > tol:
            alpha_k = alpha_next
            sec_derivative = sec_deriv(alpha_k)
            
            if sec_derivative == 0:
                return alpha_k
            else:
                alpha_next = alpha_k - deriv(alpha_k)/sec_derivative
        return alpha_next
    
    def inexact_line_search(self):
        pass
    
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
        
    
# här använder jag mig av en faktor delta som är = alpha * s, där s ska vara s = -hessian(x_k) @ g(x_k)
# Vet inte om detta skall deklareras någonstans i newton_step-metoden som ett attribut eller
# om vi ska ha det som parametrar i varje _hessian skuggning. Lämnar detta öpper för er att bestämma
# mitt förslag är att vi lägger till det som parametrar.
# I övrigt är alla metoder implementerade.
class GoodBroydenSolver(Solver):
    
    def _hessian(self, x_k):
        delta =  self.alpha*(-self.hessian@self.gradient) #osäker på denna
        gamma = self.gradient(x_k) - self.gradient(x_k-delta)
        u = delta - self.hessian@gamma
        a = 1 /u@gamma
        return self.hessian + a*outer(u, u)

#Allmänt osäker på denna bad Broyden-metoden, svårt att hitta info.       
class BadBroydenSolver(Solver):
    
    def _hessian(self, x_k):
        delta =  self.alpha*(-self.hessian@self.gradient)
        gamma = self.gradient(x_k) - self.gradient(x_k-delta)
        u = gamma - self.hessian@delta 
        a = 1 / gamma@gamma #Inte säker på om det ska vara gamma@gamma eller
        # delta@delta
        return self.hessian + a*outer(u,gamma)
        
class DFP2Solver(Solver):
    
    def _hessian(self, x_k):
        delta =  self.alpha*(-self.hessian@self.gradient)
        gamma = self.gradient(x_k) - self.gradient(x_k-delta)
        u1 = outer(delta,delta)
        a1 = delta@gamma
        u2 = outer(self.hessian@gamma,gamma)@self.hessian
        a2 = gamma@(self.hessian@gamma)
        return self.hessian + a1*u1 - a2*u2
    
class BFGS2Solver(Solver):
    
    def _hessian(self, x_k):
        delta =  self.alpha*(-self.hessian@self.gradient)
        gamma = self.gradient(x_k) - self.gradient(x_k-delta)
        hg = self.hessian@gamma
        dg = delta@gamma
        u1 = gamma@(self.hessian@gamma)
        a1 = a2 = a3 = dg
        u2 = outer(delta,delta)
        u3 = outer(hg,delta) + outer(hg,delta).T #Transponat för motsat ordning 
        return self.hessian+(1+a1*u1)*(a2*u2)-a3*u3
    
    
if __name__ == '__main__':
    #function = lambda x: (x[0]-1)**2 + x[1]**2
    function = lambda x: 100*((x[1]-x[0]**2)**2)+(1-x[0])**2
    op = OptimizationProblem(function, array([230,30]))
    s = Solver(op)
    zero = s.newton(mode='exact')
    print(zero)
    
        
    
    
