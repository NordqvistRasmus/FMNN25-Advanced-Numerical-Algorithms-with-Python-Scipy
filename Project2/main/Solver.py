
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 21:12:56 2019
@author: Mattias Lundström, Pontus Nordqvist, Johan Liljegren, Arvid Rolander, Antonio Alas
"""
from  scipy import *
from scipy.linalg import inv
from  pylab import *
from OptimizationProblem import OptimizationProblem
from scipy.optimize import minimize, rosen

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

    
class Solver:
    """
    Solver class which solves the problem.
    
    ...
    
    Attributes
    ----------
        problem: A problem which is defined by a problem class.
        function: The objective function defined by the problem.
        n: The dimension of the guess if it is given as a column vector.
        delta_grad: The small increment for calculating the finite difference
                    approximation for the gradient. Defaults to 1e-8.
        delta_values_grad: A diagonal matrix with the dimmension nxn whose 
                           elements are the delta_grad.
        delta_hess: The small increment for calculating the finite difference
                    approximation for the hessian. Defaults to 1e-4.
        delta_mat_hess: A diagonal matrix with the dimmension nxn whose 
                        elements are the delta_hess.
    """
    def __init__(self, problem):
        self.problem = problem
        self.function = problem.objective_function
        self.n = problem.x_0.shape[0]
        self.delta_grad = 1e-8
        self.delta_values_grad = diag([self.delta_grad for i in range(self.n)])
        self.delta_hess = 1e-4
        self.delta_mat_hess = diag([self.delta_hess for i in range(self.n)])
        
    
    #Plotting code partly from 'Three-Dimensional Plotting in Matplotlib' - Jake VanderPlas
    def plot(self, function, type = 'surface'): 
        
        if (type == 'surface'):
            x = linspace(-2,2,250)
            y = linspace(-1,3,250)
            X, Y = meshgrid(x, y)
            Z = rosen([X, Y])
            ax = plt.axes(projection='3d')
            ax.plot_surface(X,Y,Z, norm = LogNorm(), rstride = 5, cstride = 5, cmap = 'RdGy_r', alpha = .9, edgecolor = 'none')
            ax.set_title('Rosenbrock function - Surface plot')
            ax.set_xlabel('x_1')
            ax.set_ylabel('x_2')
            ax.set_zlabel('f(x)')

        if (type == 'contour'):
            ax = plt.axes()
            x = linspace(-5,5,500)
            y = linspace(-5,5,500)
            X, Y = meshgrid(x, y)
            Z = rosen([X, Y])
            plt.contour(X, Y, Z, 150, cmap = 'RdGy');
            ax.set_title('Rosenbrock function - Contour plot')
            ax.set_xlabel('x_1')
            ax.set_ylabel('x_2')
      
        plt.show()
    
    def newton(self, mode='default', tol = 1e-8, maxIteration = 1000):
        iterations = 0
        x_k = self.problem.x_0
        x_next = self._newton_step(x_k, mode)
        
        while(norm(x_k-x_next)>tol
              and norm(self._gradient(x_k)) > tol
              and iterations < maxIteration):
            x_k = x_next
            x_next = self._newton_step(x_k, mode)
            iterations +=1
        if norm(x_k-x_next) > tol:
            #Might lead to wrong answer?
            pass
        if maxIteration == iterations:
            print('Mode:', mode, "------ Did not converge in" ,iterations, "iterations. \n")
            return x_next
            
        print('Mode:', mode, "------ Converged in" ,iterations, "iterations. \n")
        #print('Norm of gradient:', norm(self._gradient(x_k)))
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
        """
        Added because we need to opt have gradient //M
        """
        if self.problem.gradient is None:
            gradient = array([(self.function(x_k + norm(x_k)*self.delta_values_grad[:,i])
                            -self.function(x_k - norm(x_k)*self.delta_values_grad[:,i]))
                            /(2*norm(x_k)*self.delta_grad)
                            for i in range(self.n)])
            return gradient
        else:
            return self.problem.gradient(x_k) #If we have gradient from problem (as a function)
        """
        Debug print
        """
        #print(norm(gradient))
        #return gradient
    
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
        
class QuasiNewtonSolver(Solver):
    
    
    def __init__(self, problem):
        super().__init__(problem)
        self.hessian = super()._hessian(problem.x_0)
        self.inverse_hessian = inv(super()._hessian(problem.x_0))
        
    def _inexact_line_search(self, x_k, s_k, mode, a_l=0, a_u=1e5,
                             rho=0.1, sigma=0.7, tau=0.1, chi=9):
        a_0 = (a_u-a_l)/2
        d_a = self.delta_grad
        f_alpha = lambda alpha: self.function(x_k + alpha*s_k)
        f_prime = lambda alpha: (self.function(x_k + (alpha+d_a)*s_k)
                                - self.function(x_k +(alpha)*s_k))/(d_a)
        fp_a0 = f_prime(a_0)
        fp_al = f_prime(a_l)
        f_a0 = f_alpha(a_0)
        f_al = f_alpha(a_l)
        
        if mode == 'wolfe':
            LC = fp_a0 >= sigma*fp_al
            RC = f_a0 <= f_al + rho*(a_0 - a_l)*fp_al
        if mode == 'goldstein':
            LC = f_a0 >= f_al + (1-rho)*(a_0-a_l)*fp_al
            RC = f_al <= f_al + rho*(a_0-a_l)*fp_al
        
        while not (LC and RC):
            if not (LC):
                da_0 = (a_0 - a_l)*(fp_a0/(fp_al - fp_a0))
                da_0 = max(da_0, tau*(a_0 - a_l))
                da_0 = min(da_0, chi*(a_0 - a_l))
                a_l = a_0
                a_0 += da_0
            else:
                a_u = min(a_0, a_u)
                ia_0 = (((a_0 - a_l)**2)*fp_al)/(2*(f_al - f_a0 + (a_0 - a_l)*fp_al))
                ia_0 = max(ia_0, a_l + tau*(a_u - a_l))
                ia_0 = min(ia_0, a_u - tau*(a_u - a_l))
                a_0 = ia_0
            
            fp_a0 = f_prime(a_0)
            fp_al = f_prime(a_l)
            f_a0 = f_alpha(a_0)
            f_al = f_alpha(a_l)
            if mode == 'wolfe':
                LC = fp_a0 >= sigma*fp_al
                RC = f_a0 <= f_al + rho*(a_0 - a_l)*fp_al
            if mode == 'goldstein':
                LC = f_a0 >= f_al + (1-rho)*(a_0-a_l)*fp_al
                RC = f_al <= f_al + rho*(a_0-a_l)*fp_al
            
        return a_0, f_a0
                
    
    def _update_hessian(self, x_k, x_next, alpha):
        pass
     
    def _newton_step(self, x_k, mode):
        s_k = (-1)*self.inverse_hessian@self._gradient(x_k)
        s_k = s_k/norm(s_k)
        alpha, f_alpha = self._inexact_line_search(x_k, s_k, 'wolfe')
        x_next = x_k + alpha*s_k
        self._update_hessian(x_k, alpha)
        return x_next
        
    
class GoodBroydenSolver(QuasiNewtonSolver):
    
    
    def _update_hessian(self, x_k, alpha):
        delta =  alpha*(self.inverse_hessian@self._gradient(x_k)) #osäker på denna
        #delta = x_next-x_k
        print('x_next: ', x_next,'x_k: ', x_k, 'alpha: ', alpha)
        gamma = self._gradient(x_next) - self._gradient(x_k)
        print('gammahamma',gamma)
        u = delta - self.inverse_hessian@gamma
        a = 1 /u@gamma
        self.inverse_hessian = self.inverse_hessian + a*outer(u, u)

#Allmänt osäker på denna bad Broyden-metoden, svårt att hitta info.       
class BadBroydenSolver(QuasiNewtonSolver):
    
    def _update_hessian(self, x_k, alpha):
        delta =  alpha*(self.inverse_hessian@self._gradient(x_k))
        gamma = self._gradient(x_k) - self._gradient(x_k-delta)
        u = gamma - self.inverse_hessian@delta 
        a = 1 / gamma@gamma 
        self.inverse_hessian = self.inverse_hessian + a*outer(u,gamma)
        
class DFP2Solver(QuasiNewtonSolver):
    
    def _update_hessian(self, x_k, alpha):
        delta =  alpha*(self.inverse_hessian@self._gradient(x_k))
        gamma = self._gradient(x_k) - self._gradient(x_k-delta)
        u1 = outer(delta,delta)
        a1 = 1/delta@gamma
        u2 = outer(self.inverse_hessian@gamma,gamma)@self.inverse_hessian
        a2 = 1/gamma@(self.inverse_hessian@gamma)
        self.inverse_hessian = self.inverse_hessian + a1*u1 - a2*u2
    
class BFGS2Solver(QuasiNewtonSolver):
    
    def _update_hessian(self, x_k, alpha):
        delta =  alpha*(self.inverse_hessian@self._gradient(x_k))
        gamma = self._gradient(x_k) - self._gradient(x_k-delta)
        hg = self.inverse_hessian@gamma
        dg = delta@gamma
        u1 = gamma@hg
        a1 = a2 = a3 = 1/dg
        u2 = outer(delta,delta)
        u3 = outer(delta, gamma)@self.inverse_hessian + self.inverse_hessian@outer(gamma, delta)
        #u3 = outer(dg, self.inverse_hessian) + outer(dg, self.inverse_hessian).T #Transponat för motsat ordning 
        self.inverse_hessian = self.inverse_hessian+(1+a1*u1)*(a2*u2)-a3*u3
    
    
if __name__ == '__main__':
    #function = lambda x: (x[0]-1)**2 + x[1]**2
    function = lambda x: 100*((x[1]-x[0]**2)**2)+(1-x[0])**2
    op = OptimizationProblem(function, array([5,5]))
    s = Solver(op)
    s.plot(rosen, 'surface')
    #s.plot(rosen, 'contour')
    
    #GoodBoy = GoodBroydenSolver(op)
    #BadBoy = BadBroydenSolver(op)
    #DP2 = DFP2Solver(op)
    #BFGS = BFGS2Solver(op)
    zero = s.newton(mode = 'default')
    zero1 = s.newton(mode='exact')
    #zero2 = GoodBoy.newton()
    #zero3 = BadBoy.newton()
    #zero4 = DP2.newton()
    #zero5 = BFGS.newton()
    
    print('Regular newton gives: ',zero, '\n')
    print('Newton with exact line search gives: ',zero1, '\n')
    #print('Good Broyden gives: ',zero2, '\n')
    #print('Bad Broyden gives: ',zero3, '\n')
    #print('DFP2 gives: ',zero4, '\n')    
    #print('BFGS2 gives: ',zero5, '\n')
    
