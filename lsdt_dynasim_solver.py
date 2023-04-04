#!/usr/bin/env python
# coding: utf-8

"""
Stand 16.12.2022
"""
from scipy.interpolate import interp1d
import numpy as np
from matplotlib import pyplot as plt
import csv
from re import sub
import math
   
#########################################################################
##### DGL-Solver ########################################################
#########################################################################

class ODESolver:
    """ODESolver
    Solves ODE on the form:
    u' = f(u, t), u(0) = U0
    Parameters
    ----------
    f : callable
        Right-hand-side function f(u, t)
    """

    def __init__(self, f, KNOTEN, Klasse):               
        # DYDANSIM
        if (Klasse == "class_DGL_Masse_trans_x_Masse"):
            self.Knoteni = KNOTEN[0]            
            self.KRAFT = [0,0]        
        if (Klasse == "class_DGL_Masse_trans_y_Masse"):
            self.Knoteni = KNOTEN[0]            
            self.KRAFT = [0,0]        
        if (Klasse == "class_DGL_Masse_trans_z_Masse"):
            self.Knoteni = KNOTEN[0]            
            self.KRAFT = [0,0]  
        if (Klasse == "class_DGL_Rolle_rot_y_Masse"):
            self.Knoteni = KNOTEN[0]            
            self.Knotenj = KNOTEN[1]            
            self.Knotenk = KNOTEN[2]            
            self.KRAFT = [0,0,0]              
        if (Klasse == "class_DGL_gelagerteRolle_rot_x_Masse"):
            self.Knoteni = KNOTEN[0]            
            self.KRAFT = [0,0]
        if (Klasse == "class_DGL_Fadenpendel_rot_x_Masse"):
            self.Knoteni = KNOTEN[0]            
            self.KRAFT = [0,0]
        if (Klasse == "class_DGL_Fadenpendel_mit_Seilkraft_rot_x_Masse"):
            self.Knoteni = KNOTEN[0]            
            self.KRAFT = [0,0]
            
        self.f = f
        self.number_of_eqns = 2
        U0in = [0,0]
        U0in = np.asarray(U0in)
        self.number_of_eqns = U0in.size                        
        self.U0 = U0in                     
        
    def aktualisiere_initial_conditions(self, Klasse):            
        #DYNDASIM
        if (Klasse == "class_DGL_Masse_trans_x_Masse"):
            self.KRAFT = [self.Knoteni.F_x]
            self.U0 = [0, self.Knoteni.x_p]                    
        if (Klasse == "class_DGL_Masse_trans_y_Masse"):
            self.KRAFT = [self.Knoteni.F_y]
            self.U0 = [0, self.Knoteni.y_p]        
        if (Klasse == "class_DGL_Masse_trans_z_Masse"):
            self.KRAFT = [self.Knoteni.F_z]
            self.U0 = [0, self.Knoteni.z_p]        
        if (Klasse == "class_DGL_Rolle_rot_y_Masse"):
            self.KRAFT = [self.Knoteni.F_x,self.Knotenj.F_x,self.Knotenk.F_x]
            self.U0 = [0, self.Knoteni.roty_p]        
        if (Klasse == "class_DGL_gelagerteRolle_rot_x_Masse"):
            self.KRAFT = [self.Knoteni.M_x]
            self.U0 = [0, self.Knoteni.rotx_p]        
        if (Klasse == "class_DGL_Fadenpendel_rot_x_Masse"):
            self.KRAFT = [self.Knoteni.M_x]
            self.U0 = [0, self.Knoteni.rotx_p]        
        if (Klasse == "class_DGL_Fadenpendel_mit_Seilkraft_rot_x_Masse"):
            self.KRAFT = [self.Knoteni.M_x]
            self.U0 = [0, self.Knoteni.rotx_p]                    
            
            
    def solve_timestep(self, delta_t):
        """
        Solves ODE according to given time points.
        The resolution is implied by spacing of 
        time points.
        Parameters
        ----------
        time_points : array_like
            Time points to solve for
        
        Returns
        -------
        u : array_like
            Solution
        t : array_like
            Time points corresponding to solution
        """

        time_points = [0, delta_t]
        self.t = np.asarray(time_points)
        n = self.t.size

        self.u = np.zeros((n, self.number_of_eqns))        
        
        self.u[0, :] = self.U0                
        # Integrate
        for i in range(n - 1):
            self.i = i
            self.u[i + 1] = self.advance()

        return self.u 

        def advance(self):
            """Advance solution one time step."""
            raise NotImplementedError            
        
class RungeKutta4_Masse(ODESolver):
    def advance(self):
        u, f, i, t, KRAFT = self.u, self.f, self.i, self.t, self.KRAFT
        dt = t[i + 1] - t[i]
        dt2 = dt / 2
        K1 = dt * f(u[i, :], KRAFT, t[i])
        K2 = dt * f(u[i, :] + 0.5 * K1, KRAFT, t[i] + dt2)
        K3 = dt * f(u[i, :] + 0.5 * K2, KRAFT, t[i] + dt2)
        K4 = dt * f(u[i, :] + K3, KRAFT, t[i] + dt)
        return u[i, :] + (1 / 6) * (K1 + 2 * K2 + 2 * K3 + K4)   
        
class Euler_Masse(ODESolver):
    def advance(self):
        u, f, i, t, KRAFT = self.u, self.f, self.i, self.t, self.KRAFT
        dt = t[i + 1] - t[i]
        dt2 = dt / 2
        K1 = dt * f(u[i, :], KRAFT, t[i])        
        return u[i, :] + K1   
 


