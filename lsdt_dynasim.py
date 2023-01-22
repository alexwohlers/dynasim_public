#!/usr/bin/env python
# coding: utf-8

"""
Stand 16.12.2022
"""

import numpy as np
from matplotlib import pyplot as plt
import csv
from re import sub
import math


#########################################################################
##### Einstellungen für Plots ###########################################
#########################################################################

plt.rcParams['figure.dpi']  = 300
plt.rcParams['savefig.dpi'] = 300

#########################################################################
##### Definition der Knoten #############################################
#########################################################################

class class_Knoten:
    def __init__(self):
        self.x                 = 0
        self.x_p               = 0
        self.x_pp              = 0    
        self.y                 = 0
        self.y_p               = 0
        self.y_pp              = 0    
        self.z                 = 0
        self.z_p               = 0
        self.z_pp              = 0    
        self.rotx              = 0
        self.rotx_p            = 0
        self.rotx_pp           = 0    
        self.roty              = 0
        self.roty_p            = 0
        self.roty_pp           = 0    
        self.rotz              = 0
        self.rotz_p            = 0
        self.rotz_pp           = 0    
        self.F_x               = 0     
        self.F_y               = 0     
        self.F_z               = 0     
        self.M_x               = 0     
        self.M_y               = 0     
        self.M_z               = 0
    def Reset_Knotenkraft(self):
        self.F_x               = 0
        self.F_y               = 0
        self.F_z               = 0
        self.M_x               = 0
        self.M_y               = 0
        self.M_z               = 0
    def Reset_Knoten(self):
        self.x                 = 0
        self.x_p               = 0
        self.x_pp              = 0    
        self.y                 = 0
        self.y_p               = 0
        self.y_pp              = 0    
        self.z                 = 0
        self.z_p               = 0
        self.z_pp              = 0    
        self.rotx              = 0
        self.rotx_p            = 0
        self.rotx_pp           = 0    
        self.roty              = 0
        self.roty_p            = 0
        self.roty_pp           = 0    
        self.rotz              = 0
        self.rotz_p            = 0
        self.rotz_pp           = 0
        self.F_x               = 0
        self.F_y               = 0
        self.F_z               = 0
        self.M_x               = 0
        self.M_y               = 0
        self.M_z               = 0
       
class class_Datenlogger_Knoten:
    def __init__(self):
        self.x                 = []
        self.x_p               = []
        self.x_pp              = []    
        self.y                 = []
        self.y_p               = []
        self.y_pp              = []    
        self.z                 = []
        self.z_p               = []
        self.z_pp              = []    
        self.rotx              = []
        self.rotx_p            = [] 
        self.rotx_pp           = []    
        self.roty              = []
        self.roty_p            = []
        self.roty_pp           = []    
        self.rotz              = []
        self.rotz_p            = []
        self.rotz_pp           = []    
        self.F_x               = []     
        self.F_y               = []     
        self.F_z               = []     
        self.M_x               = []     
        self.M_y               = []     
        self.M_z               = []
    def Werte_loeschen(self):
        self.x.clear()
        self.x_p.clear()
        self.x_pp.clear()
        self.y.clear()
        self.y_p.clear()
        self.y_pp.clear()
        self.z.clear()
        self.z_p.clear()
        self.z_pp.clear()
        self.rotx.clear()
        self.rotx_p.clear()
        self.rotx_pp.clear()
        self.roty.clear()
        self.roty_p.clear()
        self.roty_pp.clear()
        self.rotz.clear()
        self.rotz_p.clear()
        self.rotz_pp.clear()
        self.F_x.clear()
        self.F_y.clear()
        self.F_z.clear()
        self.M_x.clear()
        self.M_y.clear()
        self.M_z.clear()
    def Werte_anhaengen(self, Knoten):
        self.x.append(Knoten.x)
        self.x_p.append(Knoten.x_p)
        self.x_pp.append(Knoten.x_pp)
        self.y.append(Knoten.y)
        self.y_p.append(Knoten.y_p)
        self.y_pp.append(Knoten.y_pp)
        self.z.append(Knoten.z)
        self.z_p.append(Knoten.z_p)
        self.z_pp.append(Knoten.z_pp)
        self.rotx.append(Knoten.rotx)
        self.rotx_p.append(Knoten.rotx_p)
        self.rotx_pp.append(Knoten.rotx_pp)
        self.roty.append(Knoten.roty)
        self.roty_p.append(Knoten.roty_p)
        self.roty_pp.append(Knoten.roty_pp)
        self.rotz.append(Knoten.rotz)
        self.rotz_p.append(Knoten.rotz_p)
        self.rotz_pp.append(Knoten.rotz_pp)
        self.F_x.append(Knoten.F_x)
        self.F_y.append(Knoten.F_y)
        self.F_z.append(Knoten.F_z)
        self.M_x.append(Knoten.M_x)
        self.M_y.append(Knoten.M_y)
        self.M_z.append(Knoten.M_z)        

class class_Datenlogger_Element:
    def __init__(self):
        self.E_kin             = []
        self.E_pot             = []
        self.E_diss            = [] 
        self.E_kin_p           = []
        self.E_pot_p           = []
        self.E_diss_p          = []
    def Werte_loeschen(self):
        self.E_kin.clear()
        self.E_pot.clear()
        self.E_diss.clear()
        self.E_kin_p.clear()
        self.E_pot_p.clear()
        self.E_diss_p.clear()
    def Werte_anhaengen(self, Element):
        self.E_kin.append(Element.E_kin)
        self.E_pot.append(Element.E_pot)
        self.E_diss.append(Element.E_diss)            
        self.E_kin_p.append(Element.E_kin_p)
        self.E_pot_p.append(Element.E_pot_p)
        self.E_diss_p.append(Element.E_diss_p)            
        
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
        if (Klasse == "class_DGL_Masse_trans_x"):
            self.Knoteni = KNOTEN[0]            
            self.KRAFT = [0,0]        
        if (Klasse == "class_DGL_Masse_trans_y"):
            self.Knoteni = KNOTEN[0]            
            self.KRAFT = [0,0]        
        if (Klasse == "class_DGL_Masse_trans_z"):
            self.Knoteni = KNOTEN[0]            
            self.KRAFT = [0,0]  
        if (Klasse == "class_DGL_gelagerteRolle_rot_x"):
            self.Knoteni = KNOTEN[0]            
            self.KRAFT = [0,0]
        if (Klasse == "class_DGL_Federpendel_rot_x"):
            self.Knoteni = KNOTEN[0]            
            self.KRAFT = [0,0]
        if (Klasse == "class_DGL_Federpendel_mit_Seilkraft_rot_x"):
            self.Knoteni = KNOTEN[0]            
            self.KRAFT = [0,0]

        self.f = f
        self.number_of_eqns = 2
        U0in = [0,0]
        U0in = np.asarray(U0in)
        self.number_of_eqns = U0in.size                        
        self.U0 = U0in                     
        
    def aktualisiere_initial_conditions(self, Klasse):                
        if (Klasse == "class_DGL_Masse_trans_x"):
            self.KRAFT = [self.Knoteni.F_x]
            self.U0 = [0, self.Knoteni.x_p]                    
        if (Klasse == "class_DGL_Masse_trans_y"):
            self.KRAFT = [self.Knoteni.F_y]
            self.U0 = [0, self.Knoteni.y_p]        
        if (Klasse == "class_DGL_Masse_trans_z"):
            self.KRAFT = [self.Knoteni.F_z]
            self.U0 = [0, self.Knoteni.z_p]        
        if (Klasse == "class_DGL_gelagerteRolle_rot_x"):
            self.KRAFT = [self.Knoteni.M_x]
            self.U0 = [0, self.Knoteni.rotx_p]        
        if (Klasse == "class_DGL_Federpendel_rot_x"):
            self.KRAFT = [self.Knoteni.M_x]
            self.U0 = [0, self.Knoteni.rotx_p]        
        if (Klasse == "class_DGL_Federpendel_mit_Seilkraft_rot_x"):
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
            
class RungeKutta4(ODESolver):
    def advance(self):
        u, f, i, t, KRAFT = self.u, self.f, self.i, self.t, self.KRAFT
        dt = t[i + 1] - t[i]
        dt2 = dt / 2
        K1 = dt * f(u[i, :], KRAFT, t[i])
        K2 = dt * f(u[i, :] + 0.5 * K1, KRAFT, t[i] + dt2)
        K3 = dt * f(u[i, :] + 0.5 * K2, KRAFT, t[i] + dt2)
        K4 = dt * f(u[i, :] + K3, KRAFT, t[i] + dt)
        return u[i, :] + (1 / 6) * (K1 + 2 * K2 + 2 * K3 + K4)   

#########################################################################
##### Anregung ##########################################################
#########################################################################

class class_Anregung:
    def __init__(self, PARAMETERVEKTOR_ANREGUNG, KNOTEN):        
        self.Knoteni = KNOTEN[0]
        self.Input_Art_der_Anregung = PARAMETERVEKTOR_ANREGUNG[0]
        self.Input_Richtung_Anregung = PARAMETERVEKTOR_ANREGUNG[1]
        self.Input_Amplitude_Anregung = PARAMETERVEKTOR_ANREGUNG[2]
        self.Input_Dauer_der_Anregung = PARAMETERVEKTOR_ANREGUNG[3]
    def Berechnung_Anregung(self, t):
        amplitude = 0
        omega = 0
        if (t < self.Input_Dauer_der_Anregung):
            if (self.Input_Art_der_Anregung == "viertel_sinuswelle"):
                omega = 0.5*np.pi/self.Input_Dauer_der_Anregung
            if (self.Input_Art_der_Anregung == "halbe_sinuswelle"):
                omega = np.pi/self.Input_Dauer_der_Anregung            
            if (self.Input_Richtung_Anregung == "nach_oben"):
                amplitude = self.Input_Amplitude_Anregung
            if (self.Input_Richtung_Anregung == "nach_unten"):
                amplitude = -self.Input_Amplitude_Anregung        
            if (self.Input_Richtung_Anregung == "nach_rechts"):
                amplitude = self.Input_Amplitude_Anregung
            if (self.Input_Richtung_Anregung == "nach_links"):
                amplitude = -self.Input_Amplitude_Anregung        
            if ((self.Input_Richtung_Anregung == "nach_oben") or (self.Input_Richtung_Anregung == "nach_unten")):
                self.Knoteni.z    = amplitude * np.sin(omega*t)
                self.Knoteni.z_p  = omega * amplitude * np.cos(omega*t)
                self.Knoteni.z_pp = -omega * omega * amplitude * np.sin(omega*t)                    
            if ((self.Input_Richtung_Anregung == "nach_rechts") or (self.Input_Richtung_Anregung == "nach_links")):
                self.Knoteni.y    = amplitude * np.sin(omega*t)
                self.Knoteni.y_p  = omega * amplitude * np.cos(omega*t)
                self.Knoteni.y_pp = -omega * omega * amplitude * np.sin(omega*t)
        else:        
            self.Knoteni.x_p = 0
            self.Knoteni.x_pp = 0
            self.Knoteni.y_p = 0
            self.Knoteni.y_pp = 0
            self.Knoteni.z_p = 0
            self.Knoteni.z_pp = 0



#########################################################################
##### Definition der massenlosen Komponenten ############################
#########################################################################

class class_dALEMBERT_Primary_Rolle_rotx:
    def __init__(self, theta, dALEMBERT_max_fehler_rotx_pp, dALEMBERT_schrittweite_rotx_pp, KNOTEN):        
        self.Knoteni = KNOTEN[0]        
        self.theta = theta                  
        self.dALEMBERT_max_fehler_rotx_pp   = dALEMBERT_max_fehler_rotx_pp
        self.dALEMBERT_schrittweite_rotx_pp = dALEMBERT_schrittweite_rotx_pp
        self.dALEMBERT_konvergenz           = False
        self.E_kin                = 0.5 * self.theta * pow(self.Knoteni.rotx_p,2)
        self.E_pot                = 0
        self.E_diss               = 0
        self.E_kin_p              = 0
        self.E_pot_p              = 0
        self.E_diss_p             = 0
    
    def return_Konvergenz(self):
        return self.dALEMBERT_konvergenz
        
    def set_Konvergenz(self, Boolscher_Wert):        
        self.dALEMBERT_konvergenz = Boolscher_Wert
        
    def Berechnung_Kraefte(self):
        temp_value = 1
        
    def Einzelschritt_Anpassung_rotz_pp(self):
        minimierungsfunktion = self.Knoteni.rotx_pp - self.Beschleunigung_aus_Moment()        
        if (abs(minimierungsfunktion) < self.dALEMBERT_max_fehler_rotx_pp):
            self.dALEMBERT_konvergenz = True
        if (minimierungsfunktion > 0):
            self.Knoteni.rotx_pp -= self.dALEMBERT_schrittweite_rotx_pp
        if (minimierungsfunktion < 0):
            self.Knoteni.rotx_pp += self.dALEMBERT_schrittweite_rotx_pp
        
    def Beschleunigung_aus_Moment(self):        
        return self.Knoteni.M_x / self.theta
        
    def Festsetzung_Ergebnisse(self, delta_t):        
        rotx_pp = self.Knoteni.rotx_pp
        rotx_p  = self.Knoteni.rotx_p + rotx_pp * delta_t
        rotx    = self.Knoteni.rotx + rotx_p * delta_t
        self.Knoteni.rotx_pp  = rotx_pp
        self.Knoteni.rotx_p   = rotx_p
        self.Knoteni.rotx     = rotx
                
    def Berechnung_Energien(self, delta_t):                         
        E_kin             = 0.5 * self.theta * pow(self.Knoteni.rotx_p,2)
        E_pot             = 0
        E_diss            = 0
        self.E_kin_p      = (E_kin - self.E_kin) / delta_t
        self.E_pot_p      = (E_pot - self.E_pot) / delta_t
        self.E_diss_p     = (E_diss - self.E_diss) / delta_t
        self.E_kin        = E_kin
        self.E_pot        = E_pot
        self.E_diss      += E_diss
        
class class_dALEMBERT_Secondary_Masse_z:
    def __init__(self, m, g, KNOTEN):        
        self.Knoteni = KNOTEN[0]        
        self.m = m
        self.g = g
        self.Knoteni.F_z      += self.m * (self.g - self.Knoteni.z_pp)
        self.E_kin             = 0.5 * self.m * pow(self.Knoteni.z_p,2)
        self.E_pot             = self.m * -self.g * self.Knoteni.z
        self.E_diss            = 0
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0
        
    def Berechnung_Kraefte(self):        
        self.Knoteni.F_z += self.m * (self.g  - self.Knoteni.z_pp)
                
    def Berechnung_Energien(self, delta_t):                         
        E_kin             = 0.5 * self.m * pow(self.Knoteni.z_p,2)
        E_pot             = self.m * -self.g * self.Knoteni.z
        E_diss            = 0        
        self.E_kin_p      = (E_kin - self.E_kin) / delta_t
        self.E_pot_p      = (E_pot - self.E_pot) / delta_t
        self.E_diss_p     = (E_diss - self.E_diss) / delta_t
        self.E_kin        = E_kin
        self.E_pot        = E_pot
        self.E_diss      += E_diss
        
        
class class_Verbindung_dALEMBERT_Secondary_z_an_Primary_rotx:
    def __init__(self, R, Richtung, KNOTEN):        
        self.Knoteni = KNOTEN[0]
        self.Knotenj = KNOTEN[1]  
        self.R = R
        self.Richtung = Richtung        
        
    def Beidseitige_Uebertragung_Werte(self):                                        
            self.Knotenj.M_x  += self.Richtung * self.Knoteni.F_z * self.R                
            self.Knoteni.z     = self.Richtung * self.Knotenj.rotx * self.R
            self.Knoteni.z_p   = self.Richtung * self.Knotenj.rotx_p * self.R
            self.Knoteni.z_pp  = self.Richtung * self.Knotenj.rotx_pp * self.R        

class class_Knotenpunkt:
    def __init__(self, KNOTEN):
        self.Knoteni = KNOTEN[0]                
        
class class_Feder_x:
    def __init__(self, c, Fv, KNOTEN):
        self.Knoteni = KNOTEN[0]
        self.Knotenj = KNOTEN[1]        
        self.c = c 
        self.Fv = Fv
        self.xv = Fv / self.c 
        self.weg_aktueller_zeitschritt = self.Knoteni.x - self.Knotenj.x        
        self.E_kin             = 0
        self.E_pot             = 0.5 * self.c * pow(self.weg_aktueller_zeitschritt + self.xv,2)
        self.E_diss            = 0               
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0               
     
    def Berechnung_Vorspannkraft_zur_Beibehaltung_Ruhelage(self, m, g):
        self.Fv = m * -g
        self.xv = self.Fv / self.c
        
    def Berechnung_Kraefte(self):                
        self.weg_aktueller_zeitschritt = self.Knoteni.x - self.Knotenj.x              
        self.Knoteni.F_x   -= (self.c * (self.weg_aktueller_zeitschritt + self.xv))
        self.Knotenj.F_x   += (self.c * (self.weg_aktueller_zeitschritt + self.xv))
        
    def Berechnung_Energien(self, delta_t):         
        self.weg_aktueller_zeitschritt = self.Knoteni.x - self.Knotenj.x        
        E_kin             = 0
        E_pot             = 0.5 * self.c * pow(self.weg_aktueller_zeitschritt + self.xv,2)
        E_diss            = 0
        self.E_kin_p      = (E_kin - self.E_kin) / delta_t
        self.E_pot_p      = (E_pot - self.E_pot) / delta_t
        self.E_diss_p     = (E_diss - self.E_diss) / delta_t
        self.E_kin        = E_kin
        self.E_pot        = E_pot
        self.E_diss      += E_diss
        
class class_Feder_y:
    def __init__(self, c, Fv, KNOTEN):
        self.Knoteni = KNOTEN[0]
        self.Knotenj = KNOTEN[1]        
        self.c = c 
        self.Fv = Fv
        self.yv = Fv / self.c 
        self.weg_aktueller_zeitschritt = self.Knoteni.y - self.Knotenj.y        
        self.E_kin             = 0
        self.E_pot             = 0.5 * self.c * pow(self.weg_aktueller_zeitschritt + self.yv,2)
        self.E_diss            = 0               
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0               
    
    def Berechnung_Vorspannkraft_zur_Beibehaltung_Ruhelage(self, m, g):
        self.Fv = m * -g
        self.yv = self.Fv / self.c
        
    def Berechnung_Kraefte(self):                
        self.weg_aktueller_zeitschritt = self.Knoteni.y - self.Knotenj.y              
        self.Knoteni.F_y   -= (self.c * (self.weg_aktueller_zeitschritt + self.yv))
        self.Knotenj.F_y   += (self.c * (self.weg_aktueller_zeitschritt + self.yv))
        
    def Berechnung_Energien(self, delta_t):         
        self.weg_aktueller_zeitschritt = self.Knoteni.y - self.Knotenj.y
        E_kin             = 0
        E_pot             = 0.5 * self.c * pow(self.weg_aktueller_zeitschritt + self.yv,2)
        E_diss            = 0
        self.E_kin_p      = (E_kin - self.E_kin) / delta_t
        self.E_pot_p      = (E_pot - self.E_pot) / delta_t
        self.E_diss_p     = (E_diss - self.E_diss) / delta_t
        self.E_kin        = E_kin
        self.E_pot        = E_pot
        self.E_diss      += E_diss
           
class class_Feder_z:
    def __init__(self, c, Fv, KNOTEN):
        self.Knoteni = KNOTEN[0]
        self.Knotenj = KNOTEN[1]        
        self.c = c 
        self.Fv = Fv
        self.zv = Fv / self.c 
        self.weg_aktueller_zeitschritt = self.Knoteni.z - self.Knotenj.z        
        self.E_kin             = 0
        self.E_pot             = 0.5 * self.c * pow(self.weg_aktueller_zeitschritt + self.zv,2)
        self.E_diss            = 0               
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0               
    
    def Berechnung_Vorspannkraft_zur_Beibehaltung_Ruhelage(self, m, g):
        self.Fv = m * -g
        self.zv = self.Fv / self.c
        
    def Berechnung_Kraefte(self):                
        self.weg_aktueller_zeitschritt = self.Knoteni.z - self.Knotenj.z
        self.Knoteni.F_z   -= (self.c * (self.weg_aktueller_zeitschritt + self.zv))
        self.Knotenj.F_z   += (self.c * (self.weg_aktueller_zeitschritt + self.zv))
        
    def Berechnung_Energien(self, delta_t):         
        self.weg_aktueller_zeitschritt = self.Knoteni.z - self.Knotenj.z
        E_kin             = 0
        E_pot             = 0.5 * self.c * pow(self.weg_aktueller_zeitschritt + self.zv,2)
        E_diss            = 0
        self.E_kin_p      = (E_kin - self.E_kin) / delta_t
        self.E_pot_p      = (E_pot - self.E_pot) / delta_t
        self.E_diss_p     = (E_diss - self.E_diss) / delta_t
        self.E_kin        = E_kin
        self.E_pot        = E_pot
        self.E_diss      += E_diss
        
class class_Daempfer_x:
    def __init__(self, d, KNOTEN):
        self.Knoteni = KNOTEN[0]
        self.Knotenj = KNOTEN[1]    
        self.d = d 
        self.weg_aktueller_zeitschritt = 0
        self.weg_letzter_zeitschritt = 0
        self.E_kin             = 0
        self.E_pot             = 0
        self.E_diss            = 0  
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0               

    def Berechnung_Kraefte(self):                        
        daempfergeschwindigkeit_x      = self.Knoteni.x_p - self.Knotenj.x_p                
        self.Knoteni.F_x              -= self.d * daempfergeschwindigkeit_x
        self.Knotenj.F_x              += self.d * daempfergeschwindigkeit_x
    
    def Berechnung_Energien(self, delta_t): 
        daempfergeschwindigkeit_x      = self.Knoteni.x_p - self.Knotenj.x_p                
        self.weg_aktueller_zeitschritt = self.Knoteni.x - self.Knotenj.x        
        E_kin             = 0
        E_pot             = 0
        E_diss            = abs(self.d * daempfergeschwindigkeit_x) * abs(self.weg_aktueller_zeitschritt - self.weg_letzter_zeitschritt)        
        self.E_kin_p      = (E_kin - self.E_kin) / delta_t
        self.E_pot_p      = (E_pot - self.E_pot) / delta_t
        self.E_diss_p     = (E_diss - self.E_diss) / delta_t
        self.E_kin        = E_kin
        self.E_pot        = E_pot
        self.E_diss      += E_diss
        self.weg_letzter_zeitschritt   = self.weg_aktueller_zeitschritt
        
class class_Daempfer_y:
    def __init__(self, d, KNOTEN):
        self.Knoteni = KNOTEN[0]
        self.Knotenj = KNOTEN[1]    
        self.d = d 
        self.weg_aktueller_zeitschritt = 0
        self.weg_letzter_zeitschritt = 0
        self.E_kin             = 0
        self.E_pot             = 0
        self.E_diss            = 0  
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0               

    def Berechnung_Kraefte(self):                        
        daempfergeschwindigkeit_y      = self.Knoteni.y_p - self.Knotenj.y_p                
        self.Knoteni.F_y              -= self.d * daempfergeschwindigkeit_y
        self.Knotenj.F_y              += self.d * daempfergeschwindigkeit_y
    
    def Berechnung_Energien(self, delta_t): 
        daempfergeschwindigkeit_y      = self.Knoteni.y_p - self.Knotenj.y_p                
        self.weg_aktueller_zeitschritt = self.Knoteni.y - self.Knotenj.y
        E_kin             = 0
        E_pot             = 0
        E_diss            = abs(self.d * daempfergeschwindigkeit_y) * abs(self.weg_aktueller_zeitschritt - self.weg_letzter_zeitschritt)        
        self.E_kin_p      = (E_kin - self.E_kin) / delta_t
        self.E_pot_p      = (E_pot - self.E_pot) / delta_t
        self.E_diss_p     = (E_diss - self.E_diss) / delta_t
        self.E_kin        = E_kin
        self.E_pot        = E_pot
        self.E_diss      += E_diss
        self.weg_letzter_zeitschritt   = self.weg_aktueller_zeitschritt
        
class class_Daempfer_z:
    def __init__(self, d, KNOTEN):
        self.Knoteni = KNOTEN[0]
        self.Knotenj = KNOTEN[1]    
        self.d = d 
        self.weg_aktueller_zeitschritt = 0
        self.weg_letzter_zeitschritt = 0
        self.E_kin             = 0
        self.E_pot             = 0
        self.E_diss            = 0  
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0               

    def Berechnung_Kraefte(self):                        
        daempfergeschwindigkeit_z      = self.Knoteni.z_p - self.Knotenj.z_p                
        self.Knoteni.F_z              -= self.d * daempfergeschwindigkeit_z
        self.Knotenj.F_z              += self.d * daempfergeschwindigkeit_z
    
    def Berechnung_Energien(self, delta_t): 
        daempfergeschwindigkeit_z      = self.Knoteni.z_p - self.Knotenj.z_p                
        self.weg_aktueller_zeitschritt = self.Knoteni.z - self.Knotenj.z
        E_kin             = 0
        E_pot             = 0
        E_diss            = abs(self.d * daempfergeschwindigkeit_z) * abs(self.weg_aktueller_zeitschritt - self.weg_letzter_zeitschritt)        
        self.E_kin_p      = (E_kin - self.E_kin) / delta_t
        self.E_pot_p      = (E_pot - self.E_pot) / delta_t
        self.E_diss_p     = (E_diss - self.E_diss) / delta_t
        self.E_kin        = E_kin
        self.E_pot        = E_pot
        self.E_diss      += E_diss
        self.weg_letzter_zeitschritt   = self.weg_aktueller_zeitschritt
        
class class_Coulomb_x:
    def __init__(self, Fr, KNOTEN):
        self.Knoteni = KNOTEN[0]
        self.Knotenj = KNOTEN[1]        
        self.Fr = Fr
        self.weg_aktueller_zeitschritt = 0
        self.weg_letzter_zeitschritt = 0
        self.E_kin             = 0
        self.E_pot             = 0
        self.E_diss            = 0  
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0               

    def Berechnung_Kraefte(self):                        
        if (self.Knoteni.x_p > 0):
            self.Knoteni.F_x -= self.Fr            
        else:
            self.Knoteni.F_x += self.Fr 
            
    def Berechnung_Energien(self, delta_t): 
        self.weg_aktueller_zeitschritt = self.Knoteni.x        
        E_kin             = 0
        E_pot             = 0
        E_diss            = abs(self.Fr) * abs(self.weg_aktueller_zeitschritt - self.weg_letzter_zeitschritt)        
        self.E_kin_p      = (E_kin - self.E_kin) / delta_t
        self.E_pot_p      = (E_pot - self.E_pot) / delta_t
        self.E_diss_p     = (E_diss - self.E_diss) / delta_t
        self.E_kin        = E_kin
        self.E_pot        = E_pot
        self.E_diss      += E_diss
        self.weg_letzter_zeitschritt   = self.weg_aktueller_zeitschritt
        
class class_Coulomb_y:
    def __init__(self, Fr, KNOTEN):
        self.Knoteni = KNOTEN[0]
        self.Knotenj = KNOTEN[1]        
        self.Fr = Fr
        self.weg_aktueller_zeitschritt = 0
        self.weg_letzter_zeitschritt = 0
        self.E_kin             = 0
        self.E_pot             = 0
        self.E_diss            = 0  
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0               

    def Berechnung_Kraefte(self):                        
        if (self.Knoteni.y_p > 0):
            self.Knoteni.F_y -= self.Fr            
        else:
            self.Knoteni.F_y += self.Fr 
            
    def Berechnung_Energien(self, delta_t): 
        self.weg_aktueller_zeitschritt = self.Knoteni.y
        E_kin             = 0
        E_pot             = 0
        E_diss            = abs(self.Fr) * abs(self.weg_aktueller_zeitschritt - self.weg_letzter_zeitschritt)        
        self.E_kin_p      = (E_kin - self.E_kin) / delta_t
        self.E_pot_p      = (E_pot - self.E_pot) / delta_t
        self.E_diss_p     = (E_diss - self.E_diss) / delta_t
        self.E_kin        = E_kin
        self.E_pot        = E_pot
        self.E_diss      += E_diss
        self.weg_letzter_zeitschritt   = self.weg_aktueller_zeitschritt
        
class class_Coulomb_z:
    def __init__(self, Fr, KNOTEN):
        self.Knoteni = KNOTEN[0]
        self.Knotenj = KNOTEN[1]      
        self.Fr = Fr
        self.weg_aktueller_zeitschritt = 0
        self.weg_letzter_zeitschritt = 0
        self.E_kin             = 0
        self.E_pot             = 0
        self.E_diss            = 0  
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0               

    def Berechnung_Kraefte(self):                        
        if (self.Knoteni.z_p > 0):
            self.Knoteni.F_z -= self.Fr            
        else:
            self.Knoteni.F_z += self.Fr 
            
    def Berechnung_Energien(self, delta_t): 
        self.weg_aktueller_zeitschritt = self.Knoteni.z
        E_kin             = 0
        E_pot             = 0
        E_diss            = abs(self.Fr) * abs(self.weg_aktueller_zeitschritt - self.weg_letzter_zeitschritt)        
        self.E_kin_p      = (E_kin - self.E_kin) / delta_t
        self.E_pot_p      = (E_pot - self.E_pot) / delta_t
        self.E_diss_p     = (E_diss - self.E_diss) / delta_t
        self.E_kin        = E_kin
        self.E_pot        = E_pot
        self.E_diss      += E_diss
        self.weg_letzter_zeitschritt   = self.weg_aktueller_zeitschritt
        
class class_Elastischer_Stoss_x:
    def __init__(self, m1, m2, k, KNOTEN):
        self.m1_Knoteni = KNOTEN[0]              
        self.m2_Knoteni = KNOTEN[1]
        self.m1 = m1
        self.m2 = m2
        self.k = k
        self.E_kin             = 0
        self.E_pot             = 0
        self.E_diss            = 0 
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0               

    def Berechnung_Kraefte(self):                                
        #Kollision
        v1_vorher = self.m1_Knoteni.x_p
        v2_vorher = self.m2_Knoteni.x_p        
        if (self.m1_Knoteni.x < self.m2_Knoteni.x) and (self.m1_Knoteni.x_p < 0):
            v1_nachher = (self.m1 * v1_vorher + self.m2 * v2_vorher - self.m2 * (v1_vorher - v2_vorher) * self.k)/(self.m1 + self.m2)
            v2_nachher = (self.m1 * v1_vorher + self.m2 * v2_vorher - self.m1 * (v2_vorher - v1_vorher) * self.k)/(self.m1 + self.m2)
            self.m1_Knoteni.x_p = v1_nachher            
            self.m2_Knoteni.x_p = v2_nachher   
            Energie_vorher = 0.5 * self.m1 * pow(v1_vorher,2) + 0.5 * self.m2 * pow(v2_vorher,2)
            Energie_nachher = 0.5 * self.m1 * pow(v1_nachher,2) + 0.5 * self.m2 * pow(v2_nachher,2)            
            self.E_diss += self.m1 * self.m2 / (2 * (self.m1 + self.m2)) * pow(v1_vorher - v2_vorher,2) * (1-pow(self.k,2))

    def Berechnung_Energien(self, delta_t):      
        self.E_kin             = 0
        self.E_pot             = 0
        self.E_diss            = 0 
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0   

class class_Elastischer_Stoss_y:
    def __init__(self, m1, m2, k, KNOTEN):
        self.m1_Knoteni = KNOTEN[0]                  
        self.m2_Knoteni = KNOTEN[1]          
        self.m1 = m1
        self.m2 = m2
        self.k = k
        self.E_kin             = 0
        self.E_pot             = 0
        self.E_diss            = 0 
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0               

    def Berechnung_Kraefte(self):                                
        #Kollision        
        if (self.m1_Knoteni.y < self.m2_Knoteni.y) and (self.m1_Knoteni.y_p < 0):
            v1_vorher = self.m1_Knoteni.y_p
            v2_vorher = self.m2_Knoteni.y_p       
            v1_nachher = (self.m1 * v1_vorher + self.m2 * v2_vorher - self.m2 * (v1_vorher - v2_vorher) * self.k)/(self.m1 + self.m2)
            v2_nachher = (self.m1 * v1_vorher + self.m2 * v2_vorher - self.m1 * (v2_vorher - v1_vorher) * self.k)/(self.m1 + self.m2)
            self.m1_Knoteni.y_p = v1_nachher            
            self.m2_Knoteni.y_p = v2_nachher
            Energie_vorher = 0.5 * self.m1 * pow(v1_vorher,2) + 0.5 * self.m2 * pow(v2_vorher,2)
            Energie_nachher = 0.5 * self.m1 * pow(v1_nachher,2) + 0.5 * self.m2 * pow(v2_nachher,2)            
            self.E_diss += self.m1 * self.m2 / (2 * (self.m1 + self.m2)) * pow(v1_vorher - v2_vorher,2) * (1-pow(self.k,2))

    def Berechnung_Energien(self, delta_t):             
        self.E_kin             = 0
        self.E_pot             = 0
        self.E_diss            = 0 
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0   

class class_Elastischer_Stoss_z:
    def __init__(self, m1, m2, k, KNOTEN):
        self.m1_Knoteni = KNOTEN[0]                     
        self.m2_Knoteni = KNOTEN[1]               
        self.m1 = m1
        self.m2 = m2
        self.k = k
        self.E_kin             = 0
        self.E_pot             = 0
        self.E_diss            = 0 
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0               

    def Berechnung_Kraefte(self):                                
        #Kollision
        v1_vorher = self.m1_Knoteni.z_p
        v2_vorher = self.m2_Knoteni.z_p        
        if (self.m1_Knoteni.z < self.m2_Knoteni.z) and (self.m1_Knoteni.z_p < 0):
            v1_nachher = (self.m1 * v1_vorher + self.m2 * v2_vorher - self.m2 * (v1_vorher - v2_vorher) * self.k)/(self.m1 + self.m2)
            v2_nachher = (self.m1 * v1_vorher + self.m2 * v2_vorher - self.m1 * (v2_vorher - v1_vorher) * self.k)/(self.m1 + self.m2)
            self.m1_Knoteni.z_p = v1_nachher            
            self.m2_Knoteni.z_p = v2_nachher  
            Energie_vorher = 0.5 * self.m1 * pow(v1_vorher,2) + 0.5 * self.m2 * pow(v2_vorher,2)
            Energie_nachher = 0.5 * self.m1 * pow(v1_nachher,2) + 0.5 * self.m2 * pow(v2_nachher,2)            
            self.E_diss += self.m1 * self.m2 / (2 * (self.m1 + self.m2)) * pow(v1_vorher - v2_vorher,2) * (1-pow(self.k,2))

    def Berechnung_Energien(self, delta_t):             
        self.E_kin             = 0
        self.E_pot             = 0
        self.E_diss            = 0 
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0           
        
#########################################################################
##### Definition der massenbehafteten Komponenten #######################
#########################################################################
        
class class_DGL_Masse_trans_x:
    """
    Gleichungssystem
    du0/dt = u1
    du1/dt = g + Fi/m - Fj/m
    """
    """    
    """
    def __init__(self, m, g, KNOTEN):
        self.m, self.g, = m, g        
        self.Knoteni = KNOTEN[0]                 
        self.RungeKutta4 = RungeKutta4(self, KNOTEN, "class_DGL_Masse_trans_x")
        self.E_kin             = 0.5 * self.m * pow(self.Knoteni.x_p,2)
        self.E_pot             = self.m * -self.g * self.Knoteni.x
        self.E_diss            = 0      
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0               

    def Berechnung_Kraefte(self):         
        self.Knoteni.F_x += self.m * self.g
        
    def Berechnung_Energien(self, delta_t):         
        E_kin             = 0.5 * self.m * pow(self.Knoteni.x_p,2)      
        E_pot             = self.m * -self.g * self.Knoteni.x      
        E_diss            = 0
        self.E_kin_p      = (E_kin - self.E_kin) / delta_t
        self.E_pot_p      = (E_pot - self.E_pot) / delta_t
        self.E_diss_p     = (E_diss - self.E_diss) / delta_t
        self.E_kin        = E_kin
        self.E_pot        = E_pot
        self.E_diss      += E_diss        
    
    def Aktualisiere_Anfangsrandbedingungen(self):
        self.RungeKutta4.aktualisiere_initial_conditions("class_DGL_Masse_trans_x")    
        
    def Loesung_Differentialgleichung_Zeitschritt(self, delta_t):
        self.Aktualisiere_Anfangsrandbedingungen()
        Rueckgabevektor            =  self.RungeKutta4.solve_timestep(delta_t)
        self.Zuweisung_motion(Rueckgabevektor, delta_t)  
        self.Aktualisiere_Anfangsrandbedingungen()

    def __call__(self, u, KRAFT, t):
        u0, u1 = u
        m, g = self.m, self.g        
        return np.array([u1, (KRAFT[0]/m)])    

    def Zuweisung_motion(self, u, delta_t):
        x_pp = (u[1,1] - self.Knoteni.x_p) / delta_t        
        x_p  = u[1,1]        
        x    = u[1,0]        
        self.Knoteni.x_pp   = x_pp
        self.Knoteni.x_p    = x_p
        self.Knoteni.x      = self.Knoteni.x + x        
                      

class class_DGL_Masse_trans_y:
    """
    Gleichungssystem
    du0/dt = u1
    du1/dt = g + Fi/m - Fj/m
    """
    """    
    """
    def __init__(self, m, g, KNOTEN):
        self.m, self.g, = m, g        
        self.Knoteni = KNOTEN[0]                 
        self.RungeKutta4 = RungeKutta4(self, KNOTEN, "class_DGL_Masse_trans_y")
        self.E_kin             = 0.5 * self.m * pow(self.Knoteni.y_p,2)
        self.E_pot             = self.m * -self.g * self.Knoteni.y
        self.E_diss            = 0      
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0               

    def Berechnung_Kraefte(self):         
        self.Knoteni.F_y += self.m * self.g
        
    def Berechnung_Energien(self, delta_t): 
        E_kin             = 0.5 * self.m * pow(self.Knoteni.y_p,2)      
        E_pot             = self.m * -self.g * self.Knoteni.y      
        E_diss            = 0
        self.E_kin_p      = (E_kin - self.E_kin) / delta_t
        self.E_pot_p      = (E_pot - self.E_pot) / delta_t
        self.E_diss_p     = (E_diss - self.E_diss) / delta_t
        self.E_kin        = E_kin
        self.E_pot        = E_pot
        self.E_diss      += E_diss        
    
    def Aktualisiere_Anfangsrandbedingungen(self):
        self.RungeKutta4.aktualisiere_initial_conditions("class_DGL_Masse_trans_y")    
        
    def Loesung_Differentialgleichung_Zeitschritt(self, delta_t):
        self.Aktualisiere_Anfangsrandbedingungen()
        Rueckgabevektor            =  self.RungeKutta4.solve_timestep(delta_t)
        self.Zuweisung_motion(Rueckgabevektor, delta_t)  
        self.Aktualisiere_Anfangsrandbedingungen()

    def __call__(self, u, KRAFT, t):
        u0, u1 = u
        m, g = self.m, self.g        
        return np.array([u1, (KRAFT[0]/m)])    

    def Zuweisung_motion(self, u, delta_t):
        y_pp = (u[1,1] - self.Knoteni.y_p) / delta_t        
        y_p  = u[1,1]        
        y    = u[1,0]        
        self.Knoteni.y_pp   = y_pp
        self.Knoteni.y_p    = y_p
        self.Knoteni.y      = self.Knoteni.y + y     
        

class class_DGL_Masse_trans_z:
    """
    Gleichungssystem
    du0/dt = u1
    du1/dt = g + Fi/m - Fj/m
    """
    """    
    """
    def __init__(self, m, g, KNOTEN):
        self.m, self.g, = m, g        
        self.Knoteni = KNOTEN[0]                 
        self.RungeKutta4 = RungeKutta4(self, KNOTEN, "class_DGL_Masse_trans_z")
        self.E_kin             = 0.5 * self.m * pow(self.Knoteni.z_p,2)
        self.E_pot             = self.m * -self.g * self.Knoteni.z_p
        self.E_diss            = 0      
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0               

    def Berechnung_Kraefte(self):         
        self.Knoteni.F_z += self.m * self.g
        
    def Berechnung_Energien(self, delta_t): 
        E_kin             = 0.5 * self.m * pow(self.Knoteni.z_p,2)      
        E_pot             = self.m * -self.g * self.Knoteni.z      
        E_diss            = 0
        self.E_kin_p      = (E_kin - self.E_kin) / delta_t
        self.E_pot_p      = (E_pot - self.E_pot) / delta_t
        self.E_diss_p     = (E_diss - self.E_diss) / delta_t
        self.E_kin        = E_kin
        self.E_pot        = E_pot
        self.E_diss      += E_diss        
    
    def Aktualisiere_Anfangsrandbedingungen(self):
        self.RungeKutta4.aktualisiere_initial_conditions("class_DGL_Masse_trans_z")    
        
    def Loesung_Differentialgleichung_Zeitschritt(self, delta_t):
        self.Aktualisiere_Anfangsrandbedingungen()
        Rueckgabevektor            =  self.RungeKutta4.solve_timestep(delta_t)
        self.Zuweisung_motion(Rueckgabevektor, delta_t)  
        self.Aktualisiere_Anfangsrandbedingungen()

    def __call__(self, u, KRAFT, t):
        u0, u1 = u
        m, g = self.m, self.g        
        return np.array([u1, (KRAFT[0]/m)])    

    def Zuweisung_motion(self, u, delta_t):
        z_pp = (u[1,1] - self.Knoteni.z_p) / delta_t        
        z_p  = u[1,1]        
        z    = u[1,0]        
        self.Knoteni.z_pp   = z_pp
        self.Knoteni.z_p    = z_p
        self.Knoteni.z      = self.Knoteni.z + z        

class class_DGL_gelagerteRolle_rot_x:
    """
    Gleichungssystem
    du0/dt = u1
    du1/dt = g + Fi/m - Fj/m
    """
    """    
    """
    def __init__(self, Theta, KNOTEN):
        self.Theta = Theta
        self.Knoteni       = KNOTEN[0]            
        self.RungeKutta4 = RungeKutta4(self, KNOTEN, "class_DGL_gelagerteRolle_rot_x")
        self.E_kin             = 0.5 * self.Theta * pow(self.Knoteni.rotx_p,2)
        self.E_pot             = 0
        self.E_diss            = 0
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0

    def Berechnung_Kraefte(self):         
        temp_value = 1
        
    def Berechnung_Energien(self, delta_t): 
        E_kin             = 0.5 * self.Theta * pow(self.Knoteni.rotx_p,2)        
        E_pot             = 0
        E_diss            = 0
        self.E_kin_p      = (E_kin - self.E_kin) / delta_t
        self.E_pot_p      = (E_pot - self.E_pot) / delta_t
        self.E_diss_p     = (E_diss - self.E_diss) / delta_t
        self.E_kin        = E_kin
        self.E_pot        = E_pot
        self.E_diss      += E_diss
    
    def Aktualisiere_Anfangsrandbedingungen(self):
        self.RungeKutta4.aktualisiere_initial_conditions("class_DGL_gelagerteRolle_rot_x")    
        
    def Loesung_Differentialgleichung_Zeitschritt(self, delta_t):
        self.Aktualisiere_Anfangsrandbedingungen()
        Rueckgabevektor            =  self.RungeKutta4.solve_timestep(delta_t)
        self.Zuweisung_motion(Rueckgabevektor, delta_t)  
        self.Aktualisiere_Anfangsrandbedingungen()

    def __call__(self, u, KRAFT, t):
        u0, u1 = u
        Theta = self.Theta
        return np.array([u1, (KRAFT[0]/Theta)])    

    def Zuweisung_motion(self, u, delta_t):
        rotx_pp = (u[1,1] - self.Knoteni.rotx_p) / delta_t        
        rotx_p  = u[1,1]        
        rotx    = u[1,0]        
        self.Knoteni.rotx_pp   = rotx_pp                
        self.Knoteni.rotx_p    = rotx_p
        self.Knoteni.rotx      = self.Knoteni.rotx + rotx
        
        
class class_DGL_Federpendel_rot_x:    
    def __init__(self, l, m, g, KNOTEN):        
        self.Knoteni       = KNOTEN[0]            
        self.Knotenj       = KNOTEN[1]            
        self.l             = l
        self.m             = m
        self.g             = g
        self.Theta         = m * pow(l,2)
        self.RungeKutta4   = RungeKutta4(self, KNOTEN, "class_DGL_Federpendel_rot_x")
        self.E_kin             = 0.5 * self.Theta * pow(self.Knoteni.rotx_p,2)
        self.E_pot             = self.m * -self.g * self.Knotenj.z
        self.E_diss            = 0
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0    
        
    def Berechnung_Kraefte(self):         
        self.Knoteni.M_x = self.m * (self.g - self.Knoteni.z_pp) * self.l * math.sin(self.Knoteni.rotx)
        
    def Berechnung_Energien(self, delta_t): 
        E_kin             = 0.5 * self.Theta * pow(self.Knoteni.rotx_p,2)        
        E_pot             = self.m * -self.g * self.Knotenj.z
        E_diss            = 0
        self.E_kin_p      = (E_kin - self.E_kin) / delta_t
        self.E_pot_p      = (E_pot - self.E_pot) / delta_t
        self.E_diss_p     = (E_diss - self.E_diss) / delta_t
        self.E_kin        = E_kin
        self.E_pot        = E_pot
        self.E_diss      += E_diss
    
    def Aktualisiere_Anfangsrandbedingungen(self):
        self.RungeKutta4.aktualisiere_initial_conditions("class_DGL_Federpendel_rot_x")    
        
    def Loesung_Differentialgleichung_Zeitschritt(self, delta_t):
        self.Aktualisiere_Anfangsrandbedingungen()
        Rueckgabevektor            =  self.RungeKutta4.solve_timestep(delta_t)
        self.Zuweisung_motion(Rueckgabevektor, delta_t)  
        self.Aktualisiere_Anfangsrandbedingungen()

    def __call__(self, u, KRAFT, t):
        u0, u1 = u
        Theta = self.Theta
        return np.array([u1, (KRAFT[0]/Theta)])    

    def Zuweisung_motion(self, u, delta_t):
        rotx_pp = (u[1,1] - self.Knoteni.rotx_p) / delta_t        
        rotx_p  = u[1,1]        
        rotx    = u[1,0]        
        self.Knoteni.rotx_pp   = rotx_pp                
        self.Knoteni.rotx_p    = rotx_p
        self.Knoteni.rotx      = self.Knoteni.rotx + rotx
        self.Knotenj.y         = self.l * math.sin(self.Knoteni.rotx)
        self.Knotenj.z         = self.Knoteni.z - self.l * math.cos(self.Knoteni.rotx)
        
class class_DGL_Federpendel_mit_Seilkraft_rot_x:    
    def __init__(self, l, m, g, c, KNOTEN):        
        self.Knoteni       = KNOTEN[0]            
        self.Knotenj       = KNOTEN[1]            
        self.l             = l
        self.m             = m
        self.g             = g
        self.c             = c
        self.Theta         = m * pow(l,2)
        self.RungeKutta4   = RungeKutta4(self, KNOTEN, "class_DGL_Federpendel_mit_Seilkraft_rot_x")
        self.E_kin             = 0.5 * self.Theta * pow(self.Knoteni.rotx_p,2)
        self.E_pot             = self.m * -self.g * self.Knotenj.z
        self.E_diss            = 0
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0
        self.Translation_y     = 0
        self.Translation_z     = 0
        
    def Berechnung_Kraefte(self):  
        #Berechnung aktueller Winkel ohne Zuweisung (kann durch Anregung abweichen)
        rotx = math.atan((self.Knotenj.y-self.Knoteni.y)/(self.Knoteni.z-self.Knotenj.z))        
        #Berechnung Längenänderung Faden
        laengenanderung_faden = pow(pow(self.Knoteni.y - self.Knotenj.y,2) + pow(self.Knoteni.z - self.Knotenj.z,2),0.5) - self.l        
        #Berechnung Kraft im Faden (positiv = Zugkraft)
        kraft_faden           = laengenanderung_faden * self.c
        self.Knoteni.F_y      = kraft_faden * -math.sin(rotx)
        self.Knoteni.F_z      = kraft_faden * math.cos(rotx)        
        #Korrektur_Länge_Faden (Verkürzung um delta_y und delta_z, so dass Kraft minimiert wird)
        delta_y = self.Knoteni.F_y / self.c
        delta_z = self.Knoteni.F_z / self.c        
        if (self.Knotenj.y < self.Knotenj.y):
            self.Translation_y += delta_y
        else:
            self.Knotenj.y -= delta_y
        self.Translation_z += delta_z
        #Rückstellmoment aus Gravitation mit realem Winkel berechnen
        self.Knoteni.M_x = -self.m * (self.g - self.Knoteni.z_pp) * self.l * math.sin(rotx)

        
    def Berechnung_Energien(self, delta_t): 
        E_kin             = 0.5 * self.Theta * pow(self.Knoteni.rotx_p,2)        
        E_pot             = self.m * -self.g * self.Knotenj.z
        E_diss            = 0
        self.E_kin_p      = (E_kin - self.E_kin) / delta_t
        self.E_pot_p      = (E_pot - self.E_pot) / delta_t
        self.E_diss_p     = (E_diss - self.E_diss) / delta_t
        self.E_kin        = E_kin
        self.E_pot        = E_pot
        self.E_diss      += E_diss
    
    def Aktualisiere_Anfangsrandbedingungen(self):
        self.RungeKutta4.aktualisiere_initial_conditions("class_DGL_Federpendel_mit_Seilkraft_rot_x")    
        
    def Loesung_Differentialgleichung_Zeitschritt(self, delta_t):
        self.Aktualisiere_Anfangsrandbedingungen()
        Rueckgabevektor            =  self.RungeKutta4.solve_timestep(delta_t)
        self.Zuweisung_motion(Rueckgabevektor, delta_t)  
        self.Aktualisiere_Anfangsrandbedingungen()

    def __call__(self, u, KRAFT, t):
        u0, u1 = u
        Theta = self.Theta
        return np.array([u1, (KRAFT[0]/Theta)])    

    def Zuweisung_motion(self, u, delta_t):
        rotx_pp = (u[1,1] - self.Knoteni.rotx_p) / delta_t        
        rotx_p  = u[1,1]        
        rotx    = u[1,0]        
        self.Knoteni.rotx_pp   = rotx_pp                
        self.Knoteni.rotx_p    = rotx_p
        self.Knoteni.rotx      = self.Knoteni.rotx + rotx
        self.Knotenj.y         = self.Translation_y - self.l * math.sin(self.Knoteni.rotx)
        self.Knotenj.z         = self.Translation_z + self.Knoteni.z - self.l * math.cos(self.Knoteni.rotx)
        
#########################################################################
##### Erstellung aller Knoten ###########################################
#########################################################################

K1            = class_Knoten()
K2            = class_Knoten()
K3            = class_Knoten()
K4            = class_Knoten()
K5            = class_Knoten()
K6            = class_Knoten()
K7            = class_Knoten()
K8            = class_Knoten()
K9            = class_Knoten()
K10           = class_Knoten()
K11           = class_Knoten()
K12           = class_Knoten()
K13           = class_Knoten()
K14           = class_Knoten()
K15           = class_Knoten()
K16           = class_Knoten()
K17           = class_Knoten()
K18           = class_Knoten()
K19           = class_Knoten()
K20           = class_Knoten()

#########################################################################
##### Reset_Knotenkraft #################################################
#########################################################################
        
    
def Reset_Knotenkraefte():
    K1.Reset_Knotenkraft()
    K2.Reset_Knotenkraft()
    K3.Reset_Knotenkraft()
    K4.Reset_Knotenkraft()
    K5.Reset_Knotenkraft()
    K6.Reset_Knotenkraft()
    K7.Reset_Knotenkraft()
    K8.Reset_Knotenkraft()
    K9.Reset_Knotenkraft()
    K10.Reset_Knotenkraft()
    K11.Reset_Knotenkraft()
    K12.Reset_Knotenkraft()
    K13.Reset_Knotenkraft()
    K14.Reset_Knotenkraft()
    K15.Reset_Knotenkraft()
    K16.Reset_Knotenkraft()
    K17.Reset_Knotenkraft()
    K18.Reset_Knotenkraft()
    K19.Reset_Knotenkraft()
    K20.Reset_Knotenkraft()


def Reset_Knoten():
    K1.Reset_Knoten()
    K2.Reset_Knoten()
    K3.Reset_Knoten()
    K4.Reset_Knoten()
    K5.Reset_Knoten()
    K6.Reset_Knoten()
    K7.Reset_Knoten()
    K8.Reset_Knoten()
    K9.Reset_Knoten()
    K10.Reset_Knoten()
    K11.Reset_Knoten()
    K12.Reset_Knoten()
    K13.Reset_Knoten()
    K14.Reset_Knoten()
    K15.Reset_Knoten()
    K16.Reset_Knoten()
    K17.Reset_Knoten()
    K18.Reset_Knoten()
    K19.Reset_Knoten()
    K20.Reset_Knoten()

#########################################################################
##### Datenlogger #######################################################
#########################################################################

datenlogger_t                  = []

datenlogger_K1        = class_Datenlogger_Knoten()
datenlogger_K2        = class_Datenlogger_Knoten()
datenlogger_K3        = class_Datenlogger_Knoten()
datenlogger_K4        = class_Datenlogger_Knoten()
datenlogger_K5        = class_Datenlogger_Knoten()
datenlogger_K6        = class_Datenlogger_Knoten()
datenlogger_K7        = class_Datenlogger_Knoten()
datenlogger_K8        = class_Datenlogger_Knoten()
datenlogger_K9        = class_Datenlogger_Knoten()
datenlogger_K10       = class_Datenlogger_Knoten()
datenlogger_K11       = class_Datenlogger_Knoten()
datenlogger_K12       = class_Datenlogger_Knoten()
datenlogger_K13       = class_Datenlogger_Knoten()
datenlogger_K14       = class_Datenlogger_Knoten()
datenlogger_K15       = class_Datenlogger_Knoten()
datenlogger_K16       = class_Datenlogger_Knoten()
datenlogger_K17       = class_Datenlogger_Knoten()
datenlogger_K18       = class_Datenlogger_Knoten()
datenlogger_K19       = class_Datenlogger_Knoten()
datenlogger_K20       = class_Datenlogger_Knoten()

datenlogger_E1        = class_Datenlogger_Element()
datenlogger_E2        = class_Datenlogger_Element()
datenlogger_E3        = class_Datenlogger_Element()
datenlogger_E4        = class_Datenlogger_Element()
datenlogger_E5        = class_Datenlogger_Element()
datenlogger_E6        = class_Datenlogger_Element()
datenlogger_E7        = class_Datenlogger_Element()
datenlogger_E8        = class_Datenlogger_Element()
datenlogger_E9        = class_Datenlogger_Element()
datenlogger_E10       = class_Datenlogger_Element()
datenlogger_E11       = class_Datenlogger_Element()
datenlogger_E12       = class_Datenlogger_Element()
datenlogger_E13       = class_Datenlogger_Element()
datenlogger_E14       = class_Datenlogger_Element()
datenlogger_E15       = class_Datenlogger_Element()
datenlogger_E16       = class_Datenlogger_Element()
datenlogger_E17       = class_Datenlogger_Element()
datenlogger_E18       = class_Datenlogger_Element()
datenlogger_E19       = class_Datenlogger_Element()
datenlogger_E20       = class_Datenlogger_Element()
    
def Datenlogger_leeren():
    datenlogger_t.clear()

    datenlogger_K1.Werte_loeschen()
    datenlogger_K2.Werte_loeschen()
    datenlogger_K3.Werte_loeschen()
    datenlogger_K4.Werte_loeschen()
    datenlogger_K5.Werte_loeschen()
    datenlogger_K6.Werte_loeschen()
    datenlogger_K7.Werte_loeschen()
    datenlogger_K8.Werte_loeschen()
    datenlogger_K9.Werte_loeschen()
    datenlogger_K10.Werte_loeschen()
    datenlogger_K11.Werte_loeschen()
    datenlogger_K12.Werte_loeschen()
    datenlogger_K13.Werte_loeschen()
    datenlogger_K14.Werte_loeschen()
    datenlogger_K15.Werte_loeschen()
    datenlogger_K16.Werte_loeschen()
    datenlogger_K17.Werte_loeschen()
    datenlogger_K18.Werte_loeschen()
    datenlogger_K19.Werte_loeschen()
    datenlogger_K20.Werte_loeschen()

    datenlogger_E1.Werte_loeschen()
    datenlogger_E2.Werte_loeschen()
    datenlogger_E3.Werte_loeschen()
    datenlogger_E4.Werte_loeschen()
    datenlogger_E5.Werte_loeschen()
    datenlogger_E6.Werte_loeschen()
    datenlogger_E7.Werte_loeschen()
    datenlogger_E8.Werte_loeschen()
    datenlogger_E9.Werte_loeschen()
    datenlogger_E10.Werte_loeschen()
    datenlogger_E11.Werte_loeschen()
    datenlogger_E12.Werte_loeschen()
    datenlogger_E13.Werte_loeschen()
    datenlogger_E14.Werte_loeschen()
    datenlogger_E15.Werte_loeschen()
    datenlogger_E16.Werte_loeschen()
    datenlogger_E17.Werte_loeschen()
    datenlogger_E18.Werte_loeschen()
    datenlogger_E19.Werte_loeschen()
    datenlogger_E20.Werte_loeschen()

def Datenlogger_Knoten_schreiben(t):
    datenlogger_t.append(t) 
    datenlogger_K1.Werte_anhaengen(K1)           
    datenlogger_K2.Werte_anhaengen(K2)           
    datenlogger_K3.Werte_anhaengen(K3)           
    datenlogger_K4.Werte_anhaengen(K4)           
    datenlogger_K5.Werte_anhaengen(K5)           
    datenlogger_K6.Werte_anhaengen(K6)           
    datenlogger_K7.Werte_anhaengen(K7)           
    datenlogger_K8.Werte_anhaengen(K8)           
    datenlogger_K9.Werte_anhaengen(K9)           
    datenlogger_K10.Werte_anhaengen(K10)           
    datenlogger_K11.Werte_anhaengen(K11)           
    datenlogger_K12.Werte_anhaengen(K12)           
    datenlogger_K13.Werte_anhaengen(K13)           
    datenlogger_K14.Werte_anhaengen(K14)           
    datenlogger_K15.Werte_anhaengen(K15)           
    datenlogger_K16.Werte_anhaengen(K16)           
    datenlogger_K17.Werte_anhaengen(K17)           
    datenlogger_K18.Werte_anhaengen(K18)           
    datenlogger_K19.Werte_anhaengen(K19)           
    datenlogger_K20.Werte_anhaengen(K20)     
    
    
#########################################################################
##### csv-Datei schreiben ###############################################
#########################################################################

def comma_float(num):
            ''' wandelt float in string mit echtem Komma um '''
            return sub(r'\.', ',', str(num))
    
def csv_schreiben(Dateiname, Knoten):
    with open(Dateiname, 'w', newline='') as output_file:
        writer = csv.writer(output_file, delimiter=";")
        writer.writerow(["t","x","x_p","x_pp","y", "y_p","y_pp","z","z_p", "z_pp","F_x","F_y","F_z"])
        for point in range(len(datenlogger_t)):            
            writer.writerow([comma_float(datenlogger_t[point]), comma_float(Knoten.x[point]), comma_float(Knoten.x_p[point]), comma_float(Knoten.x_pp[point]), comma_float(Knoten.y[point]), comma_float(Knoten.y_p[point]), comma_float(Knoten.y_pp[point]), comma_float(Knoten.z[point]), comma_float(Knoten.z_p[point]), comma_float(Knoten.z_pp[point]), comma_float(Knoten.F_x[point]), comma_float(Knoten.F_y[point]), comma_float(Knoten.F_z[point])])
            
            
#########################################################################
##### Kennfelder ########################################################
#########################################################################

class class_Messung_Demonstrator:
    def __init__(self,Dateiname):
        self.pos = []        
        self.t = []
        self.anregung = []    
        self.antwort = []  
        self.Dateiname = Dateiname
        with open(Dateiname,'r') as csvdatei:
            csv_reader_object = csv.reader(csvdatei, delimiter=',')
            zeilennummer = 0
            for row in csv_reader_object:
                if (zeilennummer == 0):            
                    zeilennummer += 1            
                else:
                    self.pos.append(zeilennummer)
                    self.t.append(float(row[0]))
                    self.anregung.append(float(row[1])/1000)        
                    self.antwort.append(float(row[2])/1000)        
                    zeilennummer += 1
        