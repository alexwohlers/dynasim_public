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

from lsdt_dynasim_solver import*

#########################################################################
##### Einstellungen für Plots ###########################################
#########################################################################

plt.rcParams['figure.dpi']  = 300
plt.rcParams['savefig.dpi'] = 300

#########################################################################
##### Definition der Knoten #############################################
#########################################################################

class class_Knoten_Masse:
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
       
class class_Datenlogger_Knoten_Masse:
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

class class_Datenlogger_Element_Masse:
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
    def Potentielle_Energie_Referenzpunktverschiebung(self):
        offset = min(self.E_pot)
        for point in range(len(self.E_pot)):
            self.E_pot[point] -= offset

########################################################################################################################

def ball_radius(mass):
    #Berechnet aus gegebener Masse und Dichte den Radius einer Vollkugel
    density = 7850
    volume = mass/density
    radius = ((3*volume)/(4*np.pi))**(1/3)
    return radius
    
def eigenkreisfrequenz_berechnen(liste_zeit,liste_schwingung):
    #Erst einmal Mittelwert bestimmen
    temp_summe = 0    
    for index_liste in range(len(liste_schwingung)):
        temp_summe += liste_schwingung[index_liste]    
    temp_mittelwert = temp_summe / len(liste_schwingung)    
    #Nun Anzahl der Durchgänge durch Mittelwert bestimmen und aus Zeit die Periodendauer bestimmen    
    temp_anzahl_positver_nulldurchgaenge = 0
    temp_index_erster_durchgang = 0
    temp_index_letzter_durchgang = 0
    for index_liste in range(len(liste_schwingung)-1):
        if (liste_schwingung[index_liste] < temp_mittelwert) and (liste_schwingung[index_liste+1] > temp_mittelwert):
            temp_anzahl_positver_nulldurchgaenge += 1
            if (temp_index_erster_durchgang == 0):
                temp_index_erster_durchgang = index_liste
            else:
                temp_index_letzter_durchgang = index_liste
    temp_Periodendauer = (liste_zeit[temp_index_letzter_durchgang]-liste_zeit[temp_index_erster_durchgang])/(temp_anzahl_positver_nulldurchgaenge-1)    
    print("Periodendauer: "+str(round(temp_Periodendauer,3)))
    print("Frequenz: "+str(round(1/temp_Periodendauer,3)))
    print("Eigenkreisfrequenz: "+str(round(2*3.1415/temp_Periodendauer,3)))
    
   

#########################################################################
##### Element zur Anbindungen eines FluidSim-Knotens  ###################
#########################################################################
   
   
class class_FluidSim_p_nach_DynaSim_F_Druck_auf_Flaeche:
    def __init__(self, Flaeche, Kraftrichtung, KNOTEN):   
        self.Knoteni_Masse = KNOTEN[0]
        self.Knotenj_Fluid = KNOTEN[1]
        self.Flaeche       = Flaeche    
        self.Kraftrichtung = Kraftrichtung        
    def Berechnung_Kraft_aus_Druck(self):
        if (self.Kraftrichtung == 'Kraft_x_negative_Richtung'):
            self.Knoteni_Masse.F_x -= self.Knotenj_Fluid.p * self.Flaeche            
        if (self.Kraftrichtung == 'Kraft_x_positive_Richtung'):
            self.Knoteni_Masse.F_x += self.Knotenj_Fluid.p * self.Flaeche            


#########################################################################
##### Anregung ##########################################################
#########################################################################

class class_Kraftanregung_x_Masse:
    def __init__(self, Kraftamplitude, Anregungsfrequenz, KNOTEN):   
        self.Knoteni = KNOTEN[0]
        self.Kraftamplitude     = Kraftamplitude
        self.Anregungsfrequenz = Anregungsfrequenz
    def Berechnung_Anregung(self, t):
        self.Knoteni.F_x += self.Kraftamplitude * math.sin(self.Anregungsfrequenz * t)
    
class class_Anregung_Masse:
    def __init__(self, Modus, PARAMETERVEKTOR_ANREGUNG, KNOTEN):        
        self.Modus = Modus
        self.Knoteni = KNOTEN[0]
        self.Knotenj = KNOTEN[1]
        self.Input_Art_der_Anregung = PARAMETERVEKTOR_ANREGUNG[0]
        self.Input_Richtung_Anregung = PARAMETERVEKTOR_ANREGUNG[1]
        self.Input_Amplitude_Anregung = PARAMETERVEKTOR_ANREGUNG[2]
        self.Input_Dauer_der_Anregung = PARAMETERVEKTOR_ANREGUNG[3]
        self.E_kin                = 0
        self.E_pot                = 0
        self.E_diss               = 0
        self.E_kin_p              = 0
        self.E_pot_p              = 0
        self.E_diss_p             = 0        
        
    def Berechnung_Anregung(self, t):
        amplitude = 0
        omega = 0
        y_vorher   = self.Knoteni.y        
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
            if (self.Modus == 'Fadenpendel'):
                self.Knoteni.rotx   += math.atan((self.Knoteni.y-y_vorher)/abs(self.Knoteni.z-self.Knotenj.z))            
        else:        
            self.Knoteni.x_p = 0
            self.Knoteni.x_pp = 0
            self.Knoteni.y_p = 0
            self.Knoteni.y_pp = 0
            self.Knoteni.z_p = 0
            self.Knoteni.z_pp = 0
    
    def Berechnung_Energien(self, delta_t):                         
        E_kin             = 0
        E_pot             = 0        
        E_diss            = 0                
        self.E_kin_p      = (E_kin - self.E_kin) / delta_t
        self.E_pot_p      = (E_pot - self.E_pot) / delta_t
        self.E_diss_p     = (E_diss - self.E_diss) / delta_t
        self.E_kin        = E_kin
        self.E_pot        = E_pot
        self.E_diss      += E_diss        


#########################################################################
##### Definition der massenlosen Komponenten ############################
#########################################################################

class class_dALEMBERT_Primary_Rolle_rotx_Masse:
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
        
class class_dALEMBERT_Secondary_Masse_z_Masse:
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
        
        
class class_Verbindung_dALEMBERT_Secondary_z_an_Primary_rotx_Masse:
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

class class_Knotenpunkt_Masse:
    def __init__(self, KNOTEN):
        self.Knoteni = KNOTEN[0]                
              
class class_Feder_x_Masse:
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
     
    def Berechnung_Vorspannkraft_zur_Beibehaltung_Nulllage(self, m, g):
        self.Fv = m * -g
        self.xv = self.Fv / self.c
        
    def Berechnung_Vorspannkraft_zur_Beibehaltung_Referenzlage(self, Referenzfedersteifigkeit, m, g):                  
        self.xv = m * -g / self.c - m * -g / Referenzfedersteifigkeit
        
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
        
class class_Feder_y_Masse:
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
    
    def Berechnung_Vorspannkraft_zur_Beibehaltung_Nulllage(self, m, g):
        self.Fv = m * -g
        self.yv = self.Fv / self.c
        
    def Berechnung_Vorspannkraft_zur_Beibehaltung_Referenzlage(self, Referenzfedersteifigkeit, m, g):                  
        self.yv = m * -g / self.c - m * -g / Referenzfedersteifigkeit
        
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
           
class class_Feder_z_Masse:
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
    
    def Berechnung_Vorspannkraft_zur_Beibehaltung_Nulllage(self, m, g):
        self.Fv = m * -g
        self.zv = self.Fv / self.c
        
    def Berechnung_Vorspannkraft_zur_Beibehaltung_Referenzlage(self, Ruhelage_z, m, g):                  
        #self.zv = m * -g / self.c - m * -g / Referenzfedersteifigkeit
        self.zv = m * -g / self.c + Ruhelage_z
        
        
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

class class_Blattfeder_beidseitig_fest_eingespannt_x_Masse:
    def __init__(self, L, E, I, Fv, KNOTEN):
        self.Knoteni = KNOTEN[0]
        self.Knotenj = KNOTEN[1]        
        self.L                 = L     
        self.E                 = E
        self.I                 = I             
        self.c                 = 24 * self.E * self.I / pow(self.L,3)
        self.Fv                = Fv
        self.xv                = Fv / self.c 
        self.weg_aktueller_zeitschritt = self.Knoteni.x - self.Knotenj.x        
        self.E_kin             = 0
        self.E_pot             = 0.5 * self.c * pow(self.weg_aktueller_zeitschritt + self.xv,2)
        self.E_diss            = 0               
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0               
     
    def Berechnung_Vorspannkraft_zur_Beibehaltung_Nulllage(self, m, g):
        self.Fv = m * -g
        self.xv = self.Fv / self.c
        
    def Berechnung_Vorspannkraft_zur_Beibehaltung_Referenzlage(self, Referenzfedersteifigkeit, m, g):                  
        self.xv = m * -g / self.c - m * -g / Referenzfedersteifigkeit
        
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
        
class class_Daempfer_x_Masse:
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
        
class class_Daempfer_y_Masse:
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
        
class class_Daempfer_z_Masse:
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
        
class class_Coulomb_x_Masse:
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
        
class class_Coulomb_y_Masse:
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
        
class class_Coulomb_z_Masse:
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

class class_Coulomb_Mx_Masse:
    def __init__(self, Mr, KNOTEN):
        self.Knoteni = KNOTEN[0]        
        self.Mr = Mr
        self.winkel_aktueller_zeitschritt = 0
        self.winkel_letzter_zeitschritt = 0
        self.E_kin             = 0
        self.E_pot             = 0
        self.E_diss            = 0  
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0               

    def Berechnung_Kraefte(self):                        
        if (self.Knoteni.rotx_p > 0):
            self.Knoteni.M_x -= self.Mr            
        else:
            self.Knoteni.M_x += self.Mr 
            
    def Berechnung_Energien(self, delta_t): 
        self.winkel_aktueller_zeitschritt = self.Knoteni.rotx
        E_kin             = 0
        E_pot             = 0
        E_diss            = abs(self.Mr) * abs(self.winkel_aktueller_zeitschritt - self.winkel_letzter_zeitschritt)        
        self.E_kin_p      = (E_kin - self.E_kin) / delta_t
        self.E_pot_p      = (E_pot - self.E_pot) / delta_t
        self.E_diss_p     = (E_diss - self.E_diss) / delta_t
        self.E_kin        = E_kin
        self.E_pot        = E_pot
        self.E_diss      += E_diss
        self.winkel_letzter_zeitschritt   = self.winkel_aktueller_zeitschritt
        
class class_Luftreibung_Masse:
    def __init__(self, cw, rho, A, KNOTEN):
        self.Knoteni = KNOTEN[0]
        self.cw = cw
        self.rho = rho
        self.A = A        
        self.E_kin             = 0
        self.E_pot             = 0
        self.E_diss            = 0  
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0               

    def Berechnung_Kraefte(self):                                
        self.Knoteni.F_x              -= 0.5 * self.cw * self.rho * self.A * pow(self.Knoteni.x_p,2)
        self.Knoteni.F_y              -= 0.5 * self.cw * self.rho * self.A * pow(self.Knoteni.y_p,2)
        self.Knoteni.F_z              -= 0.5 * self.cw * self.rho * self.A * pow(self.Knoteni.z_p,2)
    
    def Berechnung_Energien(self, delta_t):                 
        E_kin             = 0
        E_pot             = 0
        E_diss            = 0
        self.E_kin_p      = (E_kin - self.E_kin) / delta_t
        self.E_pot_p      = (E_pot - self.E_pot) / delta_t
        self.E_diss_p     = (E_diss - self.E_diss) / delta_t
        self.E_kin        = E_kin
        self.E_pot        = E_pot
        self.E_diss      += E_diss        
        
class class_Elastischer_Stoss_x_Masse:
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
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0   

class class_Elastischer_Stoss_y_Masse:
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
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0   

class class_Elastischer_Stoss_z_Masse:
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
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0           
        
#########################################################################
##### Definition der massenbehafteten Komponenten #######################
#########################################################################
        
class class_DGL_Masse_trans_x_Masse:
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
        self.RungeKutta4_Masse = RungeKutta4_Masse(self, KNOTEN, "class_DGL_Masse_trans_x_Masse")
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
        self.RungeKutta4_Masse.aktualisiere_initial_conditions("class_DGL_Masse_trans_x_Masse")    
        
    def Loesung_Differentialgleichung_Zeitschritt(self, delta_t):
        self.Aktualisiere_Anfangsrandbedingungen()
        Rueckgabevektor            =  self.RungeKutta4_Masse.solve_timestep(delta_t)
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
                      

class class_DGL_Masse_trans_y_Masse:
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
        self.RungeKutta4_Masse = RungeKutta4_Masse(self, KNOTEN, "class_DGL_Masse_trans_y_Masse")
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
        self.RungeKutta4_Masse.aktualisiere_initial_conditions("class_DGL_Masse_trans_y_Masse")    
        
    def Loesung_Differentialgleichung_Zeitschritt(self, delta_t):
        self.Aktualisiere_Anfangsrandbedingungen()
        Rueckgabevektor            =  self.RungeKutta4_Masse.solve_timestep(delta_t)
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
        

class class_DGL_Masse_trans_z_Masse:
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
        self.RungeKutta4_Masse = RungeKutta4_Masse(self, KNOTEN, "class_DGL_Masse_trans_z_Masse")
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
        self.RungeKutta4_Masse.aktualisiere_initial_conditions("class_DGL_Masse_trans_z_Masse")    
        
    def Loesung_Differentialgleichung_Zeitschritt(self, delta_t):
        self.Aktualisiere_Anfangsrandbedingungen()
        Rueckgabevektor            =  self.RungeKutta4_Masse.solve_timestep(delta_t)
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

class class_DGL_gelagerteRolle_rot_x_Masse:
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
        self.RungeKutta4_Masse = RungeKutta4_Masse(self, KNOTEN, "class_DGL_gelagerteRolle_rot_x_Masse")
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
        self.RungeKutta4_Masse.aktualisiere_initial_conditions("class_DGL_gelagerteRolle_rot_x_Masse")    
        
    def Loesung_Differentialgleichung_Zeitschritt(self, delta_t):
        self.Aktualisiere_Anfangsrandbedingungen()
        Rueckgabevektor            =  self.RungeKutta4_Masse.solve_timestep(delta_t)
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

class class_DGL_Rolle_rot_y_Masse:
    """
    Rolle, die auf Oberfläche rollt
    Es können zwei Kräfte angreifen in Knoten j und k
    Im Knoten i dürfen von außen keine Elemente angeschlossen sein!
    """
    """    
    """
    def __init__(self, m, R, aj, ak, KNOTEN):        
        self.Knoteni           = KNOTEN[0] #Drehpunkt            
        self.Knotenj           = KNOTEN[1] #Krafteinleitung 1                   
        self.Knotenk           = KNOTEN[2] #Krafteinleitung 2           
        self.m                 = m                
        self.R                 = R
        self.Theta             = m * pow(R,2)
        self.aj                = aj        #Abstand Drehpunkt zu Krafteinleitung 1  (positiver Abstand = oben, negativer Abstand = unten)
        self.ak                = ak        #Abstand Drehpunkt zu Krafteinleitung 1  (positiver Abstand = oben, negativer Abstand = unten)
        self.RungeKutta4_Masse       = RungeKutta4_Masse(self, KNOTEN, "class_DGL_Rolle_rot_y_Masse")
        self.E_kin             = 0.5 * self.Theta * pow(self.Knoteni.roty_p,2)
        self.E_pot             = 0
        self.E_diss            = 0
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0

    def Berechnung_Kraefte(self):         
        temp_value = 1
        
    def Berechnung_Energien(self, delta_t): 
        E_kin             = 0.5 * self.Theta * pow(self.Knoteni.roty_p,2) + 0.5 * self.m * pow(self.Knoteni.x_p,2)         
        E_pot             = 0
        E_diss            = 0
        self.E_kin_p      = (E_kin - self.E_kin) / delta_t
        self.E_pot_p      = (E_pot - self.E_pot) / delta_t
        self.E_diss_p     = (E_diss - self.E_diss) / delta_t
        self.E_kin        = E_kin
        self.E_pot        = E_pot
        self.E_diss      += E_diss
    
    def Aktualisiere_Anfangsrandbedingungen(self):
        #Trägheitskraft im Drehpunkt berechnen
        self.Knoteni.F_x = -self.m * self.Knoteni.x_pp
        self.RungeKutta4_Masse.aktualisiere_initial_conditions("class_DGL_Rolle_rot_y_Masse")    
        
    def Loesung_Differentialgleichung_Zeitschritt(self, delta_t):
        self.Aktualisiere_Anfangsrandbedingungen()
        Rueckgabevektor            =  self.RungeKutta4_Masse.solve_timestep(delta_t)
        self.Zuweisung_motion(Rueckgabevektor, delta_t)  
        self.Aktualisiere_Anfangsrandbedingungen()

    def __call__(self, u, KRAFT, t):
        u0, u1 = u
        Theta = self.Theta
        return np.array([u1, ((KRAFT[0]*(self.R) + KRAFT[1]*(self.R+self.aj) + KRAFT[2]*(self.R+self.ak))/Theta)])    

    def Zuweisung_motion(self, u, delta_t):               
        roty_pp = (u[1,1] - self.Knoteni.roty_p) / delta_t        
        roty_p  = u[1,1]        
        roty    = u[1,0]        
        self.Knoteni.roty_pp   = roty_pp                
        self.Knoteni.roty_p    = roty_p
        self.Knoteni.roty      = self.Knoteni.roty + roty      
        self.Knoteni.x_pp      = self.R * self.Knoteni.roty_pp
        self.Knotenj.x_pp      = (self.R + self.aj) * self.Knoteni.roty_pp
        self.Knotenk.x_pp      = (self.R + self.ak) * self.Knoteni.roty_pp
        self.Knoteni.x_p       = self.R * self.Knoteni.roty_p
        self.Knotenj.x_p       = (self.R + self.aj) * self.Knoteni.roty_p
        self.Knotenk.x_p       = (self.R + self.ak) * self.Knoteni.roty_p
        self.Knoteni.x         = self.R * self.Knoteni.roty        
        self.Knotenj.x         = (self.R + self.aj) * self.Knoteni.roty
        self.Knotenk.x         = (self.R + self.ak) * self.Knoteni.roty
        
        
class class_DGL_Fadenpendel_rot_x_Masse:    
    def __init__(self, l, m, g, KNOTEN):        
        self.Knoteni       = KNOTEN[0]            
        self.Knotenj       = KNOTEN[1]            
        self.l             = l
        self.m             = m
        self.g             = g
        self.Theta         = m * pow(l,2)
        self.RungeKutta4_Masse   = RungeKutta4_Masse(self, KNOTEN, "class_DGL_Fadenpendel_rot_x_Masse")
        self.E_kin             = 0.5 * self.Theta * pow(self.Knoteni.rotx_p,2)
        self.E_pot             = self.m * -self.g * self.Knotenj.z
        self.E_diss            = 0
        self.E_kin_p           = 0
        self.E_pot_p           = 0
        self.E_diss_p          = 0    
        
    def Berechnung_Kraefte(self):         
        #Reibkräfte: self.Knotenj.F_y und self.Knotenj.F_z        
        self.Knoteni.M_x += self.m * (self.g - self.Knoteni.z_pp) * self.l * math.sin(self.Knoteni.rotx)
        
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
        self.RungeKutta4_Masse.aktualisiere_initial_conditions("class_DGL_Fadenpendel_rot_x_Masse")    
        
    def Loesung_Differentialgleichung_Zeitschritt(self, delta_t):
        self.Aktualisiere_Anfangsrandbedingungen()
        Rueckgabevektor            =  self.RungeKutta4_Masse.solve_timestep(delta_t)
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
        self.Knotenj.y_p       = self.l * math.sin(self.Knoteni.rotx_p)
        self.Knotenj.z_p       = self.Knoteni.z - self.l * math.cos(self.Knoteni.rotx_p)
        
class class_DGL_Fadenpendel_mit_Seilkraft_rot_x_Masse:    
    def __init__(self, l, m, g, c, KNOTEN):        
        self.Knoteni       = KNOTEN[0]            
        self.Knotenj       = KNOTEN[1]            
        self.l             = l
        self.m             = m
        self.g             = g
        self.c             = c
        self.Theta         = m * pow(l,2)
        self.RungeKutta4_Masse   = RungeKutta4_Masse(self, KNOTEN, "class_DGL_Fadenpendel_mit_Seilkraft_rot_x_Masse")
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
        self.Knoteni.M_x += -self.m * (self.g - self.Knoteni.z_pp) * self.l * math.sin(rotx)

        
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
        self.RungeKutta4_Masse.aktualisiere_initial_conditions("class_DGL_Fadenpendel_mit_Seilkraft_rot_x_Masse")    
        
    def Loesung_Differentialgleichung_Zeitschritt(self, delta_t):
        self.Aktualisiere_Anfangsrandbedingungen()
        Rueckgabevektor            =  self.RungeKutta4_Masse.solve_timestep(delta_t)
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
        self.Knotenj.y_p       = self.l * math.sin(self.Knoteni.rotx_p)
        self.Knotenj.z_p       = self.Knoteni.z - self.l * math.cos(self.Knoteni.rotx_p)
        
#########################################################################
##### Erstellung aller Knoten ###########################################
#########################################################################

K1_Masse            = class_Knoten_Masse()
K2_Masse            = class_Knoten_Masse()
K3_Masse            = class_Knoten_Masse()
K4_Masse            = class_Knoten_Masse()
K5_Masse            = class_Knoten_Masse()
K6_Masse            = class_Knoten_Masse()
K7_Masse            = class_Knoten_Masse()
K8_Masse            = class_Knoten_Masse()
K9_Masse            = class_Knoten_Masse()
K10_Masse           = class_Knoten_Masse()
K11_Masse           = class_Knoten_Masse()
K12_Masse           = class_Knoten_Masse()
K13_Masse           = class_Knoten_Masse()
K14_Masse           = class_Knoten_Masse()
K15_Masse           = class_Knoten_Masse()
K16_Masse           = class_Knoten_Masse()
K17_Masse           = class_Knoten_Masse()
K18_Masse           = class_Knoten_Masse()
K19_Masse           = class_Knoten_Masse()
K20_Masse           = class_Knoten_Masse()

#########################################################################
##### Reset_Knotenkraft #################################################
#########################################################################
        
    
def Reset_Knotenkraefte_Masse():
    K1_Masse.Reset_Knotenkraft()
    K2_Masse.Reset_Knotenkraft()
    K3_Masse.Reset_Knotenkraft()
    K4_Masse.Reset_Knotenkraft()
    K5_Masse.Reset_Knotenkraft()
    K6_Masse.Reset_Knotenkraft()
    K7_Masse.Reset_Knotenkraft()
    K8_Masse.Reset_Knotenkraft()
    K9_Masse.Reset_Knotenkraft()
    K10_Masse.Reset_Knotenkraft()
    K11_Masse.Reset_Knotenkraft()
    K12_Masse.Reset_Knotenkraft()
    K13_Masse.Reset_Knotenkraft()
    K14_Masse.Reset_Knotenkraft()
    K15_Masse.Reset_Knotenkraft()
    K16_Masse.Reset_Knotenkraft()
    K17_Masse.Reset_Knotenkraft()
    K18_Masse.Reset_Knotenkraft()
    K19_Masse.Reset_Knotenkraft()
    K20_Masse.Reset_Knotenkraft()


def Reset_Knoten_Masse():
    K1_Masse.Reset_Knoten()
    K2_Masse.Reset_Knoten()
    K3_Masse.Reset_Knoten()
    K4_Masse.Reset_Knoten()
    K5_Masse.Reset_Knoten()
    K6_Masse.Reset_Knoten()
    K7_Masse.Reset_Knoten()
    K8_Masse.Reset_Knoten()
    K9_Masse.Reset_Knoten()
    K10_Masse.Reset_Knoten()
    K11_Masse.Reset_Knoten()
    K12_Masse.Reset_Knoten()
    K13_Masse.Reset_Knoten()
    K14_Masse.Reset_Knoten()
    K15_Masse.Reset_Knoten()
    K16_Masse.Reset_Knoten()
    K17_Masse.Reset_Knoten()
    K18_Masse.Reset_Knoten()
    K19_Masse.Reset_Knoten()
    K20_Masse.Reset_Knoten()

#########################################################################
##### Datenlogger #######################################################
#########################################################################

datenlogger_t_Masse                  = []

datenlogger_K1_Masse        = class_Datenlogger_Knoten_Masse()
datenlogger_K2_Masse        = class_Datenlogger_Knoten_Masse()
datenlogger_K3_Masse        = class_Datenlogger_Knoten_Masse()
datenlogger_K4_Masse        = class_Datenlogger_Knoten_Masse()
datenlogger_K5_Masse        = class_Datenlogger_Knoten_Masse()
datenlogger_K6_Masse        = class_Datenlogger_Knoten_Masse()
datenlogger_K7_Masse        = class_Datenlogger_Knoten_Masse()
datenlogger_K8_Masse        = class_Datenlogger_Knoten_Masse()
datenlogger_K9_Masse        = class_Datenlogger_Knoten_Masse()
datenlogger_K10_Masse       = class_Datenlogger_Knoten_Masse()
datenlogger_K11_Masse       = class_Datenlogger_Knoten_Masse()
datenlogger_K12_Masse       = class_Datenlogger_Knoten_Masse()
datenlogger_K13_Masse       = class_Datenlogger_Knoten_Masse()
datenlogger_K14_Masse       = class_Datenlogger_Knoten_Masse()
datenlogger_K15_Masse       = class_Datenlogger_Knoten_Masse()
datenlogger_K16_Masse       = class_Datenlogger_Knoten_Masse()
datenlogger_K17_Masse       = class_Datenlogger_Knoten_Masse()
datenlogger_K18_Masse       = class_Datenlogger_Knoten_Masse()
datenlogger_K19_Masse       = class_Datenlogger_Knoten_Masse()
datenlogger_K20_Masse       = class_Datenlogger_Knoten_Masse()

datenlogger_E1_Masse        = class_Datenlogger_Element_Masse()
datenlogger_E2_Masse        = class_Datenlogger_Element_Masse()
datenlogger_E3_Masse        = class_Datenlogger_Element_Masse()
datenlogger_E4_Masse        = class_Datenlogger_Element_Masse()
datenlogger_E5_Masse        = class_Datenlogger_Element_Masse()
datenlogger_E6_Masse        = class_Datenlogger_Element_Masse()
datenlogger_E7_Masse        = class_Datenlogger_Element_Masse()
datenlogger_E8_Masse        = class_Datenlogger_Element_Masse()
datenlogger_E9_Masse        = class_Datenlogger_Element_Masse()
datenlogger_E10_Masse       = class_Datenlogger_Element_Masse()
datenlogger_E11_Masse       = class_Datenlogger_Element_Masse()
datenlogger_E12_Masse       = class_Datenlogger_Element_Masse()
datenlogger_E13_Masse       = class_Datenlogger_Element_Masse()
datenlogger_E14_Masse       = class_Datenlogger_Element_Masse()
datenlogger_E15_Masse       = class_Datenlogger_Element_Masse()
datenlogger_E16_Masse       = class_Datenlogger_Element_Masse()
datenlogger_E17_Masse       = class_Datenlogger_Element_Masse()
datenlogger_E18_Masse       = class_Datenlogger_Element_Masse()
datenlogger_E19_Masse       = class_Datenlogger_Element_Masse()
datenlogger_E20_Masse       = class_Datenlogger_Element_Masse()
    
def Datenlogger_leeren_Masse():
    datenlogger_t_Masse.clear()

    datenlogger_K1_Masse.Werte_loeschen()
    datenlogger_K2_Masse.Werte_loeschen()
    datenlogger_K3_Masse.Werte_loeschen()
    datenlogger_K4_Masse.Werte_loeschen()
    datenlogger_K5_Masse.Werte_loeschen()
    datenlogger_K6_Masse.Werte_loeschen()
    datenlogger_K7_Masse.Werte_loeschen()
    datenlogger_K8_Masse.Werte_loeschen()
    datenlogger_K9_Masse.Werte_loeschen()
    datenlogger_K10_Masse.Werte_loeschen()
    datenlogger_K11_Masse.Werte_loeschen()
    datenlogger_K12_Masse.Werte_loeschen()
    datenlogger_K13_Masse.Werte_loeschen()
    datenlogger_K14_Masse.Werte_loeschen()
    datenlogger_K15_Masse.Werte_loeschen()
    datenlogger_K16_Masse.Werte_loeschen()
    datenlogger_K17_Masse.Werte_loeschen()
    datenlogger_K18_Masse.Werte_loeschen()
    datenlogger_K19_Masse.Werte_loeschen()
    datenlogger_K20_Masse.Werte_loeschen()

    datenlogger_E1_Masse.Werte_loeschen()
    datenlogger_E2_Masse.Werte_loeschen()
    datenlogger_E3_Masse.Werte_loeschen()
    datenlogger_E4_Masse.Werte_loeschen()
    datenlogger_E5_Masse.Werte_loeschen()
    datenlogger_E6_Masse.Werte_loeschen()
    datenlogger_E7_Masse.Werte_loeschen()
    datenlogger_E8_Masse.Werte_loeschen()
    datenlogger_E9_Masse.Werte_loeschen()
    datenlogger_E10_Masse.Werte_loeschen()
    datenlogger_E11_Masse.Werte_loeschen()
    datenlogger_E12_Masse.Werte_loeschen()
    datenlogger_E13_Masse.Werte_loeschen()
    datenlogger_E14_Masse.Werte_loeschen()
    datenlogger_E15_Masse.Werte_loeschen()
    datenlogger_E16_Masse.Werte_loeschen()
    datenlogger_E17_Masse.Werte_loeschen()
    datenlogger_E18_Masse.Werte_loeschen()
    datenlogger_E19_Masse.Werte_loeschen()
    datenlogger_E20_Masse.Werte_loeschen()

def Datenlogger_Knoten_schreiben_Masse(t):
    datenlogger_t_Masse.append(t) 
    datenlogger_K1_Masse.Werte_anhaengen(K1_Masse)           
    datenlogger_K2_Masse.Werte_anhaengen(K2_Masse)           
    datenlogger_K3_Masse.Werte_anhaengen(K3_Masse)           
    datenlogger_K4_Masse.Werte_anhaengen(K4_Masse)           
    datenlogger_K5_Masse.Werte_anhaengen(K5_Masse)           
    datenlogger_K6_Masse.Werte_anhaengen(K6_Masse)           
    datenlogger_K7_Masse.Werte_anhaengen(K7_Masse)           
    datenlogger_K8_Masse.Werte_anhaengen(K8_Masse)           
    datenlogger_K9_Masse.Werte_anhaengen(K9_Masse)           
    datenlogger_K10_Masse.Werte_anhaengen(K10_Masse)           
    datenlogger_K11_Masse.Werte_anhaengen(K11_Masse)           
    datenlogger_K12_Masse.Werte_anhaengen(K12_Masse)           
    datenlogger_K13_Masse.Werte_anhaengen(K13_Masse)           
    datenlogger_K14_Masse.Werte_anhaengen(K14_Masse)           
    datenlogger_K15_Masse.Werte_anhaengen(K15_Masse)           
    datenlogger_K16_Masse.Werte_anhaengen(K16_Masse)           
    datenlogger_K17_Masse.Werte_anhaengen(K17_Masse)           
    datenlogger_K18_Masse.Werte_anhaengen(K18_Masse)           
    datenlogger_K19_Masse.Werte_anhaengen(K19_Masse)           
    datenlogger_K20_Masse.Werte_anhaengen(K20_Masse)     
    
    
#########################################################################
##### csv-Datei schreiben ###############################################
#########################################################################

def comma_float(num):
            ''' wandelt float in string mit echtem Komma um '''
            return sub(r'\.', ',', str(num))
    
def csv_schreiben_Masse(Dateiname, Knoten_Masse):
    with open(Dateiname, 'w', newline='') as output_file:
        writer = csv.writer(output_file, delimiter=";")
        writer.writerow(["t","x","x_p","x_pp","y", "y_p","y_pp","z","z_p", "z_pp","F_x","F_y","F_z"])
        for point in range(len(datenlogger_t_Masse)):            
            writer.writerow([comma_float(datenlogger_t_Masse[point]), comma_float(Knoten_Masse.x[point]), comma_float(Knoten_Masse.x_p[point]), comma_float(Knoten_Masse.x_pp[point]), comma_float(Knoten_Masse.y[point]), comma_float(Knoten_Masse.y_p[point]), comma_float(Knoten_Masse.y_pp[point]), comma_float(Knoten_Masse.z[point]), comma_float(Knoten_Masse.z_p[point]), comma_float(Knoten_Masse.z_pp[point]), comma_float(Knoten_Masse.F_x[point]), comma_float(Knoten_Masse.F_y[point]), comma_float(Knoten_Masse.F_z[point])])
            
            
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
        