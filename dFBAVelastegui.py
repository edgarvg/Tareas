# -*- coding: utf-8 -*-
"""
Created on Sat Dec 15 12:41:20 2018

@author: EDGAR
"""
import numpy as np
from scipy.integrate import odeint
# Libraries for the Metabolic block
import cobra
import os
from os.path import join
data_dir="."
model=cobra.io.read_sbml_model(join(data_dir, "iMT1026v3.xml"))
model.objective="Ex_biomass"
# Plotting libraries
import matplotlib.pyplot as plt

def main():
        # Kinetic block
    def kineticBlock(Met): #ocupamos el etanol pq hay inhibicion por producto
        # Parameters
        vmetmax, K_met,= 2.5, 0.88
        f_Sorb = 0.1 #dato del paper
        v_met=vmetmax*Met/(K_met+Met)
        LB_Sorb =v_met*f_Sorb 
            
        return [v_met,LB_Sorb]
    
    # Metabolic block
    def metabolicBlock(v_met,LB_Sorb):
        #Methanol
        model.reactions.get_by_id("Ex_meoh").upper_bound=-v_met
        model.reactions.get_by_id("Ex_meoh").lower_bound=-v_met
    #Glycerol
        model.reactions.get_by_id("Ex_glyc").upper_bound=0
        model.reactions.get_by_id("Ex_glyc").lower_bound=0
    #Glucosa
        model.reactions.get_by_id("Ex_glc_D").upper_bound=0
        model.reactions.get_by_id("Ex_glc_D").lower_bound=0
    #sorbitol exchange
        
        model.reactions.get_by_id("Ex_sbt_D").upper_bound=-LB_Sorb
        model.reactions.get_by_id("Ex_sbt_D").lower_bound=-LB_Sorb
    #Oxygen
        #model.reactions.get_by_id("Ex_o2").upper_bound=-LB_O2
        #model.reactions.get_by_id("Ex_o2").lower_bound=-LB_O2
        
    
        solution = model.optimize()
        u = 5*solution.f
        v_sorb = model.reactions.get_by_id("Ex_sbt_D").x
        v_fab = model.reactions.get_by_id("Ex_fab").x
        v_AOD = model.reactions.get_by_id("AOD").x
        v_DAS = model.reactions.get_by_id("DAS").x
        v_CATp = model.reactions.get_by_id("CATp").x
        return [u,v_sorb,v_fab,v_AOD,v_DAS,v_CATp]
        #DEFINIR LOWER BOUNDS
    # Dynamic block
    def f(y,t,params):
        V,VX,VMet,VSorb = y # Current values
        F,u, v_met, v_sorb = params  # unpack parameters
        Met_F = 50 
        MW_Met,MW_Sorb= [0.32,0.182] #Molecular weights
        derivs=[F,                  # dV/dt
                u*VX,                # dVX/dt
                F*Met_F-v_met*MW_Met*(VX), # dVMet/dt
                v_sorb*MW_Sorb*(VX)] # dVSorb/dt
        return derivs
    def dynamicBlock(y,params, ti,tf):
        time=np.linspace(ti,tf,100)
        soln = odeint(f,y,time,args=(params,))
        # Get solutions at the final time point (tf):
        V=soln[-1,0]
        X=soln[-1,1]/V
        Met=soln[-1,2]/V
        Sorb=soln[-1,3]/V #ojo hay q dividir para volumen para eliminar Ej XV, solo necesito X
        return [V,X,Met,Sorb]    
    
    
    # Miscelaneous functions
    # Feed flow 
    def F(t):
        F=0.3 #para F constante
        return F
    
    # Save results along the fermentation path
    u_path,V_path,X_path = [],[],[]
    Met_path,Sorb_path = [],[]
    v_met_path,LB_Sorb_path,v_fab_path,v_AOD_path,v_DAS_path,v_CATp_path= [],[],[],[],[],[]
    def savePath(u,V,X,Met,v_met,LB_Sorb,v_fab,v_AOD,v_DAS,v_CATp):
        global u_path,V_path,X_path,LB_Sorb_path
        global Met_path,Sorb_path
        global v_met_path,LB_Sorb_path,v_fab_path,v_AOD_path,v_DAS_path,v_CATp_path
        u_path += [u]
        V_path += [V]
        X_path += [X]
        Met_path += [Met]
        Sorb_path += [Sorb]
        v_met_path += [v_met]
        LB_Sorb_path += [LB_Sorb]
        v_fab_path += [v_fab]
        v_AOD_path += [v_AOD]
        v_DAS_path += [v_DAS]
        v_CATp_path += [v_CATp]
        # Initial conditions
    Met=20
    Sorb=20
    V,X=[0.5,5] #volumen, inóculo inicial
    
    # Running the simulation over time
    time=np.linspace(0,24,200) #dividido en 200 intervalos
    for i in range(len(time)):
        # KINETIC BLOCK: 
        v_met,LB_Sorb = kineticBlock(Met)
        # METABOLIC BLOCK
        u,v_sorb,v_fab,v_AOD,v_DAS,v_CATp= metabolicBlock(v_met,LB_Sorb)
        # DYNAMIC BLOCK
        
        if i==len(time)-1: continue #si el loop llega al penultimo punto, ya no haga la ecuacion diferencial y pare
        y = [V,X*V,Met*V,Sorb*V]
        params = [F(time[i]),u,v_met, v_sorb]
        V,X,Met,Sorb = dynamicBlock(y, params, time[i],time[i+1])  
        # Save results along the fermentaion path
        savePath(u,V,X,Met,Sorb,LB_Sorb,v_fab,v_AOD,v_DAS,v_CATp)   #Almacenamiento en un vector
        plt.plot(time[1:200],Met_path,'r',linewidth=2,label='Met') 
    plt.plot(time[1:200],X_path,'b',linewidth=2,label='X')
    plt.plot(time[1:200],Sorb_path,'g',linewidth=2,label='Sorb')
    #plt.plot(time[1:200],u_path,'b--',linewidth=2,label='u')
    plt.legend()
    plt.ylabel('Biomass, Methanol and Sorbitol [g/L]')
    plt.xlabel('Time [h]')
    plt.show()
    plt.plot(time[1:200],v_AOD_path,'r',linewidth=2,label='AOX')
    plt.plot(time[1:200],v_DAS_path,'r',linewidth=2,label='DAS')                                       
    plt.plot(time[1:200],v_CATp_path,'r',linewidth=2,label='CATp')
    
    #plt.plot(time[1:200],LB_Sorb_path,'r',linewidth=2,label='vsorb') 
    plt.legend()
    plt.ylabel('AOX, DAS and CATp [mmol/gDCW h]')
    plt.xlabel('Time [h]')
    plt.show()
if __name__=='__main__':
    main()