# -*- coding: utf-8 -*-
"""
The code is made to calculate wavelengths (l) at which the diffracion of light occures in 1D or 2D gratings. 

"""
import sympy as sp
import numpy as np
dx, dy, m, l, t, f, n = sp.symbols('dx dy m l t f n')
dim=int(input("1D or 2D structure (input 1 or 2) ")) 
dx_val =float(input("Period along x direction (nm) = ")) 
with open('output.txt', 'w') as file:
    file.write("Angle(deg) -Lambda(nm) +Lambda(nm)\n")
    if dim == 1:
        eq = (((2*sp.pi/l)*sp.sin(t*sp.pi/180)*sp.cos(f*sp.pi/180))+(2*sp.pi*m/dx))**2-(2*sp.pi/l)**2
        f_val = float(input("Azimuthal angle (deg) = ")) #(0 if not applicable)
        m_val = int(input("Diffraction order = ")) #Integer for the diffraction order(absolute value e.g. 1, 2, 3...)
        t_start = int(input("Start angle of light incidence (deg) = "))
        t_end = int(input("End angle of light incidence (deg) = "))
        t_step = float(input("Calculation step (deg) = "))
        t_val=sp.symbols('t_val')
        t_val_list=np.arange(t_start, t_end, t_step, dtype=float)
        solutions_neg = []
        solutions_pos = []
        for t_val in t_val_list:
            solution = sp.solve(eq.subs({t: t_val, f: f_val, m: m_val, dx: dx_val}), l)
            solutions_neg.append(abs(solution[0]))
            solutions_pos.append(abs(solution[1]))
        for i, t_val in enumerate(t_val_list):
            print(f"t = {t_val}(deg), lambda(-) = {solutions_neg[i].evalf()}(nm), lambda(+) = {solutions_pos[i].evalf()}(nm)")
            file.write(f"{t_val} {solutions_neg[i]} {solutions_pos[i]}\n")
    elif dim == 2:
        eq = (((2*sp.pi/l)*sp.sin(t*sp.pi/180)*sp.cos(f*sp.pi/180))+(2*sp.pi*m/dx))**2+(((2*sp.pi/l)*sp.sin(t*sp.pi/180)*sp.sin(f*sp.pi/180))+(2*sp.pi*n/dy))**2-(2*sp.pi/l)**2
        dy_val = float(input("Period along y direction (nm) = "))
        m_val = int(input("Diffraction order for x direction = ")) #Integer for the diffraction order(absolute value e.g. 1, 2, 3...), associated with x direction (0 if abscent)
        f_val = float(input("Azimuthal angle (deg) = ")) #(0 if not applicable)
        n_val = int(input("Diffraction order for y direction = ")) #Integer for the diffraction order(absolute value e.g. 1, 2, 3...), associated with y direction (0 if abscent)
        t_start = int(input("Start angle of light incidence (deg) = "))
        t_end = int(input("End angle of light incidence (deg) = "))
        t_step = float(input("Calculation step (deg) = "))
        t_val=sp.symbols('t_val')
        t_val_list=np.arange(t_start, t_end, t_step, dtype=float)
        solutions_neg = []
        solutions_pos = []
        for t_val in t_val_list:
            solution = sp.solve(eq.subs({t: t_val, f: f_val, m: m_val, dx: dx_val, dy: dy_val, n: n_val}), l)
            solutions_neg.append(abs(solution[0]))
            solutions_pos.append(abs(solution[1]))
        for i, t_val in enumerate(t_val_list):
            print(f"t = {t_val}(deg), lambda(-) = {solutions_neg[i].evalf()}(nm), lambda(+) = {solutions_pos[i].evalf()}(nm)")
            file.write(f"{t_val} {solutions_neg[i]} {solutions_pos[i]}\n")  


        

