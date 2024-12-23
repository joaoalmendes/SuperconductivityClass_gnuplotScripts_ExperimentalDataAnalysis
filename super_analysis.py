import numpy as np
import matplotlib.pyplot as plt
import re
import sys
import os
#   TP1
def file_to_list(file_path, extra):
    T_H_m = [[], [], []]
    if extra == True: T_H_m = [[], [], [], []]
    with open(file_path, 'r') as f:
        next(f)
        for line in f:
            # Strip any leading/trailing whitespace and split the line by spaces
            str_floats = line.strip().split()
            # Convert each string to a float and store in a list
            for num in range(len(str_floats)): T_H_m[num].append(float(str_floats[num]))
            # Append the list of floats to the main list
    return T_H_m

def calculate_qui(file, path):
    T_H_m = file_to_list(path+file, extra=False).copy()

    rho, m = 6.3, 0.01033    #  g/cm^-3; g
    for value in range(len(T_H_m[2])):  # convert moment to mag
        T_H_m[2][value] *= (rho/m)

    X = []

    for field, mag in zip(T_H_m[1], T_H_m[2]):
        X.append(float(mag/field))

    F = []
    N = 4/3 * np.pi
    for qui in range(len(X)):
        f = (-X[qui]/(1 - N*X[qui]))*4*np.pi
        F.append(f)

    os.chdir(path)

    with open('worked_'+file, 'w') as f:
        f.write('# Temperature (K)  f (adimensional)\n')
        for T, f_value in zip(T_H_m[0], F):
            f.write(f"{T} {f_value}\n")

def Hc_T_data(input_file_path, output_file_path):
    temperatures = []
    averages = []

    pattern = re.compile(r'(\d+K):\s+H_c0\s*=\s*([\d\.\-]+)\s+H_c1\s*=\s*([\d\.\-]+)')

    with open(input_file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue  # Skip empty lines

            # Use regular expression to parse the line
            match = pattern.match(line)
            if not match:
                print(f"Error parsing line: {line}")
                continue
            
            temperature, H_c0, H_c1 = match.groups()

            # Print the extracted H_c0 and H_c1 values
            print("H_c0:", H_c0, "H_c1:", H_c1)

            # Skip the line if any of the values are 'indf'
            if H_c0 == 'indf' or H_c1 == 'indf':
                continue

            # Convert H_c0 and H_c1 to float
            try:
                H_c0 = float(H_c0)
                H_c1 = float(H_c1)
            except ValueError as e:
                print(f"Error converting H_c0 or H_c1 to float: {e}")
                continue

            # Calculate the average
            average = (abs(H_c0) + abs(H_c1)) / 2

            # Append the results to the lists
            temperatures.append(temperature)
            averages.append(average)

    # Write the results to the output file
    with open(output_file_path, 'w') as output_file:
        for temp, avg in zip(temperatures, averages):
            output_file.write(f"{temp} {avg:.2f}\n")

def Jc_H_data(input_file_path, output_file_path):
    results = []

    # Regular expression to match the new pattern
    pattern = re.compile(r'H,M\s*=\s*([\d\.\-]+);\s*([\d\.\-]+)\s*-H,M\s*=\s*([\d\.\-]+);\s*([\d\.\-]+)')
    rho, m = 6.3, 0.01033    #  g/cm^-3; g
    with open(input_file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue  # Skip empty lines

            # Use regular expression to parse the line
            match = pattern.match(line)
            if not match:
                print(f"Error parsing line: {line}")
                continue
            
            # Extract values
            H1, m1, H2, m2 = match.groups()

            # Convert extracted values to float
            try:
                H1 = float(H1)
                M1 = float(m1) * (rho/m)
                H2 = float(H2)
                M2 = float(m2) * (rho/m)
            except ValueError as e:
                print(f"Error converting values to float: {e}")
                continue

            # Perform calculations
            Jc =  1.5*(M2 - M1)/0.2  # emu cm^-4
            H = (abs(H1) + abs(H2)) / 2

            # Append the results to the list
            results.append((H, Jc))
    results.sort(key=lambda x: x[0])
    # Write the results to the output file
    with open(output_file_path, 'w') as output_file:
        for avg_H, sum_M in results:
            output_file.write(f"{avg_H:.2f} {sum_M:.6f}\n")

#   TP2
def X_T(input_file_path, output_file_path):
    T_H_m = file_to_list(output_file_path+input_file_path, extra=True).copy()
    s = input_file_path.split('_')[0]
    if s == 'GD2O3': rho, m, molar_mass = 7.41, 0.01115, 362.5    #  g/cm^-3; g; g/mol
    if s == 'Ho2O3': rho, m, molar_mass = 8.41, 0.02001, 377.86    #  g/cm^-3; g
    if s == 'Nd2O3': rho, m, molar_mass = 7.24, 0.01326, 336.48    #  g/cm^-3; g
    for value in range(len(T_H_m[2])):  # convert moment to mag
        T_H_m[2][value] *= (rho/m)

    X = []
    inv_X = []
    for field, mag in zip(T_H_m[1], T_H_m[2]):
        qui = float((mag/field) * (molar_mass/rho))
        X.append(qui)
        inv_X.append(float(1/qui))
    os.chdir(output_file_path)
    with open('X&invX_T__'+s+'.txt', 'w') as f:
        f.write('# Temperature (K)  X (CGS) 1/X (CGS)\n')
        for T, qui, inv_qui in zip(T_H_m[0], X, inv_X):
            f.write(f"{T} {qui} {inv_qui}\n")

def M_H(input_file_path, output_file_path):
    T_H_m = file_to_list(output_file_path+input_file_path, extra=True).copy()
    chemical, temp = input_file_path.split('_')[0], input_file_path.split('_')[2].split('.')[0]
    rho, m, molar_mass = 8.41, 0.02001, 377.86    #  g/cm^3; g
    for value in range(len(T_H_m[2])):  # convert moment to mag
        T_H_m[2][value] *= 1e-3/((2*6.022*1e23*(m/molar_mass)) * 9.274*1e-24)
    os.chdir(output_file_path)
    with open('mag_H__'+chemical+'_'+temp+'.txt', 'w') as f:
        f.write('# H (Oe)  Mag(emu/cm^3)\n')
        for h, mag in zip(T_H_m[1], T_H_m[2]):
            f.write(f"{h} {mag}\n")

#   TP3
def copper(input_file_path, output_file_path):
    T_H_m = file_to_list(output_file_path+input_file_path, extra=False).copy()
    s = input_file_path.split('_')[0]
    m, rho, molar_mass = 0.01386, 1.882, 199.65*2 #   g; g/cm^3; g/mol
    X_d = -0.5*molar_mass*pow(10,-6)
    X_m = []
    X = []
    for moment in range(len(T_H_m[2])):
        T_H_m[2][moment] += 1.179 * pow(10, -8) * T_H_m[1][moment]
        qui_m = (T_H_m[2][moment]/T_H_m[1][moment]) *  (molar_mass/m)
        X_m.append(float(qui_m))
        X.append(float(qui_m - X_d))
    os.chdir(output_file_path)
    with open('T_Xm_X__'+s+'.txt', 'w') as f:
        f.write('# Temperature (K)  X_m (CGS) X (CGS) 1/X (CGS)\n')
        for T, qui_m, qui in zip(T_H_m[0], X_m, X):
            f.write(f"{T} {qui_m} {qui} {float(1/qui)}\n")
    
path = r'./Dados_TP2/'
current_path = os.getcwd()
file = str(sys.argv[1])

#calculate_qui(file, path)
#Hc_T_data(file, path+'Hc_T.txt')
#Jc_H_data(file, path+'Jc_H.txt')
X_T(file, path)
#copper(file, path)
#M_H(file, path)