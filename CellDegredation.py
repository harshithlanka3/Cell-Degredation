import math
# makes calculations and life easier


# f(T) - Temperature dependent part of the function
kCalRef = 3.694e-4      # Reference stress factor value at T = 298.15K and SOC = 50% in h^-0.5
EaCal = 20592           # Fitting stress factors to SOC = 100% in J / mol
Rg = 8.314              # Universal gas constant in J / (mol K)
Tref = 298.15           # Reference temperature in K

def temperature_function(T): # T in K
    exp = (-EaCal / Rg) * ((1 / T) - (1 / Tref))
    return math.e ** exp

# f(SOC) - SOC dependent part of the function
alpha = 0.384           # constant in Tafel equation
F = 96485               # Faraday's constant in C / mol
UaRef = 0.123           # Reference potential at SOC 50% in V
k0 = 0.142              # Constant offset

def anode_stochiometry(soc):
    return 8.5e-3 + (soc * (7.8e-1 - 8.5e-3))


def open_circuit_potential_data(soc): # soc in decimal. Ua() in article.
    return 0.6379 + 0.5416 * (math.e ** (-305.5309 * anode_stochiometry(soc))) + \
           (0.044 * math.tanh(-(anode_stochiometry(soc) - 0.1958)/0.1088)) - \
           (0.1978 * math.tanh((anode_stochiometry(soc) - 1.0571)/0.0854)) - \
           (0.6875 * math.tanh((anode_stochiometry(soc) + 0.0117)/0.0529)) - \
           (0.0175 * math.tanh((anode_stochiometry(soc) - 0.5692)/0.0875))

def stateOfChargeFunction(soc): # SOC in decimal
    exp = ((alpha * F) / Rg) * ((UaRef - open_circuit_potential_data(soc)) / Tref)
    return (math.e ** exp) + k0

def stressFactor(T, soc): # T in K and SOC in decimal. Calculating KCal
    return kCalRef * temperature_function(T) * stateOfChargeFunction(soc)


#-----Final Function-----#

def calendar_aging_capacity_loss(T, soc, time): # T in K and SOC in decimal and time in hours.
    return stressFactor(T, soc) * (time ** 0.5)