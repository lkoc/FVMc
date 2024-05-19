import numpy as np
import matplotlib.pyplot as plt

# Material properties: density (rho), specific heat capacity (Ce), and thermal conductivity (k)
materials = {
    "Aluminum": {"rho": 2712, "Ce": 897, "k": 237},
    "Copper": {"rho": 8960, "Ce": 385, "k": 401},
    "XLPE": {"rho": 940, "Ce": 2300, "k": 0.4},
    "Typical Soil1": {"rho": 2000, "Ce": 800, "k": 0.2},
    "Typical Soil2": {"rho": 1500, "Ce": 900, "k": 0.33}
}

# Space step range from 0.001 m to 0.032 m
h = np.linspace(0.001, 0.32, 100)

# Calculating delta t_max for each material
#delta_t_max = {material: (h**2 * props["rho"] * props["Ce"]) / (4 * props["k"]) 
#               for material, props in materials.items()}

delta_t_max = {material: (2*h* props["rho"] * props["Ce"]) / (4 * props["k"]) 
               for material, props in materials.items()}


# Plotting with logarithmic axes
plt.figure(figsize=(10, 6))
for material, dt_max in delta_t_max.items():
    plt.loglog(h, dt_max, label=material)

plt.xlabel("Space step (h) [m]")
plt.ylabel("Maximum time step ($\Delta t_{max}$) [s]")
plt.title("Maximum Time Step as a Function of Space Step (Logarithmic Scale)")
plt.legend()
plt.grid(True, which="both", linestyle='--')
plt.show()

'''
To calculate the maximum time step (\( \Delta t_{\text{max}} \)) for a
 transient random walk for the heat equation for materials like 
 aluminum, copper, XLPE, and typical soil, we'll need to know their 
 respective properties. The properties include density (\( \rho \)), 
 specific heat capacity (\( C_e \)), and thermal conductivity (\( k \)). 
 
 The formula you've given is:

\[ \Delta t_{\text{max}} = \frac{h^2 \cdot \rho \cdot C_e}{4 \cdot k} \]

We'll calculate \( \Delta t_{\text{max}} \) for a range of
 space steps (\( h \)) from 0.001 m to 0.032 m for each material and plot the
 results. Let's define the properties for each material:

- Aluminum (already provided in your question)
  - \( \rho = 2712 \, \text{kg/m}^3 \)
  - \( C_e = 897 \, \text{J/(kg.K)} \)
  - \( k = 237 \, \text{W/(m.K)} \)
  
- Copper
  - \( \rho = 8960 \, \text{kg/m}^3 \)
  - \( C_e = 385 \, \text{J/(kg.K)} \)
  - \( k = 401 \, \text{W/(m.K)} \)
  
- XLPE (Cross-linked polyethylene)
  - \( \rho = 940 \, \text{kg/m}^3 \)
  - \( C_e = 2300 \, \text{J/(kg.K)} \) (assuming a typical value)
  - \( k = 0.4 \, \text{W/(m.K)} \) (assuming a typical value)
  
- Typical Soil
  - \( \rho = 1500 \, \text{kg/m}^3 \) (assuming a typical value)
  - \( C_e = 840 \, \text{J/(kg.K)} \) (assuming a typical value)
  - \( k = 1.5 \, \text{W/(m.K)} \) (assuming a typical value)
  
Let's proceed with the calculations and the plot.

The graph above shows the maximum time step (\( \Delta t_{\text{max}} \)) as
 a function of the space step (\( h \)) for aluminum, copper, XLPE, 
 and two typical soils. 
 Each line represents a different material, illustrating how 
 \( \Delta t_{\text{max}} \) varies with \( h \) for each. 
 This visual representation aids in understanding the impact of material 
 properties on the stability condition for the heat equation's numerical 
 solution.
'''