#### Continent-content collision: 
##### Model parameters
* Model duration: µ = 50 Ma
* Horizons to track: 10 
* Horizon depths: [5, 8, 11, 13, 15, 17, 19, 21, 23, 25]

##### Lithology and heat : 
* Basement thermal conductivity: µ = 2.5
* Surface heat production: µ = 2
* Heat production depth: µ = 45 Km
* Surface heat flow: µ = 100 > (hp_dep * surf_hp)
* Number of heat flow variations (INT): 0
* Time of heat flow at change: Na
* Value of heat flow at change: Na
* Temperature at upper surface boundary (ºC):

##### Exhumation history
* Uplift rate: uniform distribution of erosion rates between 0 and 5 mm/yr
* Num uplift/erosion periods (INT): µ = 3
* Beginning of uplift period: 
* Duration of uplift period: 
* Total uplift during episode

# Lithology: 

Assumptions: constant thermal conductivity. Lithologically consistent. 

### Basement thermal conductivity (W/m*K): b_cond
Lithology dependent. Represents the thermal conductivity of the base of the model. Gabbros and granites range from 2-3.5. 

### Surface heat production (µW/$m^3$): s_hp
* Old continental crust: 1.5
* Young continental crust: 2
* Gabbros: .2

### Heat production depth (Km): hp_dep
Bimodal distribution. Archean crust around 35km, Proterozoic crust is 45km [(Source: Tewari et al.)](https://www.sciencedirect.com/topics/earth-and-planetary-sciences/crustal-thickness#:~:text=Abstract,2.7%20g%2Fcm3). 

### Surface heat flow (mW/$m^2$): hf_surf
Mean hf_surf by geological type ([Source: Global Heat Map](http://heatflow.org/province/?data_type=heat_flow&group_by=type#)): 
* Collisional orogen: 63
* Continental collisional orogen: 100
* Continental orogen: 89
* Fold and thrust belt: 65
* Orogen: 68
* Orogenic belt: 60
* Orogenic plateau: 61
* Subduction orogen: 88

### Number of heat flow variations:
Between 0 and 2: {d, i, di, id, ii, dd}

### Time of heat flow value at change (Ma):
Dependent on model time. 

### Value of heat flow at change(mW/$m^2$): see hf_surf