# Traceback
General code for performing trajectory traceback from AIU traceback code.

## Project 1: SNR G166.0+04.3
Green's catalog: https://www.mrao.cam.ac.uk/surveys/snrs/snrs.data.html

on vizier table "VII/272/snrs": https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=VII/272

### Summary SNR G166.0+04.3

https://www.mrao.cam.ac.uk/surveys/snrs/snrs.G166.0+4.3.html

| **Parameter**               | **Value**                                |
|-----------------------------|------------------------------------------|
| **Right Ascension**         | 05 26 30                                 |
| **Declination**             | +42 56                                   |
| **Size (arcmin)**          | 55×35                                    |
| **Type**                    | S                                        |
| **Flux Density at 1 GHz (Jy)** | 7                                        |
| **Spectral Index**          | 0.37                                     |
| **Radio**                   | Two arcs of strikingly different radii   |
| **Optical**                 | Nearly complete ring                     |
| **X-ray**                   | Predominantly in SW                      |
| **Distance**                | HI indicates 4.5 kpc, optical extinction suggests 3.2 kpc |

## Stars in question:

195633320791780608, 195633325090480896

## Stars in the region:

Region is defined as 1 deg around the center at: 

81.625 deg, 42.93333333 deg
has 130,096 sources

and distance range (with RUWE < 1.4 and PlxQ > 5):

2.5 kpc to 5 kpc
has 4,749 sources







## Workflow:
Aim: to trace back any 2 Gaia stars.
### Step 1: Create the input.tsv file
This file conatins the two stars to be traced back. The Gaia source IDs are required for the two stars.
Using the commands:

```python
trace = Traceback(195633320791780608, 195633325090480896)
trace.create_input_file()
```
### Step 2: Create the trace.conf file
This file contains the parameters of the traceback simulation. This requires a sample_config.conf file containing the main fields for the traceback.

Use the class ```SimulationConfig``` to create it:

```python
# Usage example
sample_config = 'sample_trace.conf'
config = SimulationConfig(sample_config)

# Change parameters with attributes
config.orbits = 100_000  
config.steps = 1000     
config.stepsize = -100
config.width = 100
config.star1 = '"two_trace/input.tsv#2"'
config.star2 = '"two_trace/input.tsv#3"'
config.outfile = '"two_trace/trace/test"'
# Save the changes
config.save_config('giessler_traceback/two_trace/trace.conf')
```

### Step 3: Run the simulation
Now that the two files: input.tsv and trace.conf are created, we are ready to run the simulation.
To run the simulation:

```python
execute_command('cd giessler_traceback ; examples/trace/traceback two_trace/trace.conf')
```
## Project 2: SNR G180.0−1.7
Green's catalog: https://www.mrao.cam.ac.uk/surveys/snrs/snrs.data.html

on vizier table "VII/272/snrs": https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=VII/272

### Summary SNR 147
[B. Dinçel et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.448.3196D/abstract) 

[https://www.mrao.cam.ac.uk/surveys/snrs/snrs.G166.0+4.3.html](https://www.mrao.cam.ac.uk/surveys/snrs/snrs.G180.0-1.7.html)

| **Parameter**               | **Value**                                |
|-----------------------------|------------------------------------------|
| **Right Ascension**         | 05 39 00                                 |
| **Declination**             | +27 50                                    |
| **Size (arcmin)**          | 180                                     |
| **Type**                    | S                                        |
| **Flux Density at 1 GHz (Jy)** | 65                                        |
| **Spectral Index**          | varies                                    |
| **Radio**                   | Large faint shell, with spectral break   |
| **Optical**                 | Wispy ring                    |
| **X-ray**                   | Possible detection                     |
|**Point sources**            | Pulsar within boundary, with faint wind nebula|
| **Distance**                | Various observations suggest about 1.2 kpc |

## Stars in question:

HD 37424, PSR J0538+2817 

3441732292729818752, n/a

| **Object**        | **Proper Motion in RA (μα*)** | **Proper Motion in Dec (μδ)** |
|-------------------|------------------------------|------------------------------|
| **HD 37424 (DR2?)**       | 10.0 ± 0.8                   | −5.9 ± 0.6                   |
| **HD 37424 (DR3)**       | 12.228 ± 0.036                   | -9.407 ± 0.017                   |
| **Pulsar (old)**         | −24.4 ± 0.1                  | 57.2 ± 0.1                   |
| **Pulsar (new)**         | −24.4 ± 0.1                  | 57.2 ± 0.1                   |
