# Traceback
General code for performing trajectory traceback from AIU traceback code.

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
3441732292729818752, 195633325090480896

| **Object**        | **Proper Motion in RA (μα*)** | **Proper Motion in Dec (μδ)** |
|-------------------|------------------------------|------------------------------|
| **HD 37424 (DR2?)**       | 10.0 ± 0.8                   | −5.9 ± 0.6                   |
| **HD 37424 (DR3)**       | 12.228 ± 0.036                   | -9.407 ± 0.017                   |
| **Pulsar (old)**         | −24.4 ± 0.1                  | 57.2 ± 0.1                   |
| **Pulsar (new)**         | −24.4 ± 0.1                  | 57.2 ± 0.1                   |


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
config.orbits = 100_000  # This updates the "Orbits" key in the [Simulation] section
config.steps = 1000      # This updates the "Steps" key
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
