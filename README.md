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
