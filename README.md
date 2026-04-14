# 🦅 Aero: NACA Airfoil Wind Tunnel

An interactive aerodynamic simulation app built with **Python** and **Streamlit**, using the **Lattice Boltzmann Method (LBM)** to simulate airflow around a NACA 0012 airfoil.

This project combines **computational fluid dynamics**, **numerical methods**, and **interactive visualization** in a compact engineering-focused dashboard.

---

## Overview

The application simulates 2D flow around a wing profile in a virtual wind tunnel. Users can adjust parameters such as:

* angle of attack
* flow velocity
* viscosity parameter

The simulation then visualizes the resulting flow field and estimates basic aerodynamic quantities such as:

* lift coefficient
* drag coefficient
* lift-to-drag ratio

---

## Features

* Interactive Streamlit dashboard
* Real-time flow simulation using the Lattice Boltzmann Method
* NACA 0012 airfoil geometry generation
* Adjustable angle of attack
* Adjustable inflow velocity and viscosity
* Vorticity field visualization
* Streamline plotting
* Estimated lift and drag metrics

---

## Method

The simulation is based on a **D2Q9 Lattice Boltzmann model**, a common approach for simulating fluid flow on a discrete grid.

### Main components

**1. LBM solver**

* Initializes density and velocity fields
* Performs streaming and collision steps
* Applies inlet, outlet, and no-slip boundary conditions
* Computes force contributions on the airfoil

**2. Airfoil geometry**

* Generates a NACA 0012 profile from its thickness distribution
* Rotates the geometry according to the chosen angle of attack
* Maps the airfoil onto the simulation grid

**3. Visualization**

* Displays the vorticity field as a background contour map
* Overlays streamlines to show flow behavior
* Computes approximate aerodynamic performance metrics

---

## Tech Stack

* Python
* NumPy
* Matplotlib
* SciPy
* Streamlit

---



## How to Run

1. Install dependencies:

```bash
pip install -r requirements.txt
```

2. Start the Streamlit app:

```bash
streamlit run app.py
```

---

## Example Parameters

The dashboard allows the user to experiment with:

* **Angle of Attack:** from -15° to 20°
* **Viscosity (Tau):** from 0.51 to 1.0
* **Flow Velocity:** from 0.02 to 0.15

These inputs change the airflow pattern and the estimated aerodynamic coefficients in real time.

---

## What This Project Demonstrates

This project highlights experience in:

* numerical simulation
* scientific computing
* fluid mechanics modeling
* interactive engineering tools
* translating theory into a working application

---

## Possible Improvements

Future extensions could include:

* more stable boundary conditions
* pressure field visualization
* additional airfoil profiles
* comparison between multiple geometries
* improved force estimation and validation against theory

---


