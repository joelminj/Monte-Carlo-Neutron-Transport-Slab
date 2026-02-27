# Monte Carlo Neutron Transport: Heterogeneous Water-Carbon Slab

This simulation models neutrons traveling through a multi-region slab consisting of Water ($H_2O$) and Graphite ($C_{12}$), employing energy-dependent cross-sections and advanced variance reduction techniques.

## 1. Project Overview

The simulation tracks the "history" of neutrons from a source at the origin ($x=0$) as they travel through:


**Region 1 (0-10 cm):** Water — an effective moderator/shield. 


**Region 2 (10-30 cm):** Graphite — relatively transparent to neutrons at the simulated energy levels. 



## 2. Theoretical Foundations

The simulation is based on the **Linear Boltzmann Equation (LBE)**, which describes the particle balance in phase space $(d^3r\, dE\, d\Omega)$. 

### Mathematical Principles

#### The Fundamental Formulation of Monte Carlo (FFMC)

To simulate random physical processes (like path length or scattering angle) on a computer, the simulation uses the FFMC to relate a random number $\eta \in [0, 1]$ to a random variable $x$ through its Cumulative Distribution Function (CDF), $P(x)$:


$$P(x) = \int_{a}^{x} p(x') \, dx' = \eta$$

#### Path Length Sampling

For a particle moving through a medium with total macroscopic cross-section $\Sigma_t$, the distance $d$ to the next collision is sampled using:


$$d = -\frac{\ln(\eta)}{\Sigma_t}$$

#### Elastic Scattering Kinematics

When a neutron of energy $E$ scatters off a nucleus of mass number $A$, the energy of the neutron after scattering $E'$ is determined by the scattering angle in the CM frame $\mu_{cm}$:


$$E' = E \frac{(1 + \alpha) + (1 - \alpha) \mu_{cm}}{2}$$

where $\alpha = \left( \frac{A-1}{A+1} \right)^2$. 

#### Flux Estimation

The simulation employs the **Track Length Estimator** to calculate the scalar flux $\Phi$ in spatial bins:


$$\Phi \approx \frac{1}{V \cdot N} \sum w \cdot d_{step}$$

where $w$ is the particle weight, $d_{step}$ is the segment length in the bin, and $N$ is the total number of source particles.

### 3. Computational Geometry & Kinematics

### Geometry Configuration and Surface Crossing
The simulation models a 3D bounding box for the heterogeneous slab. While neutrons move in 3D space $(x, y, z)$, the material heterogeneity is strictly 1D along the x-axis.

Dimensions: $X \in [0, 30]$ cm, $Y \in [-10, 10]$ cm, $Z \in [-15, 15]$ cm.
Interface: A discrete material boundary exists at $x = 10$ cm separating Water and Graphite.
To prevent floating-point precision errors (infinite logic loops) when a particle lands exactly on a boundary, the simulation employs a surface crossing nudge. If the distance to the interface $d_{interface}$ is the minimum sampled distance, the particle is advanced to the boundary and nudged slightly into the new material:

$$x_{new} = x_{interface} + (u \cdot 10^{-6})$$

where $u$ is the direction cosine along the x-axis. Energy and direction remain unchanged during a pure surface crossing.

### Source Definition (Initial Phase Space)
Every neutron history is initialized with a specific set of parameters that define its starting state. In this simulation, the source is modeled as a monoenergetic, isotropic point source emitting into the slab:

Position ($\vec{r}_0$): Neutrons are injected exactly at the origin, representing the entry boundary of the water region:

$$\vec{r}_0 = (0, 0, 0)$$

Energy ($E_0$): The source is monoenergetic, representing fast neutrons (e.g., from a fission spectrum) prior to any moderation:

$$E_0 = 2.0 \text{ MeV}$$

Direction ($\vec{\Omega}_0$): Neutrons are emitted isotropically. To achieve a uniform distribution over the unit sphere, the polar angle cosine ($\mu$) and azimuthal angle ($\phi$) are sampled using uniform random numbers $\xi_1, \xi_2 \sim \mathcal{U}(0,1)$:

$$\mu = 2\xi_1 - 1$$

$$\phi = 2\pi\xi_2$$

The initial direction cosines $(u_0, v_0, w_0)$ are then calculated as:

$$u_0 = \sin(\arccos(\mu)) \cos(\phi)$$

$$v_0 = \sin(\arccos(\mu)) \sin(\phi)$$

$$w_0 = \mu$$

Statistical Weight ($w_0$): Because the LBE is solved using implicit capture, every neutron is initialized with a statistical weight of exactly 1:

$$w_0 = 1.0$$

### Target Nucleus Selection
When a collision occurs in a compound material like Water ($H_2O$), the simulation must stochastically determine which isotope the neutron scatters off. This is sampled using the ratio of their macroscopic scattering cross-sections. A uniformly distributed random number $\xi \sim \mathcal{U}(0,1)$ is compared against the scattering probability of Hydrogen:

$$P(H) = \frac{\sigma_{s,H}(E) \cdot N_H}{\sigma_{s,H}(E) \cdot N_H + \sigma_{s,O}(E) \cdot N_O}$$

If $\xi < P(H)$, the collision is processed with Hydrogen ($A=1$); otherwise, it is processed with Oxygen ($A=16$).

### Elastic Scattering Kinematics and 3D Rotation
Scattering is assumed to be isotropic in the Center-of-Mass (CM) frame. The CM scattering angle cosine, $\mu_{cm}$, is sampled uniformly: $\mu_{cm} = 2\xi - 1$.

To continue tracking the particle in the Laboratory (Lab) frame, we must map $\mu_{cm}$ to the Lab scattering angle cosine, $\mu_{lab}$:

$$\mu_{lab} = \frac{1 + A \mu_{cm}}{\sqrt{1 + A^2 + 2A \mu_{cm}}}$$

Because the neutron is traveling in 3D space with a pre-collision direction vector $\vec{\Omega} = (u, v, w)$, we cannot simply add $\mu_{lab}$. We must perform a 3D Directional Rotation. Given the scattering polar angle $\theta = \arccos(\mu_{lab})$ and a randomly sampled azimuthal angle $\phi = 2\pi\xi$, the new direction cosines $(u', v', w')$ are calculated using a rotation matrix:

$$u' = u \mu_{lab} + \frac{\sqrt{1-\mu_{lab}^2} \cdot (u w \cos\phi - v \sin\phi)}{\sqrt{1-w^2}}$$

$$v' = v \mu_{lab} + \frac{\sqrt{1-\mu_{lab}^2} \cdot (v w \cos\phi + u \sin\phi)}{\sqrt{1-w^2}}$$

$$w' = w \mu_{lab} - \sqrt{1-\mu_{lab}^2} \cdot \sqrt{1-w^2} \cdot \cos\phi$$

(Note: If the particle is traveling nearly parallel to the Z-axis, i.e., $1 - |w| < 10^{-5}$, a simplified coordinate rotation is applied to avoid division by zero). Vector normalization is enforced post-rotation to prevent numerical drift.

### Event Selection Logic
At each step, the neutron faces three potential outcomes. The distance to each is calculated:

$d_{coll} = -\frac{\ln(\xi)}{\Sigma_t(E)}$ (Distance to next collision)
$d_{bound}$ (Distance to external leakage boundary)
$d_{interface}$ (Distance to material interface)

The actual distance traveled by the particle is the minimum of these three:

$$d_{step} = \min(d_{coll}, d_{bound}, d_{interface})$$

The spatial tally (Track Length Estimator) is updated using $d_{step}$, and the particle's physical state is updated based on which event won the race.


## 5. Variance Reduction Techniques

To improve computational efficiency and reduce the variance of the estimators, the LBE is solved using implicit techniques:

1. **Implicit Capture (Survival Biasing)**: Instead of terminating a particle upon absorption, its weight $w$ is reduced by the survival probability (albedo) at every collision:

$$w_{new} = w_{old} \times \left( \frac{\Sigma_s(E)}{\Sigma_t(E)} \right)$$

2. **Russian Roulette**: To prevent the CPU from wasting cycles tracking negligible particle weights, Russian Roulette is triggered when $w < 10^{-4}$. The particle is given a 10% chance ($P=0.1$) to survive. If it survives, its weight is boosted to $w_{new} = w_{old} / 0.1$ to maintain a "fair game" and unbiased tally. If it fails, the history is terminated.

3. **Energy Cutoff**: Neutrons thermalize as they lose energy. To bound the simulation execution time, histories are forcefully terminated when the particle's energy drops below a cutoff threshold of $10^{-5}$ MeV (10 eV).




## 4. Implementation Details

The core logic is structured into several modular Python classes:

* `NuclearConstants`: Defines atomic densities and molar masses for $H_2O$ and $C_{12}$. 
* `CrossSectionModels`: Provides energy-dependent macroscopic cross-sections. 
* `HeterogeneousGeometry`: Manages material boundaries and distance-to-surface calculations. 
* `HeterogeneousMCSimulation`: The primary engine that executes the random walk, manages tallies, and computes batch statistics. 



### Requirements

* `numpy`
* `matplotlib`
* `dataclasses`
