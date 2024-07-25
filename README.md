# action_potential_simulation
Mathematical model of action potential propagation with a study of the possibility of changing the potential amplitude under the influence of calcium influx

### Physics 
Information in our nerve system is transmitted by electric impulses - action potential - from one axon to another. Experiments show that calcium influx may cause changes in membrane potential. We use mathemtical model on base of cabel equation to explore the influence of shunting part with calcium on action potential propogation
The cabel equation 
```math
 \frac{\partial V }{\partial t} = \frac{1}{2r_i}\frac{\partial^2 V}{\partial t^2} - I_a - s(x, t)(V - V_s)
```
where $V$ is the membrane potential, $C$ is the membrane capacity, $I_a$ is the ion current, $r_i$ is the linear resistance of the intracellular fluid, $s$ is the synaptic conductance on shunt, $V_s$ is the synaptic current reversal potential
Thomas algorithm-based implicit scheme has been applied

If you interested in neuronal effects, biophysical results and more details of this experiment you may check the published paper: https://link.springer.com/chapter/10.1007/978-3-031-19032-2_23

### Required packages
- math
- matplotlib
- numpy
- time
- tkinter

### Functionality
* Simulation 
    * Initialization of all parameters and fields of variable
    * Calculation for each time step
* Interface
  * Description of all window elements
    * form to enter parameter
    * button to start simulation
    * output of simulation time
    * graphics that represents results by plotting resulting potential from axon length
and from time in observation point
* Fields
  * Operations with membrane potential and coeffisients fields. Thomas algorithm use those fields from previous time step to calculate its on new step
* EquationCoeffs
  * Equations to calculate special variable for define sodium and potassium current
* gauss_elimination
  * Thomas algorithm for tridiagonal matrix

### Interface instructions
When you run program, you will see the interface like in picture below. 
* Firstfull, you need to enter the main parameter, that is the case of all this simulation, - the synaptic conductance on shunt
    * You may rest it equal zero or change to value between 0 and 10, that's enough for notable change in graphs
* Next, you need to press the button 'Start simulation'
Lines will be change every time step, showing distrubution of membrane potential in mV along the axon and change of membrane potential in time in observation point on shunt(near the right enc of the axon)

![](/assets/images/interface.png)
