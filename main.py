"""
Contains a mathematical model that describes the shunting affect of GABA receptors
on presynaptic terminal through modulation of calcium influx

Program has interface, which contains:
- form to enter parameter
- button to start simulation
- output of simulation time
- graphics that represents results by plotting resulting potential from axon length
  and from time in observation point
"""

import tkinter as tk
from tkinter import ttk

import matplotlib.pyplot as plt
import numpy as np
import math

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import time


def gauss_elimination(a, b, d, f: np.array) -> np.array:
    """
    Thomas algorithm for tridiagonal matrix
    Args:
        a: diagonal coefficient
        b: coefficients below the diagonal
        d: coefficients above the diagonal
        f: right side of the matrix
    """
    length = len(f)
    d[0] = d[0] / b[0]
    f[0] = f[0] / b[0]
    b[0] = 1.0

    for i in range(1, length):
        b[i] = b[i] - a[i] * d[i - 1]
        f[i] = (f[i] - a[i] * f[i - 1]) / b[i]
        d[i] = d[i] / b[i]
        b[i] = 1.0

    for i in range(1, length):
        f[length - i - 1] = f[length - i - 1] - f[length - i] * d[length - i - 1]

    return f


class EquationCoeffs:
    """
    Class contains equations to calculate special variable for define sodium and
    potassium conductance (through coefficients m, n and h)
    """

    def __init__(self, V: float) -> None:
        Vrev_na = -50  # reverse potential for potassium channels

        self.alpham = 0.1 * (V + 35) / (1.0 - math.exp(-0.1 * (V + 35)) + 1e-9)
        self.betam = 4.0 * math.exp(-0.0556 * (V - Vrev_na))
        self.alphah = 0.07 * math.exp(-(V - Vrev_na) / 20)
        self.betah = 1.0 / (1 + math.exp(-0.1 * (V - Vrev_na)) + 1e-9)
        self.alphan = 0.01 * (V + 50) / (1 - math.exp(-0.1 * (V - Vrev_na)) + 1e-9)
        self.betan = 0.125 * math.exp(-(V + 60) / 80)


class Fields:
    """
    Class contains functions to deal with fields to use in next time step

    Parameters:
        V: electric potential of neuronal membrane
        m, h, n: coefficients for calculate sodium and potassium conductance
    """

    def __init__(self, V, m, h, n) -> None:
        """ """
        self.V = V
        self.m = m
        self.h = h
        self.n = n

    def save_coeffs_fields(self, m, h, n):
        """Rewrite coefficients m, h and n"""
        self.m = m
        self.h = h
        self.n = n

    def save_V_field(self, V):
        """Rewrite potential field"""
        self.V = V

    def get_fields(self):
        """Get potential and m, h, n coefficients"""
        return self.V, self.m, self.h, self.n


class Simulation:
    """
    Class contains calculation of initial parameters and fields and
    calculation of one time step of neuron simulation
    """

    def __init__(self) -> None:

        self.__params_init__()  # init all static parameters
        self.__grid_init__()  # init grid parameters

        # Start recording from time point
        self.point = int(self.length // self.dx - self.shunt_length // self.dx)
        self.__field_init__()  # init field data

        # Time step in seconds
        self.time = np.array(
            [self.dt * j for j in range(int(self.impulse_time // self.dt))]
        )
        self.point = self.grid_size - 1

    def __params_init__(self) -> None:
        """Initialize static parameters"""
        self.Vna = 50  # reverse potential for sodium channels {mV}
        self.Vk = -77  # reverse potential for potassium channels {mV}
        self.gnap = 120  # sodium conductance constant {mS/cm2}
        self.gkp = 36  # potassium conductance constant {mS/cm2}
        self.capacity = 1  # membrane capacity {mkF/cm2}
        self.gl = 0.3  # leakage conductance mS/cm2
        self.Vl = -50  # membrane current reverse potential {mV}
        self.Vus = -40  # -65 synaptic current reverse potential {mV}
        self.Vca = 135  # reverse potential for calcium channels {mV}
        self.gcap = 4  # calcium conductance constant {mS/cm2}
        self.dt = 0.05  # time step {ms}
        self.radius = 1e-4  # axon radius {cm}
        self.rm = 10  # membrane resistivity {kOm*cm2}
        self.impulse_time = 30  # 100  # Time inside programm in sec
        self.time_steps = int(self.impulse_time / self.dt + 1)  # quantity of time steps
        self.ri = 0.1 / self.radius / self.radius
        # "ri" - linear resistance of the intracellular fluid {kOm/cm}
        # 0.1 is cylindrical resistance {kOm*cm}

    def __grid_init__(self) -> None:
        """Initialize of grid parameters"""
        self.grid_size = 3160  # Number of grid nodes
        self.length = 0.1  # axon length {cm}
        self.shunt_length = 2e-3  # shunt length {cm}
        self.dx = self.length / (self.grid_size - 1)  # Step for x coordinate {cm}
        self.x = np.array(
            [
                (i - 1) * self.length / (self.grid_size - 1)
                for i in range(self.grid_size)
            ]
        )  # x coordinates field {cm}
        self.N_shunt = int(
            self.shunt_length // self.dx
        )  # Number of grid nodes for shunt

    def __field_init__(self) -> None:
        """Initialize fields before calculations"""
        self.s = np.zeros(self.grid_size)

        V_initial = np.array([-61.5] * self.grid_size)  # initial change of potential
        m = np.zeros(self.grid_size)
        n = np.zeros(self.grid_size)
        h = np.zeros(self.grid_size)

        for i in range(self.grid_size):

            # Define synaptic conductance ans initial potential on shunt
            if (i > self.point - self.N_shunt) and (i < self.point + self.N_shunt):
                # self.s[i] = self.conductance_value #7.646395  #{mS/cm2}
                V_initial[i] = -228061 * self.x[i] ** 2 + 41184 * self.x[i] - 1905.1

            # Define innitial potantial field
            # Equations choosed to reduce field settling time
            if i > self.point + self.N_shunt:
                V_initial[i] = (
                    -891373 * self.x[i] ** 3
                    + 304171 * self.x[i] ** 2
                    - 34101 * self.x[i]
                    + 1210.4
                )
            elif i < self.point - self.N_shunt:
                V_initial[i] = (
                    981824 * self.x[i] ** 4
                    - 118408 * self.x[i] ** 3
                    + 5398.2 * self.x[i] ** 2
                    - 75.642 * self.x[i]
                    - 60.999
                )

            coeffs = EquationCoeffs(V_initial[i])

            m[i] = coeffs.alpham / (coeffs.alpham + coeffs.betam)
            h[i] = coeffs.alphah / (coeffs.alphah + coeffs.betah)
            n[i] = coeffs.alphan / (coeffs.alphan + coeffs.betan)

        self.current_fields = Fields(V_initial, m, h, n)

    def conductance_field(self, conductance: float):
        """Initialize synaptic conductance from enter in interface"""
        for i in range(self.grid_size):

            # Define synaptic conductance and on shunt
            if (i > self.point - self.N_shunt) and (i < self.point + self.N_shunt):
                self.s[i] = conductance  # 7.646395  #{mS/cm2}

    def _one_step_field(self):
        """Calculate fields on one time step to define coefficients for Thomas algorithm"""
        # Coefficients for tridiagonal matrix algorithm/Thomas algorithm
        a = np.zeros(self.grid_size)
        b = np.zeros(self.grid_size)
        d = np.zeros(self.grid_size)
        f = np.zeros(self.grid_size)

        # Coefficients for potassium channels equation
        m = np.zeros(self.grid_size)
        h = np.zeros(self.grid_size)
        n = np.zeros(self.grid_size)

        [V_previous_step, m_previous_step, h_previous_step, n_previous_step] = (
            self.current_fields.get_fields()
        )

        for i in range(self.grid_size):

            coeffs = EquationCoeffs(V_previous_step[i])

            m[i] = (coeffs.alpham + m_previous_step[i] / self.dt) / (
                coeffs.alpham + coeffs.betam + 1 / self.dt
            )
            n[i] = (coeffs.alphan + n_previous_step[i] / self.dt) / (
                coeffs.alphan + coeffs.betan + 1 / self.dt
            )
            h[i] = (coeffs.alphah + h_previous_step[i] / self.dt) / (
                coeffs.alphah + coeffs.betah + 1 / self.dt
            )

            gna = self.gnap * h[i] * m[i] ** 3
            gk = self.gkp * n[i] ** 4

            a[i] = 1.0 / self.radius / 2 / self.ri / self.dx / self.dx
            b[i] = (
                -1.0 / self.radius / self.ri / self.dx / self.dx
                - self.capacity / self.dt
                - gna
                - gk
                - self.gl
                - self.s[i]
            )
            d[i] = 1.0 / self.radius / 2 / self.ri / self.dx / self.dx
            f[i] = (
                -self.capacity * V_previous_step[i] / self.dt
                - gna * self.Vna
                - gk * self.Vk
                - self.gl * self.Vl
                - self.s[i] * self.Vus
            )

        # border conditions
        a[0] = 0
        b[0] = 1 / self.ri / self.dx
        d[0] = -1 / self.ri / self.dx
        f[0] = 0  # mkA
        a[self.grid_size - 1] = -1 / self.ri / self.dx
        b[self.grid_size - 1] = 1 / self.ri / self.dx
        d[self.grid_size - 1] = 0
        f[self.grid_size - 1] = 0

        self.current_fields.save_coeffs_fields(m, h, n)

        return a, b, d, f

    def time_step(self, step: int) -> None:
        """Calculate potential on one time step"""
        [a, b, d, f] = self._one_step_field()

        # Define initial impulse in time
        if (self.dt * step > 10) and (self.dt * step < 15):
            f[0] = 0.0002

        result = gauss_elimination(a, b, d, f)

        self.current_fields.save_V_field(result)


class Interface:
    """Class contains operation with tkinter interface"""

    def __init__(self, root: tk.Tk) -> None:
        plt.ion()
        self.control_frame = ttk.Frame(root)
        self.conductance_shunt = tk.StringVar(value=0)
        self.control_frame.pack(side=tk.TOP, fill=tk.X)

        ttk.Label(self.control_frame, text="Synaptic conductance on shunt").pack(
            side=tk.LEFT
        )
        conductance_shunt_entry = ttk.Entry(
            self.control_frame, textvariable=self.conductance_shunt
        )
        conductance_shunt_entry.pack(side=tk.LEFT)

        self.label = ttk.Label(text="Time: 0.00 s", font=("Arial", 14))
        self.label.pack(side=tk.TOP)

        self.figure, self.ax = plt.subplots(1, 2, figsize=(14, 4.5), dpi=100)
        self.figure.canvas = FigureCanvasTkAgg(self.figure, root)
        self.figure.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.ax[0].set_xlabel("x, cm")
        self.ax[0].set_ylabel("Potential. mV")
        self.ax[1].set_xlabel("time, s")
        self.ax[1].set_ylabel("Potential, mV")
        

    def define_start_button(self, function: callable):
        """Define button, that start simuation"""
        start_button = ttk.Button(
            self.control_frame,
            text="Start Simulation",
            command=lambda: function(self),
        )
        start_button.pack(side=tk.LEFT)

    def read_conductance(self) -> float:
        """Get conductance value from tkinter frame"""
        return float(self.conductance_shunt.get())


def start(interface: Interface):
    """
    Интерфейс уже содержит в себе необходимые переменные, их нужно переделать

    Args:
        figure: The first parameter.
        canvas: The second parameter.
    """
    simulation = Simulation()

    conductance = interface.read_conductance()
    simulation.conductance_field(conductance)
    plot_time = []
    plot_V = []
    interface.ax[0].clear()
    interface.ax[1].clear()

    for step in range(simulation.time_steps):
        simulation.time_step(step)
        interface.label.config(text=f"Time: {np.round(simulation.dt*step, 2):0<5} s")
        interface.ax[0].clear()
        interface.ax[0].plot(simulation.x, simulation.current_fields.V)
        interface.ax[0].set_xlabel("x, cm")
        interface.ax[0].set_ylabel("Potential, mV")

        plot_time.append(np.round(simulation.dt * step, 2))
        plot_V.append(simulation.current_fields.V[simulation.point])
        interface.ax[1].clear()
        interface.ax[1].plot(plot_time, plot_V, color="green")
        interface.ax[1].set_xlabel("time, s")
        interface.ax[1].set_ylabel("Potential, mV")

        interface.figure.canvas.flush_events()
        interface.figure.canvas.draw()
        time.sleep(0.02)


if __name__ == "__main__":
    root = tk.Tk()
    root.title("Action Potential Simulation")

    interface = Interface(root)
    interface.define_start_button(start)

    root.mainloop()
