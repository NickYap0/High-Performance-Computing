# Parallelisation of 2D Shallow Water Propagation Simulation

The objective of this coursework assignment is to write a parallel numerical code to solve the 2D shallow-water
equations. These partial differential equations (PDE) can be used to model, for example, the propagation of tsunamis.

$$
\begin{align}
    \frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} + g \frac{\partial h}{\partial x} = 0\\
    \frac{\partial v}{\partial t} + u \frac{\partial v}{\partial x} + v \frac{\partial v}{\partial y} + g \frac{\partial h}{\partial y} = 0\\
    \frac{\partial h}{\partial t} + u \frac{\partial hu}{\partial x} + v \frac{\partial hv}{\partial y} = 0\\
\end{align}
$$

where $u(x, y)$ and $v(x, y)$ are the $x$- and $y$-components of velocity, $h(x, y)$ is the surface height, and $g$ is the
acceleration due to gravity. For simplicity, we have neglected friction, Coriolis and viscous forces.


# Partial Derivatives

Finite difference method is adopted in order to numerically solve this problem. For spatial derivative terms, a 6th order central difference stencil was used, the algorithm is shown as below:

$$
\frac{\partial u}{\partial x} \approx \frac{1}{\Delta x}\Big(\frac{1}{60}u_{i-3} + \frac{1}{20}u_{i-2}-\frac{3}{4}u_{i-1}+\frac{3}{4}u_{i+1}-\frac{3}{20}u_{i+2}+\frac{1}{60}u_{i+3}\Big)
$$

Similar expressions can be constructed for other derivatives. The PDE should be discretised on a periodic rectan-
gular grid of size $Nx\times Ny$ with $dx = dy = 1$ and a domain of $[0, Nx]\times[0, Ny]$.


For the time-integration we will use 4th-order Runge-Kutta explicit scheme. This can be formulated as:

$$
y_{n+1} = y_n + \frac{1}{6}(k_1+2k_2+2k_3+k_4)\Delta t
$$

with 

$$
\begin{align}
    k_1 = f(y_n)\\
    k_2 = f(y_n + \Delta tk_1/2)\\
    k_3 = f(y_n + \Delta tk_2/2)\\
    k_4 = f(y_n + \Delta tk_3)\\
\end{align}
$$


## Evaluation of $f$
Rearrange Equations (1)-(3) to put all but the time derivative terms on the right-hand side. The function f
above is therefore used to evaluate the discrete approximation to these time derivatives by computing these right-
hand side terms in discrete form. The bulk of the operations are accounted for by the numerical approximations to
the derivatives. In this simulation, the user will be able to enter command line input (--choice) to choose whether to construct the spatial derivative terms using a $\textit{for-loop appraoch}$(1) or $\textit{matrix-vector multiplication}$(2) (BLAS)

## Initial and Boundary Conditions
The simulation uses the following initial conditions:

$$
\begin{align}
u(x,y,0) = v(x,y,0) = 0\\
h(x,y,0) = h_0(x,y)
\end{align}
$$

where $h_0(x, y)$ is a function prescribing the initial surface height (four test cases included).


Periodic boundary conditions should be used on all boundaries, so that waves which propagate out of one of the
boundaries re-enter on the opposite boundary. This means that, for example, with $dx = 1$ and $Nx = 100$, the last
point is at $x = 99$, even though the domain is of length 100 in $x$. This is because the point at $x = 100$ is the same as the point at $x = 0$.

## Command-line input
The repository include a Makefile we make ease when compiling the program. It also incorporates parallelisation using openMP with highest level optimise -Wall -o3. A sample command line input is as such,

```command line
--dt 0.1 --T 80.0 --Nx 100 --Ny 100 --ic 1 --np 9 --calc 1
```

where $--np 9$ sets the environment of the number of threads computing. 
