---
layout: post
title: Lecture 4
---

- [Strings](#strings)
  - [Equation of motion](#equation-of-motion)
  - [Solutions to the 1D wave equation](#solutions-to-the-1d-wave-equation)
    - [Standing sinusoid as sum of travelling sinusoids](#standing-sinusoid-as-sum-of-travelling-sinusoids)
    - [Travelling sinusoid as sum of standing sinusoids](#travelling-sinusoid-as-sum-of-standing-sinusoids)
    - [General travelling wave solutions for $y(x,t)$ for the 1D wave equation](#general-travelling-wave-solutions-for-yxt-for-the-1d-wave-equation)
  - [Energy for the 1D wave equation](#energy-for-the-1d-wave-equation)
    - [Energy balance for an infinitely long string](#energy-balance-for-an-infinitely-long-string)
  - [Initial Conditions](#initial-conditions)
    - [D'Alembert solution with prescribed initial conditions](#dalembert-solution-with-prescribed-initial-conditions)
    - [Fourier solution with prescribed initial conditions](#fourier-solution-with-prescribed-initial-conditions)
  - [Energy of a semi-infinite string. Boundary conditions. {#sec:WaveEquationsBCs}](#energy-of-a-semi-infinite-string-boundary-conditions-secwaveequationsbcs)
    - [D'Alembert solution with prescribed boundary condition](#dalembert-solution-with-prescribed-boundary-condition)
    - [Fourier solution with prescribed boundary condition](#fourier-solution-with-prescribed-boundary-condition)
  - [Strings of finite length](#strings-of-finite-length)
    - [Boundedness of the solutions {#sec:BoundWECnt}](#boundedness-of-the-solutions-secboundwecnt)
    - [D'Alembert solution for strings of finite length](#dalembert-solution-for-strings-of-finite-length)
    - [Fourier solution for strings of finite length](#fourier-solution-for-strings-of-finite-length)
    - [Boundary condition with nonzero finite impedance](#boundary-condition-with-nonzero-finite-impedance)
  - [Viscous Loss](#viscous-loss)
  - [Dispersion. Group Velocity](#dispersion-group-velocity)

# Strings

This chapter introduces the analysis of continuous systems. The case of
a taut string will be considered firts, described by the one-dimensional
wave equation. While simple, this equation serves as a model for various
other systems, including wave propagation in acoustic tubes, and is
generally regarded as the first-ever example of a partial differential
equation. Outside of its historical importance, the wave equation can be
approached by a vast amount of analytical and numerical techniques.
Exact solutions exists in the form of travelling wavefronts (the
D'Alembert solution), as well as in the form of standing wave patterns
(the modes). Equivalence of these solutions is proven by Fourier theory.
Numerically, the solution may be computed using direct numerical
simulation, such as finite differences. A discretisation of the
D'Alembert solution forms the core of digital waveguide theory, probably
the most popular physical modelling technique to date. The solution of
the spatial eigenvalue problem using the spatial modes allows to
describe the system using a set of uncoupled oscillators, and the
numerical techniques used in previous chapters can be directly applied.

## Equation of motion

We begin by assuming a lossless, infinitely thin string lying under
tension $T(x,t)$ along the $x$-axis, whose vertical displacement is
denoted $y$, as illustrated in Figure
[1.1](#fig:infinite_string){reference-type="ref"
reference="fig:infinite_string"}.

::: center
![A infinitely long string lies oriented lengthways along the $x$-axis,
with its motion confined to the transverse
$y$-axis.](Figures/infinite_string.pdf){#fig:infinite_string
width="6.5cm"}
:::

::: center
![A infinitesimal section of the infinite string.
](Figures/stringElement.pdf){#fig:infinitesimal_string
width="0.7\\linewidth"}
:::

We are now going to derive the dynamic equilibrium equations for the
element of string depicted in Figure
[1.2](#fig:infinitesimal_string){reference-type="ref"
reference="fig:infinitesimal_string"}. For the moment, we are assuming
that the displacement $y$, the deflection angle $\theta$ and the tension
$T$ are all functions of the spatial coordinate $x$. Restricting the
attention to the vertical direction, the infinitesimal vertical force
$df$ acting on the string element is given by the difference of the
force at the right end and the force at the left end, hence
$$\begin{equation}\label{eq:balance_eq}
df = T(x+dx,t) \sin\left( \theta(x+dx,t) \right) - T(x,t) \sin\left( \theta(x,t) \right) = dx \frac{\partial}{\partial x}(T \sin\theta)$$
In dynamic equilibrium, the infinitesimal force is equal to the
infinitesimal mass times the acceleration. The infinitesimal mass is
given by the volume density $\rho$ times the volume of the string's
element, $dm = ds \rho A$, where $A$ is the cross section of the string
(supposed constant here), and where $ds$ is the infinitesimal arclength,
$ds = dx \sqrt{1+(\partial y / \partial x)^2}$. Hence
$$\begin{equation}\label{eq:inf_mass}
df = ds \rho A \frac{\partial^2 y}{\partial t^2}\end{equation}$$
We can now eliminate
$df$ in [\[eq:balance_eq\]](#eq:balance_eq){reference-type="eqref"
reference="eq:balance_eq"} using
[\[eq:inf_mass\]](#eq:inf_mass){reference-type="eqref"
reference="eq:inf_mass"}, to get $$\begin{equation}\label{eq:WEtheta}
\rho A \frac{\partial^2 y}{\partial t^2}ds = \frac{\partial}{\partial x}(T \sin\theta)dx$$
What we need now are explicit expressions for the tension $T$ and the
angle $\theta$ as a function of the coordinate $x$. From Figure
[1.2](#fig:infinitesimal_string){reference-type="ref"
reference="fig:infinitesimal_string"}, the angle is obtained
geometrically as
$$\tan\theta = \frac{\partial y}{\partial x} \quad \rightarrow \quad \sin\theta = \frac{\frac{\partial y}{\partial x}}{\sqrt{1 + (\frac{\partial y}{\partial x})^2}}$$
In the limit of small displacements, both the deflection angle and the
displacement are small, so that
$\sin \theta  \approx \partial y / \partial x$, and $dx \approx ds$. One
may also assume that, under such conditions, the tension is uniform,
i.e. $T(x,t)=T_0$. Thus, using these in
[\[eq:WEtheta\]](#eq:WEtheta){reference-type="eqref"
reference="eq:WEtheta"} yields $$\begin{equation}\label{eq:WE}
 \frac{\partial^2 y}{\partial t^2} = c^2 \frac{\partial^2 y}{\partial x^2}$$
where we have conveniently defined the wave velocity
$c = \sqrt{T_0 /\rho A}$,

This is an example of a *partial differential equation* (often
abbreviated as PDE). Such equations involve *partial* derivatives, and
relate the rate of change in time of a given quantity (here,
displacement), to its spatial variation. Equation
[\[eq:WE\]](#eq:WE){reference-type="eqref" reference="eq:WE"} is often
referred to as the *simple wave equation*, or just as *the wave
equation*. One first observation about the wave equation is that
acceleration is proportional to the second spatial derivative. Remember
that the second spatial derivative is a measure of the curvature of the
function: hence, if the curvature is positive at a given location along
the string, then force is pointing upward, and if the curvature is
negative, then the force is pointing downward. We are now going to
introduce some of the standard techniques for the solution of equations
such as [\[eq:WE\]](#eq:WE){reference-type="eqref" reference="eq:WE"}.

## Solutions to the 1D wave equation

*Separation of variables* is among the techniques often employed in the
analysis of PDEs. For that, we assume that $$\begin{equation}\label{eq:sepVar}
y(x,t) = X(x)T(t)\end{equation}$$
Although there is no specific reason why a solution
would have this form, we will assume that this is the case, and we will
check at the end whether the solution obtained via separation is
actually meaningful. Inserting
[\[eq:sepVar\]](#eq:sepVar){reference-type="eqref"
reference="eq:sepVar"} into [\[eq:WE\]](#eq:WE){reference-type="eqref"
reference="eq:WE"}, and dividing both sides by $XT$, one gets
$$\frac{T^{\prime\prime}}{T} = c^2 \frac{X^{\prime\prime}}{X}\end{equation}$$
Note
that, in the equation above, primes indicate derivatives with respect to
the function's own argument (i.e. $t$ for $T$, and $x$ for $X$.) Notice
also that the left-hand side depends only on $t$, and the right-hand
side depends only on $x$: hence, we may conveniently write

::: subnumcases
[\[eq:TXsyst\]]{#eq:TXsyst label="eq:TXsyst"} = s\^2
[\[eq:TXsyst1\]]{#eq:TXsyst1 label="eq:TXsyst1"}\
c\^2 = s\^2 [\[eq:TXsyst2\]]{#eq:TXsyst2 label="eq:TXsyst2"}
:::

where $s$ is a constant that does not depend on $t$ or $x$. The interest
of using this system is that the two equations are now in the form of
*ordinary differential equations* (ODE) of the kind already encountered
whilst studying the oscillator. We can then solve
[\[eq:TXsyst1\]](#eq:TXsyst1){reference-type="eqref"
reference="eq:TXsyst1"} using the exponential technique. This gives
$$T(t) = A_+ e^{st} + A_- e^{-st}, \quad A_+,A_- \in \mathbb{C}$$
Clearly, for this solution to be oscillating and not exponentially
growing or decaying one must set $s \in j\mathbb{R}$ and, in fact, one
may define $$s = j \omega, \quad \omega \in \mathbb{R}\end{equation}$$
We are now in
the position of interpreting
[\[eq:TXsyst1\]](#eq:TXsyst1){reference-type="eqref"
reference="eq:TXsyst1"}: this is just the equation of a simple harmonic
oscillator in time, with radian frequency given by $\omega$. The fact
that the temporal component oscillates at such frequency ensures that
the spatial component is also oscillating: this is encaspsulated in
[\[eq:TXsyst\]](#eq:TXsyst){reference-type="eqref"
reference="eq:TXsyst"}. This results in the *dispersion relation* for
the wave equation: $$\begin{equation}\label{eq:DispRelWECnt}
\omega = c \gamma,\end{equation}$$
where $\gamma$ is the spatial frequency, measured
in rad/m, and called the *wavenumber*. The solution to
[\[eq:TXsyst2\]](#eq:TXsyst2){reference-type="eqref"
reference="eq:TXsyst2"} is obtained analogously. One has
$$X(x) = B_+ e^{j\gamma x} + B_- e^{-j\gamma x}, \quad B_+,B_- \in \mathbb{C}$$
Thus, inserting these in
[\[eq:sepVar\]](#eq:sepVar){reference-type="eqref"
reference="eq:sepVar"} gives $$\begin{equation}\label{eq:genSolWE}
y(x,t) = \left( A_+ e^{j\omega t} + A_- e^{-j\omega t}\right)\left( B_+ e^{j\gamma x} + B_- e^{-j\gamma x}\right)$$
We remark two things: the first, is that the solution $y(x,t)$ is now
expressed as a combination of complex-valued functions. The second, is
that such combination depends on four arbitrary constants
$A_+,A_-,B_+,B_-$. We make the observation that these constant can be
adjusted so that the general solution
[\[eq:genSolWE\]](#eq:genSolWE){reference-type="eqref"
reference="eq:genSolWE"} satisfies the *initial conditions* in time,
plus the *boundary conditions* at the string's ends. We will come back
to these conditions later. For the moment, we will assume that the
solution [\[eq:genSolWE\]](#eq:genSolWE){reference-type="eqref"
reference="eq:genSolWE"} exists for all times, and that the string's
domain is the entire real axis (so that the boundaries are located at
$\pm \infty$): in this case, the four constants are completely
arbitrary, and we may use them to construct *real* solutions.

### Standing sinusoid as sum of travelling sinusoids

![Sum of travelling waves giving rise to a standing wave. The black wave
is travelling to the left, and the grey wave to the right. The sum of
the two is depicted as crosses (+), and is a standing wave. Here,
$\omega = 10\pi$, $\gamma = 4\pi$: the wavelength is
$\lambda = 2\pi/\gamma = 0.5$m, as visible from the
plots.](Figures/TravWavesSum.png){#fig:travStd width="\\linewidth"}

Under some arbitrary choice of the four constants in
[\[eq:genSolWE\]](#eq:genSolWE){reference-type="eqref"
reference="eq:genSolWE"}, one gets $$\begin{equation}\label{eq:SttoTr}
y(x,t) = \cos (\omega t )\cos (\gamma x )\end{equation}$$
You may think of this
function as a standing wave pattern that extends across the whole
$x$-axis. You may notice that, at the locations given by
$x = (2n+1)\pi /2 \gamma$, $n\in\mathbb Z$, the solution is always zero:
these are the *nodes* of the standing wave. On the other hand, at the
locations given by $x = n\pi/\gamma$, $n\in\mathbb Z$, the solution
reaches its maximum amplitude ($\pm 1$) at the times $t = m\pi/\omega$,
$m\in\mathbb Z$, see also Figure
[1.3](#fig:travStd){reference-type="ref" reference="fig:travStd"}. The
distance between two maxima is then $\lambda = 2\pi / \gamma$: this
distance is called the *wavelength*, and $\gamma$ is called the
$\emph{wavenumber}$ (i.e. it is the spatial radian frequency). Using
basic trigonometric identities, one may re-write
[\[eq:SttoTr\]](#eq:SttoTr){reference-type="eqref"
reference="eq:SttoTr"} as $$\begin{equation}\label{eq:SumTrav}
y(x,t) = \frac{1}{2}\left(\cos (\omega t + \gamma x) + \cos (\omega t - \gamma x ) \right)$$
Now, the term $\cos (\omega t + \gamma x)$ is a cosine function
*travelling to the left*; the term $\cos (\omega t - \gamma x )$ is a
cosine function *travelling to the right*, see also Figure
[1.3](#fig:travStd){reference-type="ref" reference="fig:travStd"}. The
waves are travelling with speed $\pm c=\pm\sqrt{T_0/ \rho A}$. Why are
these travelling waves? Just think of the points of constant phase:
these are $\omega t \pm \gamma x = \text{constant}$. Deriving the
expression above one obtains
$$v_\phi = \frac{dx}{dt} = \pm \frac{\omega}{\gamma} = \pm c\end{equation}$$
The
velocity of constant phase, indicated here as $v_\phi$, is called the
*phase velocity*. The absolute value of the phase velocity is given by
$|\omega/\gamma|$, and for the simple wave equation such ratio is always
equal to $c$.

The passage from [\[eq:SttoTr\]](#eq:SttoTr){reference-type="eqref"
reference="eq:SttoTr"} to
[\[eq:SumTrav\]](#eq:SumTrav){reference-type="eqref"
reference="eq:SumTrav"}, explained in terms of trigonometry, bears in
fact some consequences from a physical point of view: we showed that a
standing wave can be written as the sum of two travelling waves going
opposite directions.

### Travelling sinusoid as sum of standing sinusoids

![Sum of standing waves giving rise to a travelling wave. The black wave
is a standing wave of the form $\sin(\omega t)\sin(\gamma x)$, and the
grey wave is a standing wave of the form $\cos(\omega t)\cos(\gamma x)$.
The sum of the two is depicted as crosses (+), and is a travelling wave
to the right. Here, $\omega = 10\pi$,
$\gamma = 4\pi$.](Figures/StdWavesSum.png){#fig:stdTrav
width="\\linewidth"}

The converse is also true. Starting again from Equation
[\[eq:genSolWE\]](#eq:genSolWE){reference-type="eqref"
reference="eq:genSolWE"}, we may arrange the four arbitrary constants so
to obtain $$y(x,t) = \cos(\omega t - \gamma x)\end{equation}$$
This is of course the
expression of a cosine travelling to the right with phase velocity
$v_\phi = c$. Using again appropriate trigonometric identities, one has
$$y(x,t) = \cos(\omega t)\cos(\gamma x) + \sin(\omega t)\sin(\gamma x)$$
and, of course, this expression indicates that we can write the
travelling cosine as the sum of two standing waves, see also Figure
[1.4](#fig:stdTrav){reference-type="ref" reference="fig:stdTrav"}.

### General travelling wave solutions for $y(x,t)$ for the 1D wave equation

The simple wave equation allows solutions that have a more general
character than the sinusoidal solutions found so far. First, rewrite the
wave equation as $$\begin{equation}\label{eq:WEfac}
\left( \frac{\partial}{\partial t} - c  \frac{\partial}{\partial x}\right)\left( \frac{\partial}{\partial t} + c  \frac{\partial}{\partial x}\right) y(x,t) = 0$$
Second, we are now performing a change of variables. Hence, we are
mapping $(x,t) \in \mathbb{R}^2$ onto $(\eta,\zeta) \in \mathbb{R}^2$ so
that $$\eta = x + ct, \qquad  \zeta = x - ct.\end{equation}$$
We are then using the
chain rule to take partial derivatives in
[\[eq:WEfac\]](#eq:WEfac){reference-type="eqref" reference="eq:WEfac"}.
Thus

::: subequations
$$\begin{aligned}
\frac{\partial}{\partial t} &= \frac{\partial \eta}{\partial t}\frac{\partial }{\partial \eta} + \frac{\partial \zeta}{\partial t}\frac{\partial }{\partial \zeta} = c \frac{\partial }{\partial \eta} - c \frac{\partial }{\partial \zeta} \\
\frac{\partial}{\partial x} &= \frac{\partial \eta}{\partial x}\frac{\partial }{\partial \eta} + \frac{\partial \zeta}{\partial x}\frac{\partial }{\partial \zeta} = \frac{\partial }{\partial \eta} + \frac{\partial }{\partial \zeta} \end{aligned}$$
:::

Using these into [\[eq:WEfac\]](#eq:WEfac){reference-type="eqref"
reference="eq:WEfac"} gives the following PDE
$$\frac{\partial^2y}{\partial \eta \partial \zeta}  = 0\end{equation}$$
Hence, direct
integration gives $$y(\eta,\zeta) = f(\eta) + g(\zeta)\end{equation}$$
or,
$$\begin{equation}\label{eq:Dalem}
y(x,t) = f(x+ct) + g(x-ct)\end{equation}$$
In view of the discussion in the previous
sections, you will recognise that $f,g$ have the form of a leftward and
rightward travelling wave, respectively. Notice that this result is more
general than the sinusoidal solutions of the previous section: here,
$f,g$ are *any* function of $x+ct$ and $x-ct$. Examples of functions
that are solution to the wave equation: $$\begin{equation}\label{eq:WEsolsEx}
y_1(x,t) = e^{(x-ct)}, \,\,\, y_2(x,t) = \log \left( (x+ct)^2 \right), \,\,\, y_3(x,t) =  ...$$
Beware that nonlinear combinations (i.e. products) of $f,g$ are not, in
general, solutions. Examples of functions that *are not* solutions to
the wave equation:
$$p_1(x,t) = (x-ct)e^{x+ct}, \,\,\, p_2(x,t) = \log \left( (x+ct)^2 \right)\left(\cos(x-ct)\right)^2, \,\,\, p_3(x,t) = ...$$
This property of allowing solutions such as $f(x+ct)$ and $g(x-ct)$ is
not valid for other one dimensional PDEs in acoustics: the simple wave
equation is special in this sense. This equation describes the motion of
ideal strings, and also the vibration of air columns with cilindrical
geometry. These forms for $f,g$ were first derived by D'Alembert in
1747, and are known as *D'Alembert solutions*. They are at the core of
early physical modelling techniques based around *digital waveguides*.

Important note. Although function such as those in
[\[eq:WEsolsEx\]](#eq:WEsolsEx){reference-type="eqref"
reference="eq:WEsolsEx"} solve the wave equation mathematically, they
might not be \"energy-compatible\". In practice, not only does a
function need to solve [\[eq:WE\]](#eq:WE){reference-type="eqref"
reference="eq:WE"}, but it also must comply with some requirements at
the string's ends, the boundary conditions. This aspect will become
clearer as we progress through the course.

## Energy for the 1D wave equation

In order to derive energy expressions for the string, we might try to
relate the work done by the tension forces to stretch an infinitesimal
element of string, as that in Figure
[1.2](#fig:infinitesimal_string){reference-type="ref"
reference="fig:infinitesimal_string"}. The corresponding change in
potential energy is given by the stretching in the deformed
configuration, minus the stretching in the undeformed configuration,
times the tension $T_0$. Hence
$$dE_p =  T_0 (ds-dx) =  T_0  \left(\sqrt{1 + \left( \frac{\partial y}{\partial x}\right)^2} - 1 \right)dx$$
To first order, one has
$$dE_p \approx \frac{T_0}{2}\left(\frac{\partial y}{\partial x}\right)^2 dx$$
From Newtonian physics of conservative systems, we know that sum of the
changes of potential and kinetic energies is equal to zero, hence
$$dE_k + dE_p = 0\end{equation}$$
where
$dE_k = \frac{\rho A}{2} \left( \frac{\partial y}{\partial t} \right)^2 dx$.
Integrating this last expression between $x_l$ and $x_r$ ($x_l<x_r$)
gives $$\begin{equation}\label{eq:EnBal}
E_k + E_p = \int_{x_l}^{x_r} \frac{\rho A}{2} \left(\frac{\partial y}{\partial t}\right)^2 dx + \int_{x_l}^{x_r} \frac{T_0}{2} \left(\frac{\partial y}{\partial x}\right)^2 dx = E$$
where $E$ is the total energy (a constant). Notice that the total energy
is a function of time only.

### Energy balance for an infinitely long string

We can try to check whether the energy balance equation
[\[eq:EnBal\]](#eq:EnBal){reference-type="eqref" reference="eq:EnBal"}
actually provides a consistent definition for the total energy. To that
extent, take a time derivative of
[\[eq:EnBal\]](#eq:EnBal){reference-type="eqref" reference="eq:EnBal"}.
This is
$$\frac{d}{dt}\left( E_k + E_p \right) = \int_{x_l}^{x_r} \rho A \left(\frac{\partial^2 y}{\partial t^2}\right)\left(\frac{\partial y}{\partial t}\right) dx + \int_{x_l}^{x_r} T_0 \left(\frac{\partial^2 y}{\partial t\partial x}\right)\left(\frac{\partial y}{\partial x}\right) dx = 0$$
For the moment, we are only concerned with infinitely longs strings, so
that $x_l = -\infty$, $x_r = \infty$. Using integration by parts to
calculate the second integral, one can rewrite the energy balance as
$$\int_{-\infty}^{\infty} \rho A \left(\frac{\partial^2 y}{\partial t^2}\right)\left(\frac{\partial y}{\partial t}\right) dx - \int_{-\infty}^{\infty} T_0 \left(\frac{\partial^2 y}{\partial x^2}\right)\left(\frac{\partial y}{\partial t}\right) dx = - T_0 \left(\frac{\partial y}{\partial t}\right)\left(\frac{\partial y}{\partial x}\right)\bigg|_{-\infty}^{\infty}$$
We make the assumption that the right-hand side vanishes. (In fact, this
puts some limitations on the form of the solutions to the wave equation.
What can you say about the functions in
[\[eq:WEsolsEx\]](#eq:WEsolsEx){reference-type="eqref"
reference="eq:WEsolsEx"}?)

Collecting the common factor in the integrals, we finally get
$$\int_{-\infty}^{\infty} \left[ \rho A \left(\frac{\partial^2 y}{\partial t^2}\right) - T_0 \left(\frac{\partial^2 y}{\partial x^2}\right) \right]\left(\frac{\partial y}{\partial t}\right) dx  = 0$$
Comparing this last expression with
[\[eq:WE\]](#eq:WE){reference-type="eqref" reference="eq:WE"}, we see
that the integrand is in fact identically zero. Thus, we observe that
the energy expressions in
[\[eq:EnBal\]](#eq:EnBal){reference-type="eqref" reference="eq:EnBal"}
are consistent with the wave equation, provided that the solutions
behave as required at the boundary of the domain.

## Initial Conditions

We are now going to try to prove that, regardless of the approach we
wish to adopt to describe the solutions, we are going to end up with the
same solution. In order to see this, we are now going to introduce
*initial conditions* on the string, such that we may assume that the
string is at rest for $t<0$, and that some prescribed initial conditions
are given at $t=0$. For simplicity, only displacement initial conditions
will be studied here, so that $$\begin{equation}\label{eq:ICs}
y(x,0) = y_0(x), \qquad \frac{\partial y(x,0)}{\partial t} = v_0(x).$$
For simplicity, only displacement initial conditions will be studied in
the following, for which $v_0(x) = 0$. We are now going to solve the
wave equation, with these initial conditions, using both the D'Alembert
solution, and the solution in terms of sinusoids (in fact, a Fourier
method).

### D'Alembert solution with prescribed initial conditions

Using [\[eq:Dalem\]](#eq:Dalem){reference-type="eqref"
reference="eq:Dalem"} into [\[eq:ICs\]](#eq:ICs){reference-type="eqref"
reference="eq:ICs"} gives
$$f(x)+g(x) = y_0(x), \qquad f^\prime(x) - g^\prime(x) = 0\end{equation}$$
Integrating
the second condition, one gets $f(x) = g(x)$, and using this into the
first condition one gets $$f(x)=g(x)=\frac{y_0(x)}{2}\end{equation}$$
so that
$$\begin{equation}\label{eq:Tr}
y(x,t) = \frac{1}{2}\left( y_0(x-ct) + y_0(x+ct)\right)\end{equation}$$
The meaning of
this last equation is that any initial displacement will split into two
travelling wavefronts, of amplitude equal to one half that of the
initial displacement. (If the initial velocity is not identically zero,
the solution is a little more involved, but it is again separable into
two travelling wavefronts.)

### Fourier solution with prescribed initial conditions

If instead we wish to use the sinusoidal representation for the
solution, note that the most general form is given by the Fourier
integral $$\begin{equation}\label{eq:FTy}
y(x,t) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{\infty} \hat y(\gamma,t) e^{j\gamma x}d\gamma$$
where $\hat y(\gamma,t)$ is the Fourier transform of $y(x,t)$, defined
as $$\begin{equation}\label{eq:yhat}
\hat y(\gamma,t) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{\infty}  y(x,t) e^{-j\gamma x}d x$$
The meaning of [\[eq:FTy\]](#eq:FTy){reference-type="eqref"
reference="eq:FTy"} is that the solution can be represented as a
superposition of sine / cosine functions of $\gamma x$, times an
appropriate complex amplitude $\hat y(\gamma,t)$. The complex amplitude
can be calculated by taking a Fourier transform of the kind
[\[eq:yhat\]](#eq:yhat){reference-type="eqref" reference="eq:yhat"} of
the wave equation [\[eq:WE\]](#eq:WE){reference-type="eqref"
reference="eq:WE"}. This gives
$$\frac{\partial^2 \hat y(\gamma,t)}{\partial t^2} = - c^2\gamma^2 \hat y(\gamma,t),$$
which is the usual harmonic oscillator equation solved by
$$\begin{equation}\label{eq:yhh}
\hat y(\gamma,t) = A_+(\gamma) e^{jc\gamma t} + A_-(\gamma)-e^{-jc\gamma t}$$
We can now impose the initial conditions
[\[eq:ICs\]](#eq:ICs){reference-type="eqref" reference="eq:ICs"}
(transformed in the Fourier domain)
$$\hat y(\gamma,0) = \hat y_0(\gamma), \qquad \frac{\partial \hat y (\gamma,0)}{\partial t} = 0$$
which give $$A_+(\gamma) = A_-(\gamma) = \frac{\hat y_0(\gamma)}{2}$$
Using these into [\[eq:yhh\]](#eq:yhh){reference-type="eqref"
reference="eq:yhh"} and that into
[\[eq:FTy\]](#eq:FTy){reference-type="eqref" reference="eq:FTy"} gives
the solution
$$y(x,t) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{\infty}\frac{\hat y_0(\gamma)}{2} \left( e^{-j c\gamma t} + e^{j c\gamma t} \right) e^{j\gamma x} d\gamma$$
which may be written as $$\begin{equation}\label{eq:FTtoTr}
y(x,t) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{\infty}\frac{\hat y_0(\gamma)}{2}  e^{j\gamma (x - c  t)}d\gamma +  \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{\infty}\frac{\hat y_0(\gamma)}{2}  e^{j\gamma (x + c  t)}d\gamma  = \frac{1}{2}\left( y_0(x-ct) + y_0(x+ct)\right)$$
where the last equality is simply an application of the definition of
the Fourier transform as per
[\[eq:FTy\]](#eq:FTy){reference-type="eqref" reference="eq:FTy"}.
Comparing [\[eq:FTtoTr\]](#eq:FTtoTr){reference-type="eqref"
reference="eq:FTtoTr"} with [\[eq:Tr\]](#eq:Tr){reference-type="eqref"
reference="eq:Tr"}, we see that the result is the same. Note as well
that, in [\[eq:FTtoTr\]](#eq:FTtoTr){reference-type="eqref"
reference="eq:FTtoTr"}, we effectively write the solution as an infinite
combination (an integral) of travelling sinusoidal functions of
amplitude given by the intial conditions.

The net result of this discussion is that it should not matter how you
decide to write out your solution: you should always end up with the
same thing once initial conditions are imposed! Note that, so far, we
have not yet talked of boundary conditions, as we have assumed that the
string is infinite in length. The rest of today's lecture is to
understand how the introduction of bounded domains imposes further
restrictions on the nature of the solutions.

## Energy of a semi-infinite string. Boundary conditions. {#sec:WaveEquationsBCs}

The nature of boundary conditions may be understood from energy
arguments. Recall that the energy of the flexible string is
$$E = E_k + E_p = \int_{-\infty}^{\infty} \frac{\rho A}{2}\left( \frac{\partial y}{\partial t} \right)^2 dx + \int_{-\infty}^{\infty} \frac{T_0}{2}\left( \frac{\partial y}{\partial x} \right)^2 dx$$
The energy balance is expressed as $dE/dt = 0$. However, these results
are only valid when the string is infinite in length. Suppose now that
the string is only defined for $x\geq 0$. In practice we have introduced
a boundary point at $x=0.$ The energy is given by $$\begin{equation}\label{eq:EnHalf}
E = E_k + E_p = \int_{0}^{\infty} \frac{\rho A}{2}\left( \frac{\partial y}{\partial t} \right)^2 dx + \int_{0}^{\infty} \frac{T_0}{2}\left( \frac{\partial y}{\partial x} \right)^2 dx$$
However, the energy balance is now $$\begin{equation}\label{eq:EnBal}
\frac{dE}{dt} = - T_0 \frac{\partial y}{\partial x }\Bigg|_{x=0}\frac{\partial y}{\partial t }\Bigg|_{x=0}$$
As you can see, the boundary term can now influence the energy balance
in such a way as to produce or disspate power! (Note that the energy
balance is expressed as the product of a force times a velocity at the
boundary point.) One way of ensuring that $dE/dt = 0$ (as per the string
of infinite length), is then to require that either the force, or the
velocity, is zero at the boundary. Thus, one may set
$$y\big|_{x=0} = 0 \,\,\, \forall t, \,\,\, \text{or } \,\,\, \frac{\partial y}{\partial x }\Bigg|_{x=0} = 0 \,\,\, \forall t$$
The first condition corresponds to zero displacement / velocity, and is
therefore called *fixed*, or of *Dirichlet* type. The second condition
corresponds to vanishing force at the boundary, and is called *free*, as
in free of load, or of *Neumann* type.

The question is, it is not immediately obvious to see why the energy
balance should have the form
[\[eq:EnBal\]](#eq:EnBal){reference-type="eqref" reference="eq:EnBal"}.
The most straightforward way of deriving such energy balance is to start
from the equation of motion, multiplied by the velocity, and to
integrate over the string's length. Hence
$$\rho A \int_0^{\infty} \frac{\partial^2 y}{\partial t^2} \frac{\partial y}{\partial t} dx = T_0 \int_0^{\infty} \frac{\partial^2 y}{\partial x^2} \frac{\partial y}{\partial t} dx$$
The left-hand side can be expressed as
$$\int_0^{\infty}\frac{\partial^2 y}{\partial t^2} \frac{\partial y}{\partial t} dx = \frac{1}{2}\frac{d}{dt}\int_0^{\infty} \left( \frac{\partial y}{\partial t}\right)^2 dx$$
The right-hand side can be integrated by parts, and written as
$$\int_0^{\infty} \frac{\partial^2 y}{\partial x^2} \frac{\partial y}{\partial t} dx = - T_0 \frac{\partial y}{\partial x }\Bigg|_{x=0}\frac{\partial y}{\partial t }\Bigg|_{x=0} - \frac{1}{2}\frac{d}{dt}  \int_0^{\infty}\left( \frac{\partial y}{\partial x}\right)^2 dx$$
Thus, we recover the energy balance
[\[eq:EnBal\]](#eq:EnBal){reference-type="eqref" reference="eq:EnBal"}
with energy components given by
[\[eq:EnHalf\]](#eq:EnHalf){reference-type="eqref"
reference="eq:EnHalf"}.

### D'Alembert solution with prescribed boundary condition

We now wish to describe the motion on the semi-infinite string using the
D'Alembert solution. Consider first a boundary condition of *fixed*
type. Using this condition one gets
$$y(0,t) = f(-ct) + g(ct) = 0, \,\,\, \forall t\end{equation}$$
This corresponds to
setting $f(-ct) = -g(ct)$. An interpretation of this result is given in
Figure [1.5](#fig:DalFixed){reference-type="ref"
reference="fig:DalFixed"}: essentially, the solution can be thought as
the sum of the actual solution, defined for $x\geq 0$, plus a "ghost"
solution defined over $x<0$. The ghost solution is just the flipped
image of the right-travelling wave $f(x-ct)$. When the actual left-going
solution enters the negative $x$ axis, the corresponding flipped ghost
solution enters the positive $x$ axis. The net result is that the actual
solution gets reflected at the boundary $x=0$ with a sign change.

![D'Alembert solution with fixed boundary at $x=0$. Thick black line is
the actual solution for $x\geq 0$. The line with crosses (+) is the
ghost solution defined for
$x<0$.](Figures/FixedBCTravel.png){#fig:DalFixed width="\\linewidth"}

Turning now the attention to the *free* boundary condition, the
mathematical condition here is
$$\frac{\partial y}{\partial x}(0,t)=-f^\prime(-ct)+g^\prime(ct) = 0 \,\,\, \forall t$$
which can be integrated to give $f(-ct)=g(ct)$

![D'Alembert solution with free boundary at $x=0$. Thick black line is
the actual solution for $x\geq 0$. The line with crosses (+) is the
ghost solution defined for
$x<0$.](Figures/FreeBCTravel.png){#fig:DalFree width="\\linewidth"}

This is sketched in Figure [1.6](#fig:DalFree){reference-type="ref"
reference="fig:DalFree"}. Here the ghost solution has the same sign as
the right-travelling wave, so that at the boundary $x=0$ the slope of
the actual solution is zero (there is a stationary point there.)

### Fourier solution with prescribed boundary condition

Suppose now that we wish to describe the solution in terms of Fourier
components (i.e. sine and cosine functions.) In analogy to the case
developed for the initial conditions, we are now taking a Fourier
transform in the $\omega$ domain, so that the solution may be
represented as $$\begin{equation}\label{eq:ytt}
y(x,t) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{\infty} \hat y(x,\omega) e^{j\omega t} d\omega$$
where the Fourier transform $\hat y(x,\omega)$ is defined as
$$\begin{equation}\label{eq:yhattOm}
\hat y(x,\omega) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{\infty}  y(x,t) e^{-j\omega t} d t$$
(Note that, with a slight abuse of notation, we have identified with the
symbol $\hat y$ both the Fourier transform in the $\gamma$ domain, in
[\[eq:yhat\]](#eq:yhat){reference-type="eqref" reference="eq:yhat"}, and
the Fourier transform in the $\omega$ domain in
[\[eq:yhattOm\]](#eq:yhattOm){reference-type="eqref"
reference="eq:yhattOm"}. These are of course distinct functions.) Using
the transformed wave equation, we get
$$-\omega^2 \hat y(x,\omega) = c^2 \frac{\partial^2 \hat y(x,\omega)}{\partial x^2}$$
which is solved, as per usual, by $$\begin{equation}\label{eq:bcss}
\hat y(x,\omega) = B_+(\omega) e^{j\frac{\omega}{c}x} + B_-(\omega) e^{-j\frac{\omega}{c}x}$$
Using the fixed boundary condition at $x=0$ gives $B_+ = -B_-$ and thus
$$\begin{equation}\label{eq:FixedFou1}
\hat y(x,\omega) = 2jB_+(\omega) \sin \frac{\omega x}{c}\end{equation}$$
The solution
is then given by using this expression in
[\[eq:ytt\]](#eq:ytt){reference-type="eqref" reference="eq:ytt"}
$$y(x,t)=\int_{-\infty}^{\infty} C_+(\omega) \sin\frac{\omega x}{c} e^{j\omega t} d\omega$$
where we defined $C_+ = 2 j B_+/\sqrt{2\pi}$. Notice that, indeed, we
could have taken a Fourier transform in
[\[eq:ytt\]](#eq:ytt){reference-type="eqref" reference="eq:ytt"} using
$e^{-j \omega t}$ instead of $e^{j \omega t}$, so that effectively a
more general solution is given by $$\begin{equation}\label{eq:FixedFou}
y(x,t)=\int_{-\infty}^{\infty} C_+(\omega) \sin\frac{\omega x}{c} e^{j\omega t} d\omega +  \int_{-\infty}^{\infty} C_-(\omega) \sin\frac{\omega x}{c} e^{-j\omega t} d\omega$$
where $C_+,C_-$ may be determined from the two initial conditions.

Using instead the free boundary condition in
[\[eq:bcss\]](#eq:bcss){reference-type="eqref" reference="eq:bcss"}, we
get $B_+ = B_-$ and $$\begin{equation}\label{eq:FreeFou}
\hat y(x,\omega) = 2 B_+(\omega) \cos \frac{\omega x}{c}\end{equation}$$
and thus, the
general solution is $$\begin{equation}\label{eq:FreeFou}
y(x,t)= \int_{-\infty}^{\infty} D_+(\omega) \cos\frac{\omega x}{c} e^{j\omega t} d\omega +  \int_{-\infty}^{\infty} D_-(\omega) \cos\frac{\omega x}{c} e^{-j\omega t} d\omega$$
where $D_+,D_-$ are again determined by the initial conditions. These
meaning of [\[eq:FixedFou\]](#eq:FixedFou){reference-type="eqref"
reference="eq:FixedFou"} is that, no matter the time evolution of the
solution, the spatial distribution is such that the solution can be
written as a linear combination of sine functions (that are always zero
at the origin.) On the other hand, a free condition is such that the
solution can be expressed as a linear combination of cosine functions,
that have zero slope at the origin, as in
[\[eq:FreeFou\]](#eq:FreeFou){reference-type="eqref"
reference="eq:FreeFou"}.

## Strings of finite length

The analysis carried out so far is very useful, as we are now ready to
analyse the case of a string of finite length. A good place to start,
once more, is the energy balance for the string, now defined for
$0\leq x \leq L$. The energy balance is $$\begin{equation}\label{eq:EnBalFull}
\frac{dE}{dt} = T_0 \frac{\partial y}{\partial x }\Bigg|_{x=L}\frac{\partial y}{\partial t }\Bigg|_{x=L} - T_0 \frac{\partial y}{\partial x }\Bigg|_{x=0}\frac{\partial y}{\partial t }\Bigg|_{x=0}$$
The origin of this expression should now be clear to you. The idea is to
perform the same mathematical steps that lead to
[\[eq:EnBal\]](#eq:EnBal){reference-type="eqref" reference="eq:EnBal"},
but where now the limits of integration are set between $0$ and $L$.
Hence $$\begin{equation}\label{eq:EnFinite}
E = E_k + E_p = \int_{0}^{L} \frac{\rho A}{2}\left( \frac{\partial y}{\partial t} \right)^2 dx + \int_{0}^{L} \frac{T_0}{2}\left( \frac{\partial y}{\partial x} \right)^2 dx$$
In order to study the motion of strings of finite length, we may again
adopt either the D'Alembert or the Fourier solution. Notice that the
boundary conditions are exactly symmetric at the two string's ends. So,
the same interpretation of a fixed or a free boundary can be given at
$x=L$. This idea of deriving boundary conditions from energy analysis is
indeed very powerful, and will be used a lot during this course.

### Boundedness of the solutions {#sec:BoundWECnt}

Given the energy balance
[\[eq:EnFinite\]](#eq:EnFinite){reference-type="eqref"
reference="eq:EnFinite"}, it is immediate to see that
$$\int_0^L \left( \frac{\partial y}{\partial t} \right)^2 \, \mathrm{d}x \leq \frac{2E}{\rho A}, \quad \int_0^L \left( \frac{\partial y}{\partial x} \right)^2 \, \mathrm{d}x \leq \frac{2E}{T_0},$$
expressing the bounds on the norms of the velocity and the gradient.
These are generalisations to the continuous case of the bounds found for
the oscillator in Section
[\[sec:BoundsSHO\]](#sec:BoundsSHO){reference-type="ref"
reference="sec:BoundsSHO"}.

### D'Alembert solution for strings of finite length

Given this interpretation for the boundary conditions, the D'Alembert
solution should be pretty obvious to you now: at either end of the
string, the right/left travelling wave front gets reflected exactly,
with a sign change for a fixed boundary, and without a sign change for a
free boundary. The left-going and right-going wavefronts at the time
$t=0$ may be determined using [\[eq:Tr\]](#eq:Tr){reference-type="eqref"
reference="eq:Tr"}. There is nothing complicated here. This solution is
indeed powerful, both because of its simplicity, and because of the ease
of implementation in the digital domain via a digital waveguide
appraoch. However, it has the drawback of not being generalisable to
other systems (outside of cylindrical air columns).

### Fourier solution for strings of finite length

On the other hand, the Fourier approach has the advantdge of being, to
some extent, universal: you may indeed apply this method to any linear,
time-invariant system, in one or more dimenensions. The principle of
application of this method has already been illustrated in the previous
examples of strings of infinite length. Suppose now that the string is
fixed at both ends, so that $$y(0,t) = y(L,t) = 0\end{equation}$$
We have already
derived the expression for the solution with a fixed boundary at $x=0$.
This has led to [\[eq:FixedFou1\]](#eq:FixedFou1){reference-type="eqref"
reference="eq:FixedFou1"}. Now, from
[\[eq:FixedFou1\]](#eq:FixedFou1){reference-type="eqref"
reference="eq:FixedFou1"}, given that the string is fixed at $x=L$, we
should also have $$\begin{equation}\label{eq:omeM}
y(L,t) = 0 \implies \sin \frac{\omega L}{c} = 0 \implies \omega = \omega_m = \frac{\pi m c}{L}, \,\,\, m \in \mathbb{N}$$
This is a very important result: what
[\[eq:omeM\]](#eq:omeM){reference-type="eqref" reference="eq:omeM"} is
saying is that only certain vibrational frequencies are now allowed,
those that multiples of the fundamental frequency
$\omega_1 = \pi c / L$. In other words, the effect of finite string
length is to select only a discrete set of frequencies $\{ \omega_m \}$,
whilst in strings of infinite length sinusoidal components of any
frequency may be found. This has a consequence on how we express the
general solution [\[eq:FixedFou\]](#eq:FixedFou){reference-type="eqref"
reference="eq:FixedFou"}. Instead of integrals, now we have
$$\begin{equation}\label{eq:modal}
y(x,t) = \sum_{m=1}^{\infty } \left( C_+^{(m)} \sin \frac{\omega_m x}{c} e^{j \omega_m t} +  C_-^{(m)} \sin \frac{\omega_m x}{c} e^{-j \omega_m t} \right)$$
where $C_+^{(m)},C_-^{(m)}$ are determined from the initial conditions.
We call the expansion [\[eq:modal\]](#eq:modal){reference-type="eqref"
reference="eq:modal"} a *modal* expansion. The modes of the system are
defined as a set comprising both the *modal frequencies*
$\{ \omega_m \}$ and the *modal shapes* $\{\sin \frac{\omega_m x}{c}\}$.

The form of [\[eq:modal\]](#eq:modal){reference-type="eqref"
reference="eq:modal"} should be familiar to you: consider one point
along the string, say $x_p$. Then $$\begin{equation}\label{eq:Modes}
y(x_p,t) = \sum_{m=1}^{\infty } \left( \bar C_+^{(m,x_p)} e^{j \omega_m t} +  \bar C_-^{(m,x_p)}  e^{-j \omega_m t} \right)$$
with $\bar C_+^{(m,x_p)} = C_+^{(m)} \sin \frac{\omega_m x_p}{c}$,
$\bar C_-^{(m,x_p)} = C_-^{(m)} \sin \frac{\omega_m x_p}{c}$. But
[\[eq:Modes\]](#eq:Modes){reference-type="eqref" reference="eq:Modes"}
is just the sum of harmonic oscillators, whose amplitude is proportional
to the value of the modal shape at the selected point. Hence, the
vibration at each point alont the string can be thought of as resulting
from the sum of harmonic oscillators, of increasing frequency. The
oscillators are "transparent" to each other, meaning that the vibration
of each oscillator is not influenced by what any other oscillator is
doing. A picture of the modal shapes is given in Figure
[1.7](#fig:fixfix){reference-type="ref" reference="fig:fixfix"}.

![Modes of a fixed-fixed string. Here, the length of the string is
$L =1$ m, and the modes have been normalised so to have a maximum
amplitude of 1 mm.](Figures/FixedFixedModes.png){#fig:fixfix
width="\\linewidth"}

If we decide to impose different boundary conditions at either end
point, we end up with different expressions for the modal shapes.
Consider for instance the case of a *free-free* string. We already know
that a free end at $x=0$ on a semi-infite string results in a general
solution of the kind
[\[eq:FreeFou\]](#eq:FreeFou){reference-type="eqref"
reference="eq:FreeFou"}. We now bring in a second boundary point at
$x=L$, with the requirement that $\partial y / \partial x (L,t) = 0$.
This condition gives
$$\frac{\partial}{\partial x}\cos\frac{\omega x}{c}\Bigg|_{x=L} = - \frac{\omega}{c} \sin \frac{\omega L}{c} = 0 \implies \omega = \omega_m = \frac{\pi m c}{L}, \,\,\, m \in \mathbb{N}_0$$
These are the same frequencies that we obtained for the fixed-fixed
case, with the small difference that here the solution for $m=0$ (zero
frequency) gives rise to a modal shape that is *not* identically zero:
this special mode is a free-body mode, in which the string does not
oscillate but "floats" is space. (This is physically meaningful, as the
string is not attached at the ends and is therefore able to float. In
reality, though, this kind of boundary condition cannot be realised, as
the tension in the string arises from being attached at both boundaries.
If either one of the boundaries is not fixed, then no tension is
applied.) A picture of the modal shapes is given in Figure
[1.8](#fig:freefree){reference-type="ref" reference="fig:freefree"}.

Another case of interest is the fixed-free case. Again, we should start
from the expression of the semi-infinite string with a fixed boundary at
$x=0$, [\[eq:FixedFou\]](#eq:FixedFou){reference-type="eqref"
reference="eq:FixedFou"}. Now, the second condition to impose is
$\partial y / \partial x (L,0) = 0$, which is
$$\frac{\partial}{\partial x}\sin\frac{\omega x}{c}\Bigg|_{x=L} =  \frac{\omega}{c} \cos \frac{\omega L}{c} = 0 \implies \omega = \omega_m = \frac{\pi (2m-1) c}{2L}, \,\,\, m \in \mathbb{N}$$
The frequencies are now odd multiples of the fundamental frequency
$\omega_1 = \pi c / 2 L$, which is one octave lower than the fundamental
frequencies of the fixed-fixed and free-free cases. A picture of the
modal shapes is given in Figure [1.9](#fig:fixfree){reference-type="ref"
reference="fig:fixfree"}.

![Modes of a free-free string. Here, the length of the string is $L =1$
m, and the modes have been normalised so to have a maximum amplitude of
1 mm. Note that the first mode is a rigid body mode at zero
frequency.](Figures/FreeFreeModes.png){#fig:freefree
width="\\linewidth"}

The last case of interest, the free-fixed case, is left to you as an
exercise. The interpretation of the modal expansion is that, for *any*
initial conditions, the final solution can be expressed as a combination
of standing waves, of shape and frequency determined by the boundary
conditions. Of course, as we have discussed at some length previously,
the net result of such combination of standing modes need not be a
standing wave: in fact, it will almost always be a travelling wavefront.

![Modes of a fixed-free string. Here, the length of the string is $L =1$
m, and the modes have been normalised so to have a maximum amplitude of
1 mm.](Figures/FixedFreeModes.png){#fig:fixfree width="\\linewidth"}

### Boundary condition with nonzero finite impedance

![Mass-spring-dashpot boundary condition. The string exerts a force on
the mass, $f(t) = -T_0 \partial y / \partial x (L,t)$. The equation of
motion of the mass is therefore $d^2 z / dt^2 = -kz -R dz/dt + f(t)$.
Here, $z(t) = y(L,t)$.](Figures/boundary.pdf){#fig:impBc
width="\\linewidth"}

The boundary conditions given above are a useful idealisation. The study
conducted so far has allowed us to understand the general physical
properties of ideal strings. Of course, real musical strings do not
behave in this idealised manner. More realistic string models will be
taken into account later in the course, when we include stiffness
effects (i.e. when we consider the effects of finite string thickness.)
Here, we want to analyse the effect of introducing more realistic
boundary conditions at one of the string's end. Let us consider a
fixed-free string and suppose that, at the right boundary, we attach
some kind of bridge. A useful model for a bridge is that of a
mass-spring system: in practice, the bridge is characterised by a mass
$m$, a stiffness $k$, and a loss parameter $R$. This kind of boundary
condition is sketched in Figure [1.10](#fig:impBc){reference-type="ref"
reference="fig:impBc"}.

We want to set up the equation of motion of the mass-spring system. Of
course, this is in the form of an oscillator equation, i.e.
$$\begin{equation}\label{eq:ImpBSc}
m \frac{d^2 z}{dt^2} = - k z - R \frac{d z}{dt} + f(t)\end{equation}$$
Here, the
displacement $z(t)$ of the bridge is equal to the displacement of the
string's end, so that $z(t) = y(L,t)$ (the bridge is attached to the
string at all times!) The external force $f(t)$ is the force exterted by
the string's end onto the mass. What is this force? We know this
expression from the energy balance
[\[eq:EnBalFull\]](#eq:EnBalFull){reference-type="eqref"
reference="eq:EnBalFull"}: $$\begin{equation}\label{eq:ft}
f(t) = - T_0 \frac{\partial y}{\partial x}(L,t)$$

Notice that the minus sign is here due to the fact that this is the
force of that the string exerts *onto* the mass. In practice, the
boundary displacement $z(t)$ now satisfies the differential equation
[\[eq:ImpBSc\]](#eq:ImpBSc){reference-type="eqref"
reference="eq:ImpBSc"}, and is not simply prescribed. The question is
whether the string plus mass-spring system satisfies an energy balance.
This is, if course, the case. Considering that the string is fixed at
$x=0$, the boundary term at $x=0$ vanishes in
[\[eq:EnBalFull\]](#eq:EnBalFull){reference-type="eqref"
reference="eq:EnBalFull"}. All we need to do is to check what happens at
$x=L$. To that extent, take a product of
[\[eq:ImpBSc\]](#eq:ImpBSc){reference-type="eqref"
reference="eq:ImpBSc"} with $dz/dt$. This gives
$$m \frac{d^2 z}{dt^2}\frac{d z}{dt} = - k z\frac{d z}{dt} - R \left(\frac{d z}{dt}\right)^2 + f(t)\frac{d z}{dt}$$
Using the usual energy identities, we get
$$\frac{d}{dt}\left( \underbrace{\frac{m}{2} \left(\frac{d z}{dt}\right)^2}_{e_k} + \underbrace{\frac{k z^2}{2}}_{e_p} \right) = - R \left(\frac{d z}{dt}\right)^2 + f(t)\frac{d z}{dt}$$
We are now taking the sum of
[\[eq:EnBalImp\]](#eq:EnBalImp){reference-type="eqref"
reference="eq:EnBalImp"} with
[\[eq:EnBalFull\]](#eq:EnBalFull){reference-type="eqref"
reference="eq:EnBalFull"}, to get
$$\frac{dE}{dt} = \frac{d}{dt}\left(E_k + E_p + e_k + e_p \right) = - R \left(\frac{d z}{dt}\right)^2 + f(t)\frac{d z}{dt} + T_0 \frac{\partial y}{\partial x}\Bigg|_{x=L} \frac{\partial y}{\partial t} \Bigg|_{x=L}$$
Of course, remembering the definitions of $f(t)$ in
[\[eq:ft\]](#eq:ft){reference-type="eqref" reference="eq:ft"} and of
$z(t)$, one has
$$f(t)\frac{d z}{dt} + T_0 \frac{\partial y}{\partial x}\Bigg|_{x=L} \frac{\partial y}{\partial t} \Bigg|_{x=L} = 0$$
and thus, the energy balance for this kind of boundary is
$$\begin{equation}\label{eq:EnBalImp}
\frac{dE}{dt} = - R \left(\frac{d z}{dt}\right)^2 \leq 0\end{equation}$$
This is a
comforting (and expected!) result. The energy in the system is
decreasing because of the loss at the boundary, and is perfectly
conserved if $R=0$ (conservative boundary.) The added mass and spring
are energy-storing devices, and therefore the total energy now includes
the corresponding kinetic and potential energies. A boundary condition
such as [\[eq:ImpBSc\]](#eq:ImpBSc){reference-type="eqref"
reference="eq:ImpBSc"} has a nonzero, finite impedance, as opposed to
the fixed condition (infinite impedance), or the free condition (zero
impedance). We can now analyse the effect of various cases.

![Graphical solution of trascendental frequency equations. Left:
[\[eq:bc1\]](#eq:bc1){reference-type="eqref" reference="eq:bc1"}. Right:
[\[eq:bc2\]](#eq:bc2){reference-type="eqref"
reference="eq:bc2"}](Figures/tanSol.png){#fig:trasc width="\\linewidth"}

Case a): $m=0,R=0$. Here, we are assuming that the mass and loss
parameter of the bridge are very small compared to the stiffness. The
boundary condition [\[eq:ImpBSc\]](#eq:ImpBSc){reference-type="eqref"
reference="eq:ImpBSc"} in this case is $$\begin{equation}\label{eq:imp1}
k z = f(t)\end{equation}$$
We want to determine the modal frequencies and modal shapes
of the string under this condition. Remembering that $z(t) = y(L,t)$,
and remembering that the modal shapes for fixed boundary at $x=0$ are in
the form of $\sin \frac{\omega x}{c}$, condition
[\[eq:imp1\]](#eq:imp1){reference-type="eqref" reference="eq:imp1"}
gives
$$k y(L,t) = -T_0 \frac{\partial y}{\partial x}(L,0) \implies k \sin \frac{\omega L}{c} = - T_0 \frac{\omega}{c} \cos \frac{\omega L}{c}$$
First, two sanity checks. It is easy to see that, for $k=0$, we get back
to the fixed-free case ($\cos \frac{\omega L}{c} = 0$) analysed before.
Also, if $k\rightarrow \infty$, we obtain again the fixed-fixed
condition ($\sin \frac{\omega L}{c} = 0$). If $k \neq 0$, the boundary
condition gives $$\begin{equation}\label{eq:bc1}
\tan \frac{\omega L}{c} = - \frac{T_0}{ck}\omega\end{equation}$$
This is a
trascendental equation, that can be solved using numerical root-finding
algorithms. A picture of the solution to this equation is given in
Figure [1.11](#fig:trasc){reference-type="ref" reference="fig:trasc"}.
Notice that, when $\omega \rightarrow \infty$, the equation reduces to
$\tan\frac{\omega L}{c}\rightarrow \infty$ which is solved by
$\omega_m = (2m-1)\pi c / 2L$, i.e. the same frequencies as the
fixed-free case. In practice, at high frequencies the effects of
stiffness are negligible compared to the case of free boundary. The most
affected modes are at low frequencies!

Case b): $k=0,R=0$. A second case of interest is given by a bridge where
the mass is much larger than the stiffness and loss parameter. The
boundary condition [\[eq:ImpBSc\]](#eq:ImpBSc){reference-type="eqref"
reference="eq:ImpBSc"} in this case is $$m \frac{d^2z}{dt^2} = f(t)$$
Thus,
$$-m \omega^2 y(L,t) = -T_0 \frac{\partial y}{\partial x}(L,0) \implies m  \sin \frac{\omega L}{c} =   \frac{T_0}{c\omega} \cos \frac{\omega L}{c}$$
Another sanity check here indicates that, when $m=0$ we obtain again the
solution to the fixed-free boundary; when $m\rightarrow \infty$ we
obtain again the frequencies of the fixed-fixed boundary (as expected!
If the mass is really small / large, the boundary is effectively free /
fixed.) For all other cases, one has $$\begin{equation}\label{eq:bc2}
\tan \frac{\omega L}{c} = \frac{T_0}{c m}\omega^{-1}\end{equation}$$
which is another
trascendental equation solvable via numerical root-finding procedures.
Notice again that, for large frequencies, we obtain
$\tan\frac{\omega L}{c}= 0$ which is solved by $\omega_m = m\pi c / L$:
at high frequencies the bridge appears rigid! See also Figure
[1.11](#fig:trasc){reference-type="ref" reference="fig:trasc"}.

## Viscous Loss

Conservative wave propagation is in practice never realised, at least
over time scales of the order of seconds (i.e. the timescale of a
musical signal, for example.) Thus, we are now going to study the
effects that losses have on wave propagation. In analogy with the
oscillator, we may introduce a loss factor propotional to the velocity
of the string. Hence $$\begin{equation}\label{eq:WEloss}
\frac{\partial^2 y}{\partial t^2} = c^2 \frac{\partial^2 y}{\partial x^2} - 2\sigma   \frac{\partial y}{\partial t}$$
Note that the energy analysis here leads to $$\begin{equation}\label{eq:EnBalWE}
\frac{d}{dt}\left( \int_{0}^{L} \frac{1}{2}\left( \frac{\partial y}{\partial t} \right)^2 dx + \int_{0}^{L} \frac{c^2}{2}\left( \frac{\partial y}{\partial x} \right)^2 dx \right) = -2\sigma \int_0^L \left(\frac{\partial y}{\partial t}\right)^2 dx,$$
and hence energy is not increasing whenever $\sigma \geq 0$, which will
be assumed in the remainder.

A first thing to notice here, is that the D'Alembert solution ceases to
be valid when $\sigma  \neq 0$. Thus, we cannot use such solution in the
current analysis. The other tools that we developed previously though
(separation of variables, analysis in the Fourier domain) are still
valid, and we are going to make use of them now. What we could do, is to
try to separate the solution into components depending only on time and
space, i.e. $$y(x,t) = T(t)X(x)\end{equation}$$
We could then use this expression into
[\[eq:WEloss\]](#eq:WEloss){reference-type="eqref"
reference="eq:WEloss"}, and follow the same steps as we did for the wave
equation without loss, in last week's lecture. The net result is that we
end up with
$$y(x,t) = T(t)X(x) = \left(A_+ e^{st} + A_- e^{-st} \right)\left(B_+ e^{j\gamma x} + B_- e^{-j \gamma x} \right)$$
where the four complex constants $A_+,A_-,B_+,B_-$ are arbitrary so long
as initial and boundary conditions are not specified. Above,
$s = -\sigma + j \omega$. In practice, the solution is a combination of
typical terms such as $$\begin{equation}\label{eq:ddp}
y(x,t) = e^{st}e^{j\gamma x}\end{equation}$$
In actual fact, when analysing PDEs such
as [\[eq:WEloss\]](#eq:WEloss){reference-type="eqref"
reference="eq:WEloss"}, we may immediately try to use solutions such as
[\[eq:ddp\]](#eq:ddp){reference-type="eqref" reference="eq:ddp"}. We
make an *assumption* that solutions have that form, and then
characterise the solutions according to the particular relationship
between the frequency $\omega$ and wavenumber $\gamma$. Such assumption
in textbooks is often called an *ansatz* (a german word.)

The use of solutions such as
[\[eq:ddp\]](#eq:ddp){reference-type="eqref" reference="eq:ddp"} in PDEs
is one of the tools you should understand and master in this course. The
reason why we want to use such expression for the solution, as noted, is
to find a relation between $\omega$ and $\gamma$ known as the
*dispersion relation*. Doing so for
[\[eq:WEloss\]](#eq:WEloss){reference-type="eqref"
reference="eq:WEloss"} gives
$$s^2 + 2 \sigma s + c^2 \gamma^2 = 0 \implies s = - \sigma \pm j \sqrt{c^2\gamma^2 - \sigma^2}$$
Thus, the dispersion relation is $$\begin{equation}\label{eq:dispRelLoss}
\omega = \sqrt{c^2\gamma^2 - \sigma^2}\end{equation}$$
Of course, when $\sigma=0$ one
recovers the dispersion of the simple wave equation $\omega = c\gamma$.
The interesting fact about
[\[eq:dispRelLoss\]](#eq:dispRelLoss){reference-type="eqref"
reference="eq:dispRelLoss"} is that now the system is *dispersive*.
Dispersion is a property of vibratory system in which waves at different
frequencies travel at different speeds. The wave equation without loss
is *not* dispersive, because all waves travel at the same speed $c$. The
wave equation with loss, on the other hand, has dispersion. Consider a
simple harmonic solution $$y(x,t) = \cos(\omega t - \gamma x)\end{equation}$$
Using
the definition of *phase velocity* $v_\phi$ (encountered in last week's
lecture), one has
$$v_\phi = \frac{\omega}{\gamma} = \frac{\sqrt{c^2 \gamma^2 - \sigma^2}}{\gamma} \approx c - \frac{\sigma^2}{2c\gamma^2} = c \left(1 - \frac{\sigma^2}{2(\omega^2+\sigma^2)} \right)$$
and thus we see that the phase velocity now depends on the frequency of
the wave! In particular, waves of high frequency are non-dispersive, as
the limit of the dispersion relation as $\omega \rightarrow \infty$ is
$c$.

## Dispersion. Group Velocity

The concept of dispersion is fairly important for our study on wave
propagation, because it brings about the question of what we mean by
*wave velocity*. For a non-dispersive system such as the simple lossless
wave equation, the answer is trivial as all wavefronts travel with the
same speed $c$. This allows to write solutions in the form given by
D'Alembert. However, the introduction of a viscous loss term has the
consequence of introducing a frequency-dependent phase velocity. So,
what is the speed of a wave?

It is convenient, in dispersive systems, to not consider sines or
cosines in isolation, but to consider *groups* of them. This is because,
from Fourier theory, a disturbance may always be decomposed onto a
linear combination of sines and cosines. Thus, let such disturbance be
$$\begin{equation}\label{eq:gvv}
y(x,t) = \sum_{i=1}^N \cos(\omega_i t - \gamma_i x + \phi_i)\end{equation}$$
(here,
for simplicity, we are considering a discrete sum, though we may replace
the sum with a Fourier integral as seen before.) The form of $y(x,t)$ is
that of a group of waves. We are also making the assumption that such
group is centered around some reference frequency $\omega_0$ and
wavenumber $\gamma_0$, so that effectively
$$\omega_i = \omega_0 + \delta \omega_i, \quad \gamma_i = \gamma_0 + \delta \gamma_i, \quad i \in [1,N]$$
with $\delta \omega_i$, $\delta \gamma_i$ being small. Thus, the change
in phase $dP_i$ for a small displacement $dx$ in the time $dt$ for each
of the waves is given by $$dP_i = \omega_i dt - \gamma_i dx\end{equation}$$
Now, this
is the very important point: if the group moves "without changing in
shape too much" then *all* the phases in the sum
[\[eq:gvv\]](#eq:gvv){reference-type="eqref" reference="eq:gvv"} have
changed by the same amount. Thus we set $dP_i - dP_j \approx 0$, or
$$(\omega_i-\omega_j)dt - (\gamma_i-\gamma_j) \approx 0 \implies \frac{dx}{dt} \approx \frac{\delta \omega}{\delta \gamma}$$
Thus, is convient to define the *group velocity*
$$v_{g} = \frac{d \omega }{d \gamma}\end{equation}$$
that has the interpretation of
the velocity of a group of waves, centered around some mean values
$\omega_0$, $\gamma_0$, that travel with mostly the same speed, so that
a packet of such waves retains its shape whilst travelling through the
medium.
