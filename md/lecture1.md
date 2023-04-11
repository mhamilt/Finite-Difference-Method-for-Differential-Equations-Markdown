---
layout: post
title: Lecture 1
---

$$
\newcommand{\dif}{\text{d}}
\newcommand{\al}{\alpha}
\newcommand{\bt}{\beta}
\newcommand{\om}{\omega}
\newcommand{\gm}{\gamma}
\newcommand{\lb}{\lambda}
\newcommand{\lm}{\lambda_-}
\newcommand{\lp}{\lambda_+}
\newcommand{\thw}{\theta_w}
\newcommand{\thp}{\theta_\phi}
\newcommand{\tht}{\theta}

\newcommand{\etp}{e_{t+}}
\newcommand{\etm}{e_{t-}}

\newcommand{\esp}{e_{x+}}
\newcommand{\esm}{e_{x-}}

\newcommand{\dsp}{\delta_{x+}}
\newcommand{\dsm}{\delta_{x-}}
\newcommand{\dsd}{\delta_{x\cdot}}
\newcommand{\dss}{\delta_{xx}}
\newcommand{\dssss}{\delta_{xxxx}}

\newcommand{\dtp}{\delta_{t+}}
\newcommand{\dtm}{\delta_{t-}}
\newcommand{\dtd}{\delta_{t\cdot}}
\newcommand{\dtt}{\delta_{tt}}
\newcommand{\dtttt}{\delta_{tttt}}
\newcommand{\dspm}{\delta_{x\pm}}

\newcommand{\mtp}{\mu_{t+}}
\newcommand{\mtm}{\mu_{t-}}
\newcommand{\mtd}{\mu_{t\cdot}}
\newcommand{\mtt}{\mu_{tt}}
$$

- [Introduction](#introduction)
  - [Newton's law for mechanical oscillations](#newtons-law-for-mechanical-oscillations)
    - [Energy Analysis {#sec:EnAnGen}](#energy-analysis-secenangen)
    - [Bounds on solution growth {#sec:BoundsSHO}](#bounds-on-solution-growth-secboundssho)
    - [Periodicity of orbits via Fourier series](#periodicity-of-orbits-via-fourier-series)
  - [Time Difference Operators](#time-difference-operators)
    - [Shift, time and averaging operators](#shift-time-and-averaging-operators)
  - [Frequency domain analysis and stability of LTI systems {#sec:FreqDomAn}](#frequency-domain-analysis-and-stability-of-lti-systems-secfreqdoman)
    - [Laplace and $z$ transforms](#laplace-and-z-transforms)
    - [Fourier and discrete time Fourier transforms](#fourier-and-discrete-time-fourier-transforms)
    - [Frequency-domain intepretation of time difference operators {#sec:FDtransformations}](#frequency-domain-intepretation-of-time-difference-operators-secfdtransformations)
    - [Recursion polynomials](#recursion-polynomials)


# Introduction

This course is about the analysis and simulation of mechanical
vibrations. Since a kind of oscillation is involved, the systems of
interest here are time-dependent These systems appear in many branches
of the sciences, including mechanical engineering, acoustics,
vibroacoustics, and others. In this chapter, Newton's law for a system
comprising one degree of freedom will be reviewed, along with useful
concepts used in later chapters. Discrete calculus will also be
introduced.

## Newton's law for mechanical oscillations

We shall begin the discussion by studying the evolution of a simple
vibrating object, comprising a mass $m$ and a subjected to a force field
originating from a suitable potential,

$$\begin{equation}\label{eq:PhiF}
    F = - \frac{d\phi}{dx}.
\end{equation}$$

Here
$\phi = \phi(x): \mathbb{R} \rightarrow \mathbb{R} \in \mathcal{C}^1$ is
the potential (assumed to be differentiable), and
$x = x(t): \mathbb{R}^+_0 \rightarrow \mathbb{R} \in \mathcal{C}^2$ is
the displacement (assumed to be differentiable twice), measured from
some convenient reference position. Here, $t \geq 0$ is time. Newton's
law gives

$$\begin{equation}\label{eq:SHO}
    m \frac{d^2 x}{dt^2} = - \frac{d\phi}{dx}.
\end{equation}$$

Equation
\eqref{eq:SHO} is an
example of *ordinary differential equation* (ODE). The equation is
*autonomous*, intended as the absence of explicit dependence on time
$t$, except as via the argument of $x$. The equation is, in general,
*nonlinear*. Linearity is expressed here as a superposition principle:
if $x_1(t)$ and $x_2(t)$ are both solutions to
\eqref{eq:SHO}, then,
under linear conditions, $x_3(t) = a_1 x_1(t) + a_2 x_2(t)$ (with $a_1$,
$a_2$ being constants) is also a solution. Since the second time
derivative on the left-hand side of
\eqref{eq:SHO} is
linear, linearity is obtained when
$$\frac{d\phi(x_3)}{d x_3} = a_1 \frac{d\phi(x_1)}{dx_1} + a_2 \frac{d\phi(x_2)}{dx_2}$$

which is solved for $\phi = c x^\alpha$, where $c$ is a constant, and
$\alpha \in \{0,2 \}$ (the case $\alpha =0$ being the trivial case of
zero force).

To be complete, \eqref{eq:SHO} requires the specification of two *initial
conditions*, usually given as the initial displacement and velocity,
such that

$$\begin{equation}\label{eq:SHO_ICs}
    x(t=0) = x_0, \,\, \frac{dx}{dt}(t=0) = v_0,
\end{equation}$$

and since the
independent variable is the time $t$,
\eqref{eq:SHO} plus
\eqref{eq:SHO_ICs} specify an *initial value problem* (IVP). When
initial conditions are imposed, the system possesses one *unique
solution*, if appropriate Lipschitz-continuity arguments are satisfied.
(There is no need to study these in detail, we will just assume that the
examples provided here have a unique solution).

### Energy Analysis {#sec:EnAnGen}

As anticipated, all through this course, we will rely heavily on energy
arguments. The fact that systems present some kind of energy balance is
a fundamental principle of physics, bearing significant consequences in
the analysis of both the continuous systems, and the numerical
approximations used to simulate them. In the simple, one dimensional
case \eqref{eq:SHO},
the work $W$ done by a force pushing on $m$ is
$$W = \int_0^x F \, \text{d} x
$$

where it is assumed that the initial
position of the body, at the time $t=0$, is 0. Using this definition in
\eqref{eq:SHO} gives
$$\int_0^x m \frac{d^2 x}{dt^2} \text{d} x = - \int_0^x \frac{d\phi}{dx} \text{d} x.$$

To integrate, one uses $dx = \frac{dx}{dt} dt$, giving

$$\begin{equation}\label{eq:EnAnGen}
    \int_0^t m \frac{d^2 x}{dt^2} \frac{d x}{dt} \text{d} t = - \int_0^t \frac{d\phi}{dx} \frac{dx}{dt} \text{d} t
\end{equation}$$

and, using simple identities, one gets

$$\begin{equation}\label{eq:En1}
    \int_0^t \frac{d}{dt}\left( \frac{m}{2}\left(\frac{dx}{dt}\right)^2 + \phi  \right) \text{d} t = 0.\end{equation}$$

Since $x(t) \in \mathcal{C}^2$, $\phi(x) \in \mathcal{C}^1$, the
equation is solved by taking

$$\begin{equation}\label{eq:En2}
    \frac{m}{2}\left(\frac{dx}{dt}\right)^2 + \phi = H_0,
\end{equation}$$

where $H_0$
(a constant) is the total energy of the system. It is convenient to
identify the *kinetic* and *potential* components of the energy, given,
respectively, by
$$E_k \triangleq \frac{m}{2}\left(\frac{dx}{dt}\right)^2, \quad E_p \triangleq \phi.$$

The expression for $H_0$ is determined by the initial conditions, and
hence

$$\begin{equation}\label{eq:EnCons}
    H_0 = \frac{m v_0^2}{2} + \phi(x_0),
\end{equation}$$

and the identity $H(t) = H_0$
holds $\forall t \geq 0$.

### Bounds on solution growth {#sec:BoundsSHO}

Energy, as seen, is conserved. It is remarked that the energy is also
*non-negative*, i.e. $H(t) \geq 0$ $\forall t$, when $\phi$ itself is
non-negative. This fact leads to the important result of *boundedness of
the solutions*. In practice, since the kinetic and potential energies
are *both* non-negative, one has

$$\begin{equation}\label{eq:bnds}
    0 \leq E_k \leq H_0, \quad 0\leq E_p\leq H_0.
\end{equation}$$

For the kinetic
term, one has, in all cases,

$$\begin{equation}\label{eq:EgyBoundVel}
    \left|\frac{dx(t)}{dt}\right| \leq \sqrt{2 H_0/m}\,\,\, \forall t,\end{equation}$$

and thus the velocity of the system is *always* bounded in terms of the
initial energy. The displacement $x(t)$ itself may or may not be
bounded. Consider the case of a quadratic potential, as in the left
panel of Fig. [1.1](#fig:phis){reference-type="ref"
reference="fig:phis"}: $\phi = \frac{K x^2}{2}$. This case (as shown
above) corresponds to the linear case, i.e. the simple harmonic
oscillator with stiffness constant $K$. From
\eqref{eq:bnds}, one
has $|x|\leq \sqrt{2 H_0/K}$. Hence, the displacement is bounded. When
the potential is not quadratic, a nonlinear system is obtained. We shall
consider the nonlinear oscillator later on. Here, qualitatively, two
nonlinear potentials are shown in the center and right panels of Fig.
[1.1](#fig:phis){reference-type="ref" reference="fig:phis"}. Bounded
motion is obtained whenever the potential attains infinity for large
values of $|x|$. However, if the potential $\phi$ is bounded by some
finite constant, unbounded motion can be observed.

![Potential functions corresponding to a quadratic function
$\phi  \propto x^2$ (left panel), and non-quadratic functions (center
and right panels). Total energy is represented as dashed or dash-dotted
lines, the latter corresponding to unbounded motion.
](figures/phiGraphs.png){#fig:phis}

To understand these claims, it can be useful to visualise the
trajectories in the *phase plane*. This is a plane whose axes are
$x = x(t)$ and $y = dx(t)/dt$, see Fig.
[\[fig:phasePorts\]](#fig:phasePorts){reference-type="ref"
reference="fig:phasePorts"}. For the quadratic case, trajectories are
ellipses, and are hence symmetric about both the $x$ and $y$ axes.
Periodicity is evident by inspection of the phase portraits, which
appear as closed loops. The energy components also oscillate
periodically, in this case at one single frequency. For the central
panel of Fig. [1.1](#fig:phis){reference-type="ref"
reference="fig:phis"}, trajectories are now symmetric only about the $x$
axis: the nonlinearity is such that symmetry is broken in $x$, as
evident from panel 2 of Fig. [1.1](#fig:phis){reference-type="ref"
reference="fig:phis"}. However, motion is still periodic, via the
combination of a number frequencies with a common factor (the
fundamental frequency). For the right panel of Fig.
[1.1](#fig:phis){reference-type="ref" reference="fig:phis"}, motion may
be bounded or unbounded, depending on the value of $H_0$, as seen in
panel 3 of Fig.
[\[fig:phasePorts\]](#fig:phasePorts){reference-type="ref"
reference="fig:phasePorts"}.


When motion is periodic, the period can be estimated from the energy,
since

$$\begin{equation}\label{eq:dxdtEn}
    \frac{dx}{dt} = \pm\sqrt{\frac{2}{m}}\sqrt{H_0-\phi}.
\end{equation}$$

The plus or
minus signs depends on which side of the phase portrait the particle is
(i.e. whether it is found in the $y\geq0$ or $y<0$ half-plane).
Considering the plus sign, and inverting, one has
$$\text{d} t = \sqrt{\frac{m}{2}}\frac{\text{d} x}{\sqrt{H_0-\phi}}.
$$

Thus,
half the period is obtained integrating between $x_1$ and $x_2$, which
are the points in the phase plane where $dx/dt = 0$. Denoting the period
$\tau$, one has
$$\tau = \sqrt{2m} \int_{x_1}^{x_2}\frac{\text{d} x}{\sqrt{H_0-\phi}}.
$$

In
general, this equation does not have a closed-form solution. Some cases
are exceptional, including the case of simple harmonic motion. For that,
$x_1 = -x_2$, and $H_0=H(t) =\frac{K x_2^2}{2}$, where the last equality
holds since $dx/dt|_{x=x_2} = 0$. Hence, the period in obtained as the
integral over one quadrant (one quarter of a loop):
$$\tau = 2\sqrt{2m} \int_0^{x_2}\frac{\text{d} x}{\sqrt{\frac{K x_2^2}{2}-\frac{K x^2}{2}}} = 4 \sqrt{\frac{m}{K}} \arctan\left( \frac{x}{\sqrt{x_2^2-x^2}}\right)\big|_0^{x_2} = 2\pi \sqrt{\frac{m}{K}}.$$

It is convenient to introduce the *radian frequency*
$\omega_0 = \sqrt{K/m} = 2\pi f_0$, where $f_0$ is a linear frequency
(in Hz). Hence, one has

$$\begin{equation}\label{eq:tau}
    \tau = \frac{2\pi}{\omega_0} = \frac{1}{f_0}.\end{equation}$$

### Periodicity of orbits via Fourier series

For Case 1 of the previous section (quadratic potential), periodicity of
the solution is expressed as

$$\begin{equation}\label{eq:Fou1}
    x(t) = a \cos(\omega_0 t) + b \sin(\omega_0 t) = A \cos(\omega_0 t - \varphi)  = C_+ e^{j\omega_0 t} + C_- e^{-j \omega_0 t}\end{equation}$$

where

$$\begin{equation}\label{eq:Fou2}
    A^2 = a^2 + b^2, \,\, a = A \cos\varphi, \,\, b = A\sin\varphi, \,\, C_+ = \frac{1}{2}(a+jb), \,\, C_- = \frac{1}{2}(a-jb)\end{equation}$$

Note that the complex exponential notation in
\eqref{eq:Fou1}
yields indeed a real solution, when $C_+,C_-$ are selected as in
\eqref{eq:Fou2}. The
expression for $x(t)$ is periodic, with period $\tau$ given in
\eqref{eq:tau}.

When the potential is not quadratic, the dynamics is nonlinear, as
observed in Case 2 and (partly) 3 of the previous subsection. Since
periodic motion still exists, this can be obtained as a combination of
frequencies, all multiples of a fundamental frequency. A generalisation
of \eqref{eq:Fou1}
can then be given as

$$\begin{equation}\label{eq:Fou3}
x(t) = \sum_{m=0}^M \left( a_m \cos(m \omega_0 t) + b_m \sin(m \omega_0 t) \right),\end{equation}$$

which can be turned into amplitude-phase or complex exponential
expressions analogous to those in
\eqref{eq:Fou1},
using identities similar to
\eqref{eq:Fou2}. The
upper bound $M$ in the sum is in theory infinity, though one chooses a
finite $M$ in any practical application. Expression
\eqref{eq:Fou3} is
called a *Fourier series* expansion, which is valid and uniquely
determined for any periodic signal $x(t)$, with period
$\tau = 2\pi / \omega_0$. If one happens to know $x(t)$ for a given
system (for example, via a measurement), then the Fourier components can
be extracted as

$$\begin{equation}\label{eq:FouComps}
a_0 = \frac{1}{\tau}\int_0^\tau x(t) \, \text{d} t, \,\, a_{m\neq 0} = \frac{2}{\tau}\int_0^\tau x(t)\cos(m\omega_0 t) \, \text{d} t, \,\, b_{m} = \frac{2}{\tau}\int_0^\tau x(t)\sin(m\omega_0 t) \, \text{d} t\end{equation}$$

These expressions are proven easily when considering the *orthogonality*
of the Fourier components (easily obtained via direct integration), i.e.

$$
\begin{aligned}
\int_{0}^\tau \cos(m\omega_0 t)\cos(n \omega_0 t) \text{d} t &= \frac{\tau}{2}\delta_{m,n} \,\,\, (m>0, n\geq 0) \label{eq:OrthoFou1} \\
\int_{0}^\tau \sin(m\omega_0 t)\sin(n \omega_0 t) \text{d} t &= \frac{\tau}{2}\delta_{m,n} \,\,\, (m>0, n\geq 0) \label{eq:OrthoFou2} \\
\int_{0}^\tau \cos(m\omega_0 t)\sin(n \omega_0 t) \text{d} t &= 0  \label{eq:OrthoFou3}
\end{aligned}
$$


Multiplying \eqref{eq:Fou3} by, say, $\sin(n\omega_0 t)$, integrating, and
using \eqref{eq:OrthoFou2}, one obtains the last expression in
\eqref{eq:FouComps}. The other expressions in
\eqref{eq:FouComps} are obtained analogously. A picture of the
Fourier components for Cases 1,2 and 3 is given in Fig.
[\[fig:Fou1\]](#fig:Fou1){reference-type="ref" reference="fig:Fou1"},
where one sees that simple harmonic motion is indeed characterised by a
single frequency of vibration.

![image](figures/Fou1.png){width="0.33\\linewidth"}
![image](figures/Fou2.png){width="0.33\\linewidth"}
![image](figures/Fou3.png){width="0.33\\linewidth"}

It may be useful to express the Fourier series in complex exponential
form. In that case, one has

$$\begin{equation}\label{eq:FourierSeriesComplex}
x(t) = \sum_{m=-M}^M c_m e^{jm\omega_0 t},
\end{equation}$$

where $c_0 = a_0$,
$c_m = \frac{1}{2}\left( a_m - j b_m\right)$ for $m>0$,
$c_m = \frac{1}{2}\left( a_m + j b_m\right)$ for $m<0$.

## Time Difference Operators

This section introduces the notation and the principles of difference
calculus, upon which the method of *finite differences* is constructed.
Though finite differences are used as a method for computer-aided
solution of differential equations, this method has surprisingly long
roots, reaching as far as the work of Courant, Friedrichs, and Lewy in
1928 (before digital computers were even available). The underlying
principle of finite differences is straightforward, and it involves the
discrete approximation of differential operators. In digital
applications, whether e.g. performing a measurement, or in
computer-aided simulation, time is most often discretised by means of a
*sample rate*, $f_s$. In practice, one usually knows (or is interested
in knowing) the state of a system at discrete time intervals, of length
$k$ seconds, the time step. The relationship between sample rate and
time step is simply $$k f_s = 1,
$$

The aim of finite differences is to
compute a *time series* $x^n$ approximating the "true" solution $x(t)$
of a given model problem. The index $n \in \mathbb{N}_0$ in $x^n$ is
shorthand for $t = t_n \triangleq nk$, i.e. $n$ is the time index,
approximating the true solution $x(t)$ at the time $t_n = nk$. The
fundamental relationship between the approximate time series $x^n$ and
the true solution $x(t)$ is as follows:

$$\begin{equation}\label{eq:ErrDef}
    x(t_n) - x^n  = E^n ,
\end{equation}$$

where $E^n$ is a time series defining the
*absolute error* of the approximation. In general, $E^n \neq 0$, and,
generally, it will not be possible to obtain an exact expression for
$E^n$ (since $x(t)$ is generally unknown!) However, the study of finite
differences is almost entirely devoted to the design of schems for which
the error remains provably *bounded* by a given power of $k$, and we
shall of course spend considerable effort in studying such schemes.

### Shift, time and averaging operators

Given the time series $x^n$, the identity, forward and backward shift
operators are given as
$${1}{x}^n = { x}^n, \quad \etp {x}^n = { x}^{n+1}, \quad \etm { x}^n = { x}^{n-1}$$

From these, one may define the time difference operators, all
approximating the first time derivative, as

::: subequations
$$\begin{aligned}
        \dtp & = \frac{(\etp - 1)}{k} \approx \frac{d}{dt},     \\
        \dtm & = \frac{(1 - \etm)}{k}\approx \frac{d}{dt},      \\
        \dtd & = \frac{(\etp - \etm)}{2k}\approx \frac{d}{dt} .
    \end{aligned}$$

:::

An approximation to the second time derivative is constructed from the
above as $$\delta_{tt}  = \dtp\dtm\approx \frac{d^2}{dt^2}.
$$

Averaging
operators (all approximating the identity) are also used throughout the
text, and are

$$
\begin{aligned}
        \mtp & = \frac{(\etp + 1)}{2} \approx 1,    \\
        \mtm & = \frac{(1 + \etm)}{2} \approx 1,    \\
        \mtd & = \frac{(\etp + \etm)}{2} \approx 1.
\end{aligned}
$$


Whilst these expressions look at least reasonable, Taylor series
arguments can be used to infer the *order* of the approximation. Hence,
the difference operators are applied to the continuous function $x(t)$,
and a Taylor expansion is applied. For the forward difference operator,
one gets 
$$\begin{aligned}
    \dtp x(t_n) = \frac{x(t_n+k) - x(t_n)}{k} & \approx \frac{x(t_n)+ k \frac{dx(t)}{dt}|*{t=t_n}+\frac{k^2}{2}\frac{d^2x(t)}{dt^2}|*{t=t_n}-x(t_n)}{k} \\&= \frac{dx(t_n)}{dt}+\frac{k}{2}\frac{d^2x(t_n)}{dt^2}\end{aligned}
$$

Since ${d^2x(t_n)}/{dt^2}$ is a value independent of $k$, in the limit
of high sample rate the expression above reduces to ${dx(t_n)}/{dt}$.
The rate at which such approximation is satisfied is linear in $k$, so
that applying definition
\eqref{eq:ErrDef} one gets

$$\begin{equation}\label{eq:Errs1}
    \frac{dx(t_n)}{dt} - \dtp x(t_n) = O(k).
\end{equation}$$

The "big-Oh" notation
$O(k^p)$ means that the rate of the approximation goes as $k^p$. For the
current case, $p=1$. Using similar arguments, it is easy to show that

$$\begin{aligned}
     \frac{dx(t_n)}{dt} - \dtm x(t_n) &= O(k), \\
     \frac{dx(t_n)}{dt} - \dtd x(t_n) &= O(k^2),\\
     \frac{d^2x(t_n)}{dt^2} - \delta_{tt} x(t_n) &= O(k^2),
     \end{aligned}
     $$

and so on. These identities define the *truncation errors* of the finite
difference approximations. It is seen, then, that some operators have a
higher accuracy than others. Generally, these employ a larger *stencil*,
that is, the footprint in time of a finite difference operator including
the coefficients. So, $\dtm$, $\dtp$, $\mtm$, $\mtp$ all have a stencil
of width 1. $\dtd$, $\mtd$, $\delta_{tt}$ have a stencil of width 2. See also
Fig. [1.2](#fig:stencil){reference-type="ref" reference="fig:stencil"}.

It is of course possible to construct operators with higher accuracy. As
an example, consider the following central approximation to the first
time derivative,
$$\bar \delta_{t\cdot} = \frac{\etm^2 - 8 \etm + 8\etp - \etp^2}{12k}.$$

It is easy (though a little tedious) to show that this approximation is
$O(k^4)$. At this point of the discussion, one may be tempted to think
that the construction of higher-order schemes for the solution of a
model problem amounts to merely employing difference operators with the
appropriate stencil. Things are, of course, more complicated than this:
primarily, operators with wide stencils in time tend to yield *unstable*
simulations, as will be seen in forthcoming examples. Other problems
exist: for instance, it is not clear how to initialise an operator with
a stencil of width $M$, when the model problem is of order $N < M$. A
similar problem arises when discretising differential operators in space
for boundary-values problems, where one must set values for the "ghost
points" (i.e. points located outside the boundary of the grid).

Constructing higher-order schemes is of course possible, and sometimes
desirable, and these words of caution suffice for the moment as cursory
understanding of the underlying difficulties.

When solving differential equations, different definitions of errors are
employed: these are the *local truncation error* (LTE) and the *global
error*. The truncation errors given by
\eqref{eq:Errs1}
and \eqref{eq:Errs}
are only useful to determine the order of the approximation of a
difference operator, but ensuring that the differential operators of a
difference scheme are approximated to an order $p$, does not
automatically ensure that the scheme is $p^{th}$ order accurate, if in
fact convergent at all.

![Stencil width and coefficients (units of $k$) of various difference
operators approximating the first time
derivative](figures/stencils.pdf){#fig:stencil width="\\linewidth"}

## Frequency domain analysis and stability of LTI systems {#sec:FreqDomAn}

Traditional analysis techinques for time-dependent problems rely on
frequency domain analysis, both in the continuous and discrete cases.

### Laplace and $z$ transforms

For the continuous function $x(t)$, the Laplace transform is obtained as

$$\begin{equation}\label{eq:LapT}
    \hat x(s) = \int_{\{-\infty,0\}}^{\infty} x(t) e^{-st} \text{d} t \triangleq \mathcal{L}\{x\}(s),\end{equation}$$

where $s = j\omega + \sigma \in \mathbb{C}$ is a complex variable. The
lower bound in the integral means that the transform can be defined to
be two-sided (starting from $-\infty$), or one-sided (starting from 0),
the latter allowing the incorporation of initial conditions. The closely
related $z$ transform is a discrete version of
\eqref{eq:LapT},
i.e.

$$\begin{equation}\label{eq:ZT}
    \hat x(z) = \sum_{n = \{-\infty,0 \}}^{\infty} x^n z^{-n} \triangleq \mathcal{Z}\{x \}(z)\end{equation}$$

for a complex number $z \in \mathbb{C}$. In the analysis of linear, time
invariant (LTI) systems, both continuous and discrete, one usually
computes $s$ and $z$ from the given model problem, via direct
substitution of the transforms. Then, the stability of the underlying
system may be inferred by direct inspection of $s$ and $z$, as will be
stated below. Applying the Laplace transform to the time derivative of
$x(t)$ gives:

$$\begin{equation}\label{eq:LPtemp1}
    \mathcal{L}\left\{\frac{d x}{dt } \right\} = \int_{-\infty}^\infty\frac{dx(t)}{dt}e^{-st} \text{d} t = s \int_{-\infty}^{\infty} x(t) e^{-st} \text{d} t =  s\mathcal{L}\{x\},\end{equation}$$

where integration by parts was used, and where it was assumed that
$x(t)$ dies out fast enough as $t\rightarrow \pm \infty$. The shift
operators under the two-sided $z$ transform become

$$\begin{equation}\label{eq:LPtemp2}
    \mathcal{Z}\left\{e_{t\pm} x^n \right\} = \sum_{n = -\infty}^{\infty} x^{n\pm 1}z^{-n} = z^{\pm 1} \sum_{n = -\infty}^{\infty} x^{n\pm 1}z^{-(n\pm 1)} = z^{\pm 1} \mathcal{Z}\{x^n \}.\end{equation}$$

Using analogous arguments as for
\eqref{eq:LPtemp1} and
\eqref{eq:LPtemp2}, one can show that higher-order derivatives and
differences transform as

$$\begin{equation}\label{eq:TransLZ}
    \mathcal{L}\left\{\frac{d^p x}{dt^p } \right\} = s^p \mathcal{L}\{x\}, \quad \mathcal{Z}\left\{e^p_{t\pm} x^n \right\} = z^{\pm p} \mathcal{Z}\left\{x^n \right\}.\end{equation}$$

In practice, it is often useful to substitute simpler expressions than
the transforms, via an *ansatz*. So, one may employ the test solutions

$$\begin{equation}\label{eq:ansatz}
    x(t) = \hat x e^{st}, \quad x^n = \hat x z^n,
\end{equation}$$

for appropriate
constant complex amplitudes $\hat x$. It is immediate to verify that
derivatives and differences of these test solutions transform in the
same way as the derivatives and differences of the Laplace and $z$
transforms in \eqref{eq:TransLZ}. For this reason, for LTI systems, the use of
the test solutions yields ultimately the same qualitative analysis as
the substitution of the full transforms, and we shall make use of such
solutions accordingly. Since time $t$ (and, accordingly, the time index
$n$) increases toward infinity, boundedness of $x(t)$ and $x^n$ in
\eqref{eq:ansatz} can be achieved if and only if
$\text{Re}(s) \leq 0$, and $|z|\leq 1$, respectively. This is, in
essence, the idea of stability in the frequency domain.

The fact that the transforms of derivatives and differences amount to
mere multiplications in the frequency domain is a powerful analysis tool
allowing for great simplifications. In theory, once the problem is
solved in the frequency domain, one can transform back to the time
domain, by computing the inverse transforms. This approach, however,
presents some drawbacks. Most notably, inverse transforms are difficult
to compute, and closed-form solutions are available only in a few cases.
Second, these techniques do not generalise to nonlinear, time variant
systems (though some exceptions exisit, such as Volterra kernels).
Frequency domain techniques remain a popular analysis tool nonetheless,
since most nonlinear systems reduce to linear under suitable conditions,
and one may infer (at least qualitatively) some useful properties of the
model systems under study.

Laplace transforms (and inverses) may be difficult to compute, and
usually one resorts to table look-ups. A couple of useful transorm
pairs, used later on, are listed here as

$$\begin{equation}\label{eq:LaplTtable}
    \mathcal{L}\left\{e^{-ct}\sin(a t)u(t)\right\}(s) = \frac{a}{(s+c)^2 + a^2}, \quad \mathcal{L}\left\{e^{-ct}\cos(a t)u(t)\right\}(s) = \frac{s+c}{(s+c)^2 + a^2}.\end{equation}$$

The function $u(t)$ is the step function at time $t=0$ (i.e.
$u(t<0)=0, u(t\geq 0) = 1)$, and thus the identities above apply equally
to the two-sided and one-sided transforms.

### Fourier and discrete time Fourier transforms

When one considers $\sigma=0$ in
\eqref{eq:LapT}, and
$z = e^{j\omega k}$ in \eqref{eq:ZT} (in the two-sided form), the continuous *Fourier
transform* and the *discrete time Fourier transform* (DTFT) are
obtained. These are useful to compute the magnitude and the phase of the
solutions (i.e. the spectrum), and will also be used throughout. The
definitions of these transforms are as:
$$\hat x(\omega) = \int_{-\infty}^{\infty} x(t) e^{-j\omega t} \text{d} t \triangleq \mathcal{F}\{x\}(\omega),\qquad \hat x(\omega) = \sum_{n = -\infty}^{\infty} x^n e^{-j\omega k n} \triangleq \mathcal{X}\{x^n \}(\omega).$$

The notation was used in both definitions, but the meaning is different.
We remark, however, that the discrete time Fourier transform is a
continuous function of the radian frequency $\omega$, and is thus not to
be confused with the discrete Fourier transform (which is evaluated at
discrete frequency bins). A couple of useful transform pairs are given
here as

$$\begin{equation}\label{eq:FouTtable}
    \mathcal{F}\left\{x(t)\sin(at)\right\}(\omega) = \frac{\hat x(\omega-a)-\hat x(\omega+a)}{2j}, \quad \mathcal{F}\left\{e^{-ct}u(t)| c \geq 0\right\}(\omega) = \frac{1}{\sqrt{2\pi}(c+j\omega)}.\end{equation}$$

In the second transform, $u(t)$ is again the step function as in
\eqref{eq:LaplTtable}.

### Frequency-domain intepretation of time difference operators {#sec:FDtransformations}

The action of the difference operators on a time series $x^n$ may be
interpreted in the $z$ domain as a *transformation*. In
\eqref{eq:LPtemp2}, it was seen that
$\mathcal{Z}\left\{e_{t\pm} x^n\right\} = z^{\pm 1} \mathcal{Z}\left\{x^n\right\}$.
The action of other difference operators may be constructed from such
identity. To evaluate the frequency response, $z$ is limited to the unit
circle, i.e. $z=e^{j\omega k}$, i.e. the DTFT
$\mathcal X\left\{x^n \right\}$. Consider, for example, the second
difference operator. This is

$$\begin{equation}\label{eq:DTFTdtt}
\mathcal{X}\left\{\delta_{tt} x^n\right\} = \frac{e^{j\omega k}-2+e^{-j\omega k}}{k^2}\mathcal{X}\left\{x^n\right\}=\frac{2}{k^2}\left(\cos(\omega k)-1 \right)\mathcal{X}\left\{x^n\right\} = -\frac{4}{k^2} \sin^2\left(\frac{\omega k}{2}\right)\mathcal{X}\left\{x^n\right\}\end{equation}$$

Analogously, one has

::: subequations
$$\begin{aligned}
\mathcal{X}\left\{\mtd x^n\right\} &= \cos(\omega k) \mathcal{X}\left\{x^n \right\}, \\
\mathcal{X}\left\{\mtt x^n\right\} &= \frac{1}{2}\left(\cos(\omega k) + 1\right)\mathcal{X}\left\{x^n \right\} \\
\mathcal{X}\left\{\dtd x^n\right\} &= \frac{j}{k}\sin(\omega k) \mathcal{X}\left\{x^n \right\}\end{aligned}$$

:::

These identities will prove useful in the study of the stability of LTI
systems via *von Neumann* analysis.

### Recursion polynomials

![The Schur-Cohn stability region, i.e. the region of the $(b,c)-$plane
for which the roots of
\eqref{eq:TempSchurCohn} have magnitude less than
unity.](figures/SchurCohn.pdf){#fig:SchurCohn width="0.3\\linewidth"}

Finite difference schemes, as will be seen shortly, produce a recursion
in time. When such recursion has constant coefficients, substitution of
*ansatz* \eqref{eq:ansatz} results in a complex equation of the form
$$\sum_{n=0}^M a_n z^n = 0,
$$

where the $a_n$'s are constant
coefficients, and $M$ is the problem order. Often, for the problems of
interest here, $M=2$, and thus the summation above reduces to (after
rescaling by $a_2 \neq 0$)

$$\begin{equation}\label{eq:TempSchurCohn}
    z^2 + b z + c = 0.
\end{equation}$$

The solutions (i.e. the poles) are given by
$z_\pm = \frac{-b \pm \sqrt{b^2-4 c}}{2}$. Useful bounds on the
coefficients $b,c$ can be derived when one considers $|z_\pm|<1$ (i.e.
when the poles are within the unit circle). These are

$$\begin{equation}\label{eq:SchurCohnStab}
    |c| = |z_+ z_-| =  |z_+| |z_-| < 1, \quad |b|<1+|c|,
\end{equation}$$

which are
necessary and sufficient conditions to enforce $|z_\pm|<1$, and are
known as *Schur-Cohn stability test*, see also Fig.
[1.3](#fig:SchurCohn){reference-type="ref" reference="fig:SchurCohn"}.
