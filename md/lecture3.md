---
layout: post
title: Lecture 3
---

- [Coupled Oscillators](#coupled-oscillators)
  - [Two masses, three springs](#two-masses-three-springs)
    - [Energy considerations](#energy-considerations)
    - [Eigenvalue decomposition](#eigenvalue-decomposition)
    - [Loss and Forcing](#loss-and-forcing)
    - [An Explicit Finite Difference Scheme](#an-explicit-finite-difference-scheme)
    - [Numerical Eigenvalues](#numerical-eigenvalues)
    - [A Family of Finite Difference Schemes](#a-family-of-finite-difference-schemes)
      - [Loss and Forcing](#loss-and-forcing-1)

# Coupled Oscillators

In this chapter we are going to study the problem of a system of coupled
oscillators, in both linear and nonlinear regimes. In the linear case,
extensions of the frequency-domain techniques already in use for the
oscillator in isolation are possible, leading to the idea of *modes of
vibration*. For the nonlinear case, analytic approximations based on
perturbation methods are possible, though the algebra quickly becomes
unwieldy. However, time-domain simulation techniques for this multimodal
case (and, particularly, energy methods) are in fact not too dissimilar
from those already encountered previously, and we shall make use of them
accordingly.

## Two masses, three springs

![Two-mass system. For small vertical displacements, the two masses move
vertically. The displacements from the rest positions are denoted $x_1$,
$x_2$.](Figures/TwoMassesGeneric.pdf){#fig:TwoMassesGeneric
width="\\linewidth"}

Consider the system sketched in Fig.
[1.1](#fig:TwoMassesGeneric){reference-type="ref"
reference="fig:TwoMassesGeneric"}. Here, the system possesses two
degrees of freedom, which may conveniently be identified with the
displacements $x_1(t)$, $x_2(t)$ of the two masses $m_1$, $m_2$, from
their rest position. The equations of motion for this system are as

::: subequations
[\[eq:TwoMass1\]]{#eq:TwoMass1 label="eq:TwoMass1"} $$\begin{aligned}
m_1 \frac{d^2 x_1}{dt^2} = -K_{11}x_1 - K_{12}(x_1 - x_2), \label{eq:Coupl1} \\
m_2 \frac{d^2 x_2}{dt^2} = -K_{22}x_2 + K_{12}(x_1 - x_2). \label{eq:Coupl2}\end{aligned}$$
:::

For the moment, it was assumed that motion is completely linear (the
amplitude of vibration is small, so that nonlinear effects can be
neglected).

The frequency-domain techniques described in Sec.
[\[sec:FreqDomAn\]](#sec:FreqDomAn){reference-type="ref"
reference="sec:FreqDomAn"} can be extended here. To that end, it is
convenient to adopt a matrix-vector formulation of the system. Though in
this example the system comprises $2$ masses, this formulation proves
advantageous since it can be used to describe systems comprising $N$
degrees of freedom. Furthermore, the matrix-vector formulation is the
natural framework to write the discrete space-time schemes that we will
encounter in later chapters. Hence,
[\[eq:TwoMass1\]](#eq:TwoMass1){reference-type="eqref"
reference="eq:TwoMass1"} is rewritten compactely as
$$\label{eq:MultiModalMatr}
{\bf M} \, \frac{d^2{\bf x}}{dt^2} = -{\bf K} \, {\bf x},$$ where
clearly
$${\bf M} = \begin{bmatrix}m_1 & 0 \\ 0 & m_2 \end{bmatrix}, \,\, {\bf K} = \begin{bmatrix}K_{11}+K_{12} & -K_{12} \\ -K_{12} & K_{22}+K_{12} \end{bmatrix}, \,\, {\bf x} = \begin{bmatrix}x_1  \\  x_2 \end{bmatrix}.$$
Extension of the frequency-domain techniques of Sec.
[\[sec:FreqDomAn\]](#sec:FreqDomAn){reference-type="ref"
reference="sec:FreqDomAn"}, one can conveniently define a test solution
of the form $${\bf x} =  {\bf X} e^{j\omega t}.$$ Here, ${\bf X}$ is a
constant complex amplitude, $j\omega$ was substituted for $s$, since for
lossless systems one has $\sigma = 0$ in $s = j\omega + \sigma$.
Substituting the test solution into
[\[eq:MultiModalMatr\]](#eq:MultiModalMatr){reference-type="eqref"
reference="eq:MultiModalMatr"}, one gets $$\label{eq:OmegaEigen}
\left( {\bf K } - \omega^2 {\bf M} \right)  {\bf X} = 0,$$ which shows
that the frequency $\omega$ is defined as the eigenvalue of the
eigenvalue problem
[\[eq:OmegaEigen\]](#eq:OmegaEigen){reference-type="eqref"
reference="eq:OmegaEigen"}. Nontrivial solutions are obtained when the
determinant of ${\bf K } - \omega^2 {\bf M}$ is zero. This defines the
polynomial whose roots yield the eigenvalues. For this two-dimensional
case, one has
$$(K_{11}+K_{12} - \omega^2 m_1)(K_{22}+K_{12} - \omega^2 m_2) - K_{12}^2 = 0.$$
Solving for $\omega^2$, one gets
$$\omega_\pm^2 = \frac{m_1(K_{22}+K_{12})+m_2(K_{11}+K_{12}) \pm \sqrt{\left(m_1(K_{22}+K_{12}) - m_2(K_{11}+K_{12})\right)^2 + 4 m_1 m_2 K_{12}^2 }}{2m_1 m_2}.$$
This solution shows that both the eigenvalues
$\omega^2_\pm \in \mathbb{R}$ (since the expression under the square
root sign is always positive); however, there is no guarantee that
$\omega^2_\pm$ are also *positive*. A negative $\omega^2$ would result
in exponentially growing behaviour of the test solution. Checking
positivity of the roots can be quite laborious in the general case. To
make things a little easier, we make the assumption $m_1=m_2=m$, for
which one has
$$2m\omega_\pm^2 = K_{11}+K_{22}+2K_{12} \pm \sqrt{\left( K_{11}-K_{22}  \right)^2 + 4  K_{12}^2}.$$
The positivity of *both* solutions is enforced when
$$\label{eq:TempIneq}
K_{11}+K_{22}+2K_{12} \geq 0, \,\, \text{ and } \,\, \left(K_{11}+K_{22}+2K_{12}\right)^2 \geq \left( K_{11}-K_{22}  \right)^2 + 4  K_{12}^2.$$
Solving for $K_{12}$, one has

::: subnumcases
K\_12 -       [\[eq:IneqK12\]]{#eq:IneqK12 label="eq:IneqK12"} K\_12 -
      K\_11+K\_22 \> 0, [\[eq:IneqK12a\]]{#eq:IneqK12a
label="eq:IneqK12a"}\
K\_12 -       K\_11+K\_22 \< 0. [\[eq:IneqK12b\]]{#eq:IneqK12b
label="eq:IneqK12b"}
:::

If instead $K_{11}+K_{22}=0$, at least one eigenvalue is surely
negative, as seen immediately from
[\[eq:TempIneq\]](#eq:TempIneq){reference-type="eqref"
reference="eq:TempIneq"}. As an example, consider the case
$K_{11}=K_{22}=1$. The conditions above give $K_{12}\geq -1/2$: this is
an interesting case, since one may allow the coupling to have negative
stiffness, whilst guaranteeing oscillating solutions overall. As a
second example, consider $K_{11}=K_{22}=-1$: here, there is no range
allowable for $K_{12}$.

Once the eigenvalues are computed, one may compute the eigenvector
${\bf X}$ from
[\[eq:OmegaEigen\]](#eq:OmegaEigen){reference-type="eqref"
reference="eq:OmegaEigen"}. Assuming for instance ${ X}_1 = a_{\pm}$
(where $a_\pm$ are just useful normalisation constants), one has
$$X_2 = a_{\pm}\frac{K_{11} + K_{12} - \omega^2_\pm m_1}{K_{12}},$$ such
that the general solution is given by
$${\bf x} = a_+\begin{bmatrix} 1 \\ \frac{K_{11} + K_{12} - \omega^2_+ m_1}{K_{12}}\end{bmatrix}\left(A_+ e^{j\omega_+ t} + A_- e^{- j\omega_+ t} \right) + a_-\begin{bmatrix} 1 \\ \frac{K_{11} + K_{12} - \omega^2_- m_1}{K_{12}}\end{bmatrix}\left(B_+ e^{j\omega_- t} + B_- e^{- j\omega_- t} \right),$$
where $A_\pm$, $B_\pm$ are four complex constants depending on the
intial conditions ${\bf x}(t=0)$, $\frac{d{\bf x}(t=0) }{dt}$. The
formula above is revealing: we showed that the motion of the system can
be written as the the sum of two harmonic motions, independent of each
other, one with frequency $\omega_+$, the other with frequency
$\omega_-$.

As an example, consider the case $K_{11},K_{22},K_{12},m_1,m_2=1$. In
this case, one has $\omega_+ = \sqrt{3}$, $\omega_- = 1$. Then
$$\label{eq:TwoMassesDecomposed}
{\bf x} = \frac{1}{\sqrt{2}}\begin{bmatrix} 1 \\ -1\end{bmatrix}\left(A_+ e^{j\sqrt{3} t} + A_- e^{- j\sqrt{3} t} \right) - \frac{1}{\sqrt{2}}\begin{bmatrix} 1 \\ 1\end{bmatrix}\left(B_+ e^{j  t} + B_- e^{- j t} \right),$$
Here the normalisation constants $a_+$, $a_-$ where chosen so that the
norm of the eigenvectors is 1. (Remember that these constants are
entirely arbitrary). If one chooses ${\bf x}(t=0)=[1,-1]^\intercal$,
$d{\bf x}(t=0)/dt=[0,0]^\intercal$, the constants are set as
$A_+=A_-=1/\sqrt{2}$, $B_+ = B_- = 0$, giving ultimately
$$\label{eq:Mode2}
{\bf x} = \begin{bmatrix} 1 \\ -1\end{bmatrix}\cos\left(\sqrt{3}t\right).$$
This shows that a two-mode system may collapse to a single mode, when
the system is started in that mode. There is no trace of the other mode
of vibration! One may of course start the system in the other mode, and
observe it oscillating in that mode only. Figs.
[1.2](#fig:TwoMassesMode1){reference-type="ref"
reference="fig:TwoMassesMode1"} and
[1.3](#fig:TwoMassesMode2){reference-type="ref"
reference="fig:TwoMassesMode2"} show the motion of the masses when the
system is started in either one of the two modes.

![Two-mass system. Visualisation of the first mode of
vibration.](Figures/TwoMassesMode1.png){#fig:TwoMassesMode1
width="0.75\\linewidth"}

![Two-mass system. Visualisation of the second mode of vibration,
corresponding to [\[eq:Mode2\]](#eq:Mode2){reference-type="eqref"
reference="eq:Mode2"}.](Figures/TwoMassesMode2.png){#fig:TwoMassesMode2
width="0.75\\linewidth"}

### Energy considerations

It is remarked that the allowable ranges for the stiffness constants
$K_{ij}$ are such that the potential energy of the system is
*non-negative*, so that one may bound the growth of the solutions in
some manner. The total energy of the system can be found by multiplying
[\[eq:Coupl1\]](#eq:Coupl1){reference-type="eqref"
reference="eq:Coupl1"} by $\frac{dx_1}{dt}$, and
[\[eq:Coupl2\]](#eq:Coupl2){reference-type="eqref"
reference="eq:Coupl2"} by $\frac{dx_2}{dt}$, and summing. Using the
usual identities, this gives $$\label{eq:EnBalCnt2Masses}
\frac{d}{dt}\left( \frac{m_1}{2}\left( \frac{dx_1}{dt} \right)^2 + \frac{m_2}{2}\left( \frac{dx_2}{dt} \right)^2 + \frac{K_{11} x_1^2}{2} +  \frac{K_{22} x_2^2}{2} + \frac{K_{12} (x_1-x_2)^2}{2}\right) \triangleq \frac{dH(t)}{dt} = 0,$$
showing that $H(t) = H(t=0) = H_0$, that is, energy is conserved. The
form of the energy comprises the energy of the two harmonic oscillators
in isolation, plus the coupling energy, proportional to $K_{12}$. The
coupling is a function of the relative distance between the masses, so
at any time such that $x_1-x_2=0$, the coupling force is null. It is
convenient to write the energy in matrix form, as
$$H(t) = \frac{1}{2}\left(\frac{d{\bf x}^\intercal}{dt} {\bf M} \frac{d{\bf x}}{dt} + {\bf x}^\intercal {\bf K} {\bf x}\right).$$
Since both ${\bf K}$, ${\bf M}$ are positive-definite[^1], one has
$$0 \leq \frac{d{\bf x}^\intercal}{dt} {\bf M} \frac{d{\bf x}}{dt} \leq H_0, \quad 0 \leq {\bf x}^\intercal {\bf K} {\bf x} \leq H_0,$$
that is, the norms of the solution remain bounded by some energy
constant incorporating the intial conditions. This is, in essence, the
idea of stability for this vector case, generalising bound
[\[eq:EgyBoundVel\]](#eq:EgyBoundVel){reference-type="eqref"
reference="eq:EgyBoundVel"} (and the analogous bound on $x$) of the
single scalar case.

### Eigenvalue decomposition

The discussion above suggests that the motion system
[\[eq:MultiModalMatr\]](#eq:MultiModalMatr){reference-type="eqref"
reference="eq:MultiModalMatr"} may in fact be decomposed onto linearly
independent blocks, called the *modes*. Since the mass matrix is usually
a diagonal matrix with positive entries, it is convenient to rewrite the
system as $$\label{eq:EigenDecTwoMasses1}
\frac{d^2{\bf x}}{dt^2} = -{ \bf M}^{-1}{\bf K} \, {\bf x},$$ where here
${{\bf M}^{-1}\bf K}$ is a positive-definite matrix. One is then able to
write $$\label{eq:eigendempos}
{\bf M}^{-1}{\bf K} = {\bf P} \, {\bf \Omega}^2 \, {\bf P}^{\intercal},$$
where ${\bf P}$ is a matrix comprising the column eigenvectors of
${\bf M}^{-1}{\bf K}$, and where ${\bf \Omega}$ is a diagonal matrix
containing the (positive) eigenvalues (that is, the resonant frequencies
of the system). Since ${\bf M}^{-1}{\bf K}$ is positive-definite, one
has that the matrix ${\bf P}$ is in fact *orthonormal*, that is,
${\bf P}^{-1}={\bf P}^\intercal$. Hence, multiplying
[\[eq:EigenDecTwoMasses1\]](#eq:EigenDecTwoMasses1){reference-type="eqref"
reference="eq:EigenDecTwoMasses1"} on the left by ${\bf P}^\intercal$
gives $$\label{eq:diagonalEigen}
\frac{d^2{\bf u}}{dt^2} = -{\bf \Omega}^2 \, {\bf u},$$ with
${\bf u} = {\bf P}^{\intercal}{\bf x}$. This is a completely diagonal
system, where the degrees of freedom are independent of each other.
Coming back again to the example
[\[eq:TwoMassesDecomposed\]](#eq:TwoMassesDecomposed){reference-type="eqref"
reference="eq:TwoMassesDecomposed"}, here one has
$${\bf P} =  \frac{1}{\sqrt{2}}\begin{bmatrix}-1 & -1 \\ -1 & 1 \end{bmatrix}, \,\, {\bf \Omega} =  \begin{bmatrix}1 & 0 \\ 0 & \sqrt{3} \end{bmatrix}.$$
One may decide to work with the diagonal system
[\[eq:diagonalEigen\]](#eq:diagonalEigen){reference-type="eqref"
reference="eq:diagonalEigen"}, and then switch back to the "physical"
coordinates $\bf x$, using ${\bf x} = {\bf P}\,{\bf u}$. Here, it is
immediate to verify that
${\bf P} \, {\bf P}^{\intercal} = {\bf P}^\intercal \, {\bf P} = {\bf I}$,
where ${\bf I}$ is the identity matrix. It is also immediate to check
that [\[eq:eigendempos\]](#eq:eigendempos){reference-type="eqref"
reference="eq:eigendempos"} is verified.

### Loss and Forcing

System [\[eq:TwoMass1\]](#eq:TwoMass1){reference-type="eqref"
reference="eq:TwoMass1"} may be generalised so to include losses and
external forcing. The system reads

::: subequations
[\[eq:TwoMass2\]]{#eq:TwoMass2 label="eq:TwoMass2"} $$\begin{aligned}
m_1 \frac{d^2 x_1}{dt^2} &= -K_{11}x_1 - K_{12}(x_1 - x_2) -2m_1c_{11}\frac{dx_{1}}{dt} + m_1 F_1 f(t), \label{eq:TwoMass2a} \\
m_2 \frac{d^2 x_2}{dt^2} &= -K_{22}x_2 + K_{12}(x_1 - x_2) - 2m_2 c_{22} + m_2 F_2 f(t) . \label{eq:TwoMass2b}\end{aligned}$$
:::

Here, the $c$'s coefficients are loss coefficients (measured in
s$^{-1}$), and $F$ is a force per unit mass. (These units are consistent
with the case of the single mass, in
[\[eq:SHOForced\]](#eq:SHOForced){reference-type="eqref"
reference="eq:SHOForced"}). The system may be written as
$$\label{eq:tempTwoModes1}
\frac{d^2{\bf x}}{dt^2} = -{\bf M}^{-1}{\bf K} \, {\bf x} - 2 {\bf C} \, \frac{d{\bf x}}{dt} + {\bf F}\, f(t).$$
Here, $f(t)$ is any suitable input time signal. Energy analysis may be
performed by left-multiplying the system by
$(\delta_{t\cdot}{\bf x})^\intercal {\bf M}$, leading to
$$\frac{1}{2}\frac{d}{dt}\left(\frac{d{\bf x}^\intercal}{dt} {\bf M} \frac{d{\bf x}}{dt} + {\bf x}^\intercal {\bf K} {\bf x}\right) = -Q(t) + P(t),$$
where
$Q = 2(\delta_{t\cdot}{\bf x})^\intercal \,{\bf M}{\bf C} \, \delta_{t\cdot}{\bf x} \geq 0$
is the dissipated power, and where
$P = (\delta_{t\cdot}{\bf x})^\intercal \, {\bf M}{\bf F} \, f(t)$ is
the injected power.

Conveniently, we may want to study the case where the input is a
time-harmonic signal, such as $f(t) = e^{j\omega t}$. Motion will
undergo an initial transient, before settling into the steady state.
There, the system will oscillate at the same frequency as the forcing.
Hence, in the steady-state, one has
$${\bf x} = {\bf X}\, e^{j\omega t},$$ and inserting this expression in
[\[eq:tempTwoModes1\]](#eq:tempTwoModes1){reference-type="eqref"
reference="eq:tempTwoModes1"} one gets
$${\bf X} = \left( -\omega^2 {\bf I} + {\bf M}^{-1}{\bf K} + 2 j \omega {\bf C} \right)^{-1}{\bf F},$$
that is the generalisation to the vector case of
[\[eq:TransXF\]](#eq:TransXF){reference-type="eqref"
reference="eq:TransXF"}. Here, one may compute the transfer functions
between any of the two inputs and ouputs. It may be convenient to define
the transfer functions are matrices, so that, say, the mechanical
admittance is
$${\bf Y} = \begin{bmatrix} Y_{11} & Y_{12} \\ Y_{21} & Y_{22} \end{bmatrix},$$
where $Y_{ij} = \frac{j \omega X_i}{ F_j}$. Analogously, the impedance
is defined as $Zij = \frac{F_i}{j \omega X_j}$. Fig
[1.4](#fig:TransFunctionsTwoMass){reference-type="ref"
reference="fig:TransFunctionsTwoMass"} presents the admittance and
impedance plots for the test case $K_{11},K_{22},K_{12},m_1,m_2=1$.

![Transfer functions for the two-mass system, with
$K_{11},K_{22},K_{12},m_1,m_2=1$. Here, the loss matrix is
${\bf C} = \text{diag}([0.02,0.01])$, and the forcing vector is
${\bf F} = [1,0]^\intercal$. For the admittance $Y$, dots (.) is
$Y_{11}$ and hats (\^) is $Y_{21}$. For the impedance, dots (.) is
$Z_{11}$, and cicles (o) is $Z_{12}$.
](Figures/TwoMassTransfFunctions.png){#fig:TransFunctionsTwoMass
width="0.85\\linewidth"}

### An Explicit Finite Difference Scheme

A discrete-time version of
[\[eq:MultiModalMatr\]](#eq:MultiModalMatr){reference-type="eqref"
reference="eq:MultiModalMatr"} can be obtained considering

::: subequations
[\[eq:TwoMassFD1\]]{#eq:TwoMassFD1 label="eq:TwoMassFD1"}
$$\begin{aligned}
m_1 \, \delta_{tt}x_1^n = -K_{11}x_1^n - K_{12}(x_1^n - x_2^n), \label{eq:CouplFD1} \\
m_2 \, \delta_{tt}x_2^n = -K_{22}x_2^n + K_{12}(x_1^n - x_2^n), \label{eq:CouplFD2}\end{aligned}$$
:::

which may be written compactly as $$\label{eq:TwoMassFD1Compact}
{\bf M}\, \delta_{tt}{\bf x}^n = - {\bf K}\, {\bf x}^n$$ Expanding out
the difference operators, one gets $$\label{eq:updateExplicitScheme}
{\bf M}\, {\bf x}^{n+1} = \left(2{\bf M}-k^2{\bf K}\right){\bf x}^n - {\bf M}{\bf x}^{n-1}.$$
Since ${\bf M}$ is fully diagonal, this scheme is *explicit*, and the
update may be computed by merely multiplying both sides by the diagonal
matrix ${\bf M}^{-1}$, and by performing the trivial matrix-vector
operations on the right-hand side.

Stability analysis may be performed via energy arguments. To that end,
multiply [\[eq:CouplFD1\]](#eq:CouplFD1){reference-type="eqref"
reference="eq:CouplFD1"} by $\delta_{t\cdot}x^n_1$, and
[\[eq:CouplFD2\]](#eq:CouplFD2){reference-type="eqref"
reference="eq:CouplFD2"} by $\delta_{t\cdot}x^n_2$, and sum. Using the
usual identities, one gets $$\label{eq:EnBalFD2Masses}
\delta_{t+}\left(\underbrace{\frac{m_1(\delta_{t-}x^n_1)^2}{2} + \frac{m_2(\delta_{t-}x^n_2)^2}{2} + \frac{K_{11}x_1^n \, e_{t-}x_1^n}{2}  + \frac{K_{22}x_2^n \, e_{t-}x_2^n}{2} + \frac{K_{12}(x_1^n-x_2^n) \, e_{t-}(x_1^n-x_2^n)}{2}}_{{\mathfrak h}^{n-1/2}}\right)  = 0,$$
which gives the discrete energy balance. Clearly,
[\[eq:EnBalFD2Masses\]](#eq:EnBalFD2Masses){reference-type="eqref"
reference="eq:EnBalFD2Masses"} discretises
[\[eq:EnBalCnt2Masses\]](#eq:EnBalCnt2Masses){reference-type="eqref"
reference="eq:EnBalCnt2Masses"} to second-order accuracy, however there
is no guarantee that the discrete energy is in fact positive. It may be
useful to rewrite the energy as a quadratic form using the matrix
notation. One has
$$\mathfrak{h}^{n-1/2} = \frac{1}{2}\left({\delta_{t-}{\bf x}}^{\intercal}\,{\bf M}\, {\delta_{t-}{\bf x}} + e_{t-}{\bf x}^\intercal \,{\bf K}\, { {\bf x}}\right) = \frac{1}{2}\left( \delta_{t-}{\bf x}^{\intercal}\,\left({\bf M}-\frac{k^2}{4}{\bf K}\right)\, \delta_{t-}{\bf x} + \mu_{t-}\left({\bf x}^{\intercal}\,{\bf K}\, {\bf x}\right) \right),$$
where the last equality follows after application of the identity
$e_{t-}{\bf u}^\intercal  {\bf A} {\bf u}= \mu_{t-}\left({\bf u}^\intercal \,{\bf A}\, {\bf u}\right) - \frac{k^2}{4}\delta_{t-}{\bf u}^\intercal \,{\bf A}\, \delta_{t-}{\bf u}$,
for a positive-definite, symmetric matrix $\bf A$. Now, since
$\mu_{t-}\left({\bf x}^{\intercal}\,{\bf K}\,  {\bf x} \right) \geq 0$,
one has that the discrete energy above is positive definite if and only
if $$\text{eig}\left({\bf M}-\frac{k^2}{4}{\bf K}\right) \geq 0.$$
Finding conditions in the general case can be quite tedious. For the
case $K_{11},K_{22},K_{12},m_1,m_2=1$, the eigenvalue equation for
eigenvalue $\lambda$ is obtained as
$$\left(\frac{k^2}{2} + \lambda - 1\right)^2 - \frac{k^2}{16}  \geq 0.$$
The quadratic has two solution, $\lambda_+ = 1 - \frac{k^2}{4}$,
$\lambda_- = 1 - \frac{3 k^2}{4}$, and they are both positive if and
only if $\lambda_+ > 0$, i.e. $$\label{eq:TwoMassesStabCond}
k \leq \frac{2}{\sqrt{3}},$$ which may be interpreted as a stability
condition here.

### Numerical Eigenvalues

Frequency-domain analysis may be performed here too, generalising the
results of Sec.
[\[sec:FreqDomSHO\]](#sec:FreqDomSHO){reference-type="ref"
reference="sec:FreqDomSHO"}. To that end, one re-writes the scheme
[\[eq:TwoMassFD1\]](#eq:TwoMassFD1){reference-type="eqref"
reference="eq:TwoMassFD1"} as
$$\delta_{tt}{\bf x}^n = - {\bf M}^{-1}{\bf K}\, {\bf x}^n.$$ Then, as
suggested in Sec.
[\[sec:FDtransformations\]](#sec:FDtransformations){reference-type="ref"
reference="sec:FDtransformations"}, one substitutes the test solution
${\bf x}^n = \hat{\bf x}e^{j\omega kn}$, and remembering
[\[eq:DTFTdtt\]](#eq:DTFTdtt){reference-type="eqref"
reference="eq:DTFTdtt"} one gets
$$-\frac{4}{k^2}\sin^2\left( \frac{\omega k}{2}\right)\hat{\bf x} = - {\bf M}^{-1}{\bf K}\, \hat {\bf x}^n,$$
showing that the eigenvalues of the matrix ${\bf M}^{-1}{\bf K}$ are
$$\frac{4}{k^2}\sin^2\left( \frac{\omega k}{2}\right) = \omega^2 - \frac{k^2\omega^4}{12} + O(k^4).$$
Hence, the numerical eigenvalues are second-order accurate with respect
to the eigenvalues of the continuos system, and the absolute value of
error grows approximately as $\omega^4$. Thus, the eigenfrequencies will
be less and less well computed as the stability limit is approached.
Note that one may decompose the right-hand using the eigendecomposition
[\[eq:eigendempos\]](#eq:eigendempos){reference-type="eqref"
reference="eq:eigendempos"}, and, defining again
${\bf u} = {\bf P}^{\intercal}{\bf x}$, the system becomes
$$\delta_{tt}{\bf u}^n = -{\bf \Omega}^2 \,{\bf u}^n,$$ showing that
motion can be completely uncoupled in the two modes of the system.
Stability analysis here is immediate, since the energy is now expressed
as the sum of *independent* harmonic oscillators. Hence, stability
condition [\[eq:StabCondSHO\]](#eq:StabCondSHO){reference-type="eqref"
reference="eq:StabCondSHO"} translates here directly, as
$$k < \frac{2}{\text{max}\left(\text{diag}({\bf \Omega})\right)}.$$ For
the test case $K_{11},K_{22},K_{12},m_1,m_2=1$, one recovers of course
[\[eq:TwoMassesStabCond\]](#eq:TwoMassesStabCond){reference-type="eqref"
reference="eq:TwoMassesStabCond"}.

### A Family of Finite Difference Schemes

Scheme
[\[eq:TwoMassFD1Compact\]](#eq:TwoMassFD1Compact){reference-type="eqref"
reference="eq:TwoMassFD1Compact"} is but one among many different
possible realisations. Consider the following discretisation, depending
on a parameter $\alpha \in [0,1]$. $$\label{eq:Impl2Masses}
{\bf M}\, \delta_{tt}{\bf x}^n = - {\bf K}\, (\alpha + (1-\alpha)\mu_{t\cdot}) {\bf x}^n.$$
Remebering that $\mu_{t\cdot}= 1 + \frac{k^2}{2}\delta_{tt}$, one can
rewrite the scheme as
$$\left({\bf M} + (1-\alpha)\frac{k^2}{2}{\bf K} \right) \, {\bf x}^{n+1} = \left(2 {\bf M} - \alpha k^2 {\bf K} \right){\bf x}^n - \left({\bf M} + (1-\alpha)\frac{k^2}{2}{\bf K} \right) \, {\bf x}^{n-1}.$$
Here, the scheme is *implicit*, since the matrix multiplying the update
vector ${\bf x}^{n+1}$ is now (generally) not diagonal! The update
equation is in the form of a linear system. This is, of course, more
laborious than the update
[\[eq:updateExplicitScheme\]](#eq:updateExplicitScheme){reference-type="eqref"
reference="eq:updateExplicitScheme"} for the fully-explicit scheme. Of
course, setting $\alpha = 1$ in the above, one recovers the explicit
scheme
[\[eq:TwoMassFD1Compact\]](#eq:TwoMassFD1Compact){reference-type="eqref"
reference="eq:TwoMassFD1Compact"}. Energy analysis is performed in the
usual way, i.e. by left-multiplying
[\[eq:Impl2Masses\]](#eq:Impl2Masses){reference-type="eqref"
reference="eq:Impl2Masses"} by $(\delta_{t\cdot}{\bf x}^n)^\intercal$.
By means of the usual identities, this leads to the discrete energy
balance
$$\delta_{t+}\frac{1}{2}\left( \delta_{t-}{\bf x}^{\intercal}\,{\bf M}\, \delta_{t-}{\bf x} + \alpha e_{t-}{\bf x}^{\intercal}\,{\bf K}\,  {\bf x} + (1-\alpha) \mu_{t-}\left({\bf x}^{\intercal}\,{\bf K}\,  {\bf x}\right) \right) \triangleq  \delta_{t+}{\mathfrak h}^{n-1/2}  = 0.$$
An interesting case is obtained when $\alpha =0$, yielding a form for
the discrete energy that is non-negative *in all cases*. Hence, scheme
[\[eq:Impl2Masses\]](#eq:Impl2Masses){reference-type="eqref"
reference="eq:Impl2Masses"} becomes *unconditionally stable* when
$\alpha = 0$. One should be aware that unconditional stability is but
one of the many aspects to be considered whilst designing an appropriate
scheme. Whilst scheme
[\[eq:Impl2Masses\]](#eq:Impl2Masses){reference-type="eqref"
reference="eq:Impl2Masses"} will remain stable, regardless of the sample
rate used, multiple other issues may affect the quality of the resulting
simulated dynamics, particularly with respect to Fourier and aliasing
considerations. Fig. [1.5](#fig:TwoMassOne){reference-type="ref"
reference="fig:TwoMassOne"} presents a numerical investigation of the
two mass system, for the test case $K_{11},K_{22},K_{12},m_1,m_2=1$.

A stability condition for scheme
[\[eq:Impl2Masses\]](#eq:Impl2Masses){reference-type="eqref"
reference="eq:Impl2Masses"} may be obtained by writing out the
expression for the energy as a quadratic form in
${\bf p} = [\left({\bf x}^n\right)^{\intercal}, \left({\bf x}^{n-1}\right)^{\intercal}]^{\intercal}$,
as
$\mathfrak{h}^{n-1/2} = {\bf p}^{\intercal}\,{\bf A}(\alpha,k)\,{\bf p}$,
where $${\bf A}(\alpha,k) = 
\begin{bmatrix}
\frac{{\bf M}}{k^2} + \frac{(1-\alpha)}{2}{\bf K} & -\frac{{\bf M}}{k^2} +\frac{\alpha}{2}{\bf K} \\ -\frac{{\bf M}}{k^2} + \frac{\alpha}{2}{\bf K} & \frac{{\bf M}}{k^2} + \frac{(1-\alpha)}{2}{\bf K}
\end{bmatrix}.$$ Positive-definiteness is obtained if and only if
$$\label{eq:StabCondMasses}
\text{eig}\left({\bf A}(\alpha,k)\right) > 0.$$ This can be interpreted
a stability condition for the system.

![Free oscillations of the two-mass system. Here,
$K_{11},K_{22},K_{12},m_1,m_2=1$. The sample rate is chosen as
$f_s = 50$ Hz, and $\alpha = 0.5$. The initial conditions are given as
$x_1(0) = 1$, $\frac{dx_1(0)}{dt} =0$, $x_2(0) = 0$,
$\frac{dx_2(0)}{dt} =0$. (a): time evolution of $m_1$ (dots) and $m_2$
(circles). (b): spectra of the solutions. (c): energy components:
kinetic (\^), potential (), and total (solid line). (d): numerical
energy error, defined as
$\Delta H = 1 - \mathfrak{h}^{n-1/2}/\mathfrak{h}^{1/2}$.](Figures/TwoMassOne.png){#fig:TwoMassOne
width="0.85\\linewidth"}

#### Loss and Forcing

Inclusion of losses and external forcing is immediate from the template
given above. Considering the continuous system
[\[eq:tempTwoModes1\]](#eq:tempTwoModes1){reference-type="eqref"
reference="eq:tempTwoModes1"}, a discrete-time version is obtained as
$$\delta_{tt}{\bf x}^n = -{\bf M}^{-1}{\bf K}\, (\alpha + (1-\alpha)\mu_{t\cdot}) {\bf x}^n - 2{\bf C}\delta_{t\cdot}{\bf x}^n + {\bf F} f^n.$$
The energy balance is obtained after left-multiplying the equation by
$(\delta_{t\cdot}{\bf x}^n)^\intercal{\bf M}$, yielding
$$\frac{1}{2}\delta_{t+}\left( \delta_{t-}{\bf x}^{\intercal}\,{\bf M}\, \delta_{t-}{\bf x} + \alpha e_{t-}{\bf x}^{\intercal}\,{\bf K}\,  {\bf x} + (1-\alpha)\mu_{t-}\left({\bf x}^{\intercal}\,{\bf K}\,  {\bf x}\right)  \right)   = -Q^n+P^n,$$
where
$Q^n = 2(\delta_{t\cdot}{\bf x}^n)^\intercal \, {\bf C} \, \delta_{t\cdot}{\bf x}^n \geq 0$
is the dissipated power, and where
$P^n = ({\delta_{t\cdot}{\bf x}^n})^\intercal \, {\bf F} \, f^n$ is the
injected power. Here, since ${\bf C}$ is a symmetric, positive-definite
matrix, the stability analysis is unaffected, and one may still impose
[\[eq:StabCondMasses\]](#eq:StabCondMasses){reference-type="eqref"
reference="eq:StabCondMasses"}. Of course, the scheme above allows to
study the system under general forcing conditions, and to simulate the
evolution of the system including transients.

[^1]: Remember that a positive-definite $N\times N$ symmetric matrix
    ${\bf A}$ is a matrix with real entries such that
    ${\bf x}^\intercal {\bf A} {\bf x} > 0$,
    $\forall {\bf x} \neq {\bf 0} \in \mathbb{R}^N$.
