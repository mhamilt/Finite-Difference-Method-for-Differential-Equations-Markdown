---
layout: post
title: Lecture 2
---

-   [The Simple Harmonic Oscillator](#chap:SHO)
    -   [A finite difference scheme](#a-finite-difference-scheme)
        -   [Stability via frequency domain analysis](#sec:FreqDomSHO)
        -   [Frequency warping and modified equation
            techniques](#eq:ModEqTechniques)
    -   [Loss](#sec:LossSHO)
        -   [Energy analysis](#energy-analysis)
        -   [A finite difference scheme](#a-finite-difference-scheme-1)
        -   [Higher-order schemes](#higher-order-schemes)
    -   [Forced Oscillations](#forced-oscillations)
        -   [Response to harmonic input
            forcing](#response-to-harmonic-input-forcing)


# The Simple Harmonic Oscillator {#chap:SHO}

Simple harmonic motion is obtained under linear condition, that is
$\$\n\begin{equation}\phi = \frac{Kx^2}{2}
\end{equation}\$\$ 

Under such choice,
[\[eq:SHO\]](#eq:SHO){reference-type="eqref" reference="eq:SHO"} becomes
$$\begin{equation}\label{eq:SHOreal}
    \frac{d^2 x}{dt^2} = - \omega_0^2 x,\end{equation}$$ 

where the radian natural
frequency $\omega_0 = \sqrt{K/m}$ was introduced. There are multiple
reasons to be wanting to study the simple harmonic oscillator here: a
variety of mechanical systems can be approximated by
[\[eq:SHOreal\]](#eq:SHOreal){reference-type="eqref"
reference="eq:SHOreal"}, see e.g. Fig.
[1.1](#fig:oscExamples){reference-type="ref"
reference="fig:oscExamples"}. Furthermore, as will be seen in later
lectures, distributed systems may themselves be approximated as a bank
of parallel oscillators, each one corresponding to one "mode" of
vibration. Finally, the simple harmonic oscillator possesses exact (i.e.
closed-form) solutions, both analytically and numerically, and is thus a
very useful test case for the numerical schemes that will be introduced
below. It was seen that the solution is expressed as the sum of two
periodic functions, i.e.
$$\begin{equation}
x(t) = a \cos(\omega_0 t) + b \sin{\omega_0 t}.
\end{equation}
$$ 

The constants $a,b$
are uniquely determined from the intial conditions. Setting $x(t=0)=x_0$
and $\frac{dx}{dt}(x=0)=v_0$, one gets 

$$
\begin{equation}\label{eq:SHOexact}
    x(t) = x_0 \cos(\omega_0 t) + \frac{v_0}{\omega_0}\sin(\omega_0 t)
\end{equation}
$$
The energy components are obtained explictly as
$$
E_k = \frac{m}{2}\left(x_0\omega_0\sin(\omega_0 t) - v_0 \cos(\omega_0 t) \right)^2, \quad E_p = \frac{K}{2}\left( x_0 \cos(\omega_0 t) + \frac{v_0}{\omega_0}\sin(\omega_0 t) \right)^2.
$$
Summing the two expressions together, one gets
$$\begin{equation}
H(t) = \frac{m v_0^2}{2} + \frac{K x_0^2}{2} = H_0,
\end{equation}
$$ 

i.e.
[\[eq:EnCons\]](#eq:EnCons){reference-type="eqref"
reference="eq:EnCons"}.

![Examples of harmonic oscillators. The mass-spring system (a), and the
series RLC circuit (b) are usually given as examples of harmonic
oscillators. The case of a cantilever beam with a point mass (c), and
the pendulum (d) may be approximated as harmonic oscillators in the case
of small vibrations.](Figures/OscillatorsExamples.pdf){#fig:oscExamples
width="\\linewidth"}

## A finite difference scheme

Consider the time series $x^n$, approximating the true solution $x(t)$
of [\[eq:SHOreal\]](#eq:SHOreal){reference-type="eqref"
reference="eq:SHOreal"}. As a first example of a working finite
difference scheme, consider $$\begin{equation}\label{eq:Scheme1}
    \delta_{tt}x^n = -\omega_0^2 x^n.\end{equation}$$ 

Expanding out the operator, one
gets $$\begin{equation}\label{eq:FDbasic}
    x^{n+1} = x^n(2-\omega_0^2 k^2) - x^{n-1},\end{equation}$$ 

and hence the update
requires one multiply and one sum. Clearly, it is convenient to store
the value $2-\omega_0^2 k^2$ offline, so that one need not recompute
this value at each time step. The update
[\[eq:FDbasic\]](#eq:FDbasic){reference-type="eqref"
reference="eq:FDbasic"} is *explicit*: by this term, we denote a scheme
such that 
$$\begin{equation}
x^{n+1} = f(x^n,x^{n-1}), \qquad \text{(explicit scheme)}
\end{equation}
$$
We will soon encounter schemes that are instead *implicit*, i.e.
$$\begin{equation}\label{eq:ImplSchemeDef}
    g(x^{n+1})= f(x^n,x^{n-1}), \qquad \text{(implicit scheme)}\end{equation}$$ 

where
$g$ is generally a nonlinear function in $x^{n+1}$.

Though scheme [\[eq:Scheme1\]](#eq:Scheme1){reference-type="eqref"
reference="eq:Scheme1"} looks reasonable, there is no guarantee (for the
moment) that the computed solutions are indeed an approximate form of
the true solution $x(t)$. In some cases, as will be seen shortly, the
time series computed by $\eqref{eq:Scheme1}$ diverges, in some other
cases, it remains bounded. The next few sections will explain the idea
of *convergence*, and the closely linked idea of *stability*. A notion
of stability may be formalised as follows. We say that a scheme is
stable in stability region $\Lambda \subseteq \mathbb{R}^+$ if, for any
positive time constant $\tau \leq n k$, there exist a positive index $M$
such that $$\begin{equation}\label{eq:StabDef}
    |x^n| \leq C_\tau \sum_{m = 0}^M |x^m|\end{equation}$$ 

for a constant $C_\tau > 0$
independent of $k \in \Lambda$. In practice, stability is defined as a
bound on the time series including the first $M+1$ steps.

### Stability via frequency domain analysis {#sec:FreqDomSHO}

As anticipated earlier, frequency-domain techinques may be employed in
the analysis of stability of linear, time-invariant discrete systems,
such as [\[eq:Scheme1\]](#eq:Scheme1){reference-type="eqref"
reference="eq:Scheme1"}. For that, an *ansatz* of the form
[\[eq:ansatz\]](#eq:ansatz){reference-type="eqref"
reference="eq:ansatz"} is substituted: 
$$\begin{equation}
x^n = \hat x z^{n}.
\end{equation}
$$ 

In this
equation, notation is a little mixed-up, since on the left-hand side $n$
denotes the time index, whereas on the right-hand side it is an
exponent! In practice, we keep the same apex notation in both cases, for
the sake of notation, but the meaning is very different. We get

$$\begin{equation}
\hat x z^n\left( z - (2-\omega_0 k^2) + z^{-1}\right) = 0, \quad \rightarrow \quad z_{\pm} = \frac{2 - \omega_0^2 k^2 \pm \omega_0^2 k^2\sqrt{1 - \frac{4}{\omega_0^2 k^2}}}{2}.
\end{equation}
$$
Hence, the solution is $$\begin{equation}\label{eq:sol_z_SHO}
    x^n =  A_+ z^n_+ + A_- z^n_-,\end{equation}$$ 

for complex constants $A_\pm$. We
assume that the scheme is started using two starting values $x^0, x^1$
(obtained from $x_0$, $v_0$ of the continuous problem). Then

$$\begin{equation}
x^0 = A_+ + A_-, \quad x^1 = A_+ z_+ + A_- z_-.
\end{equation}
$$ 

From these, the
complex constants are obtained as 

$$\begin{equation}\label{eq:ApAm}
    A_+ = \frac{x^0z_- - x^1}{z_- - z_+},\,\,\, A_- = \frac{x^1 - x^0z_+ }{z_- - z_+}.
\end{equation}
$$

If the square root in $z_\pm$ is a real number, than $z_-$ has magnitude
larger than unity, and the solution $x^n$ will therefore grow
exponentially over time: this is *instability*. On the other hand, when
the square root in imaginary, then $z_\pm$ become oscillating, and
$z_\pm$ are complex conjugates. This condition is obtained as

$$\begin{equation}\label{eq:StabCondSHO}
    k < \frac{2}{\omega_0},\end{equation}$$ 

which is an upper bound on the time step,
once the natural frequency of the oscillator is set. In the case of
oscillating solutions, one has 
$$\begin{equation}
z_\pm = r e^{\pm j \theta},
\end{equation}
$$ 

with
$r = 1$, and
$\tan\theta = \left(\omega_0^2k^2\sqrt{\frac{4}{\omega_0^2 k^2}-1}\right)/(2-\omega_0^2 k^2)$.

To check stability, definition
[\[eq:StabDef\]](#eq:StabDef){reference-type="eqref"
reference="eq:StabDef"} is applied, to give 

$$\begin{equation}\label{eq:SHObound} 
    |x^n| = |A_+ r^n e^{j \theta n} + A_- r^n e^{-j \theta n}| \leq |A_+| + |A_-| \leq \frac{|x^0|+|x^1|}{|\sin \theta|} \triangleq C_\tau \sum_{m=0}^1 |x^m|,
    \end{equation}$$

and thus the absolute value of the solution at the time $n>1$ is bounded
in terms of the values at $n=0,1$, with bounding constant
$C_\tau = 1/|\sin\theta|$. (The first inequality in the above was
obtained via the triangle inequality. Then, the fact that
$|r|=|e^{\pm j\theta}|=1$ was used, and finally the values from
[\[eq:ApAm\]](#eq:ApAm){reference-type="eqref" reference="eq:ApAm"} were
substituted in). A numerical check of the current bound is given in Fig.
[1.2](#fig:SHObounds){reference-type="ref" reference="fig:SHObounds"}.

![Simple Harmonic Oscillator. Numerical check on bound
[\[eq:SHObound\]](#eq:SHObound){reference-type="eqref"
reference="eq:SHObound"}. The sample rate used in the examples is
$f_s = 2000$ Hz. Starting values are $x^0 = x^1 = 1$. Dashed line line
is bound [\[eq:SHObound\]](#eq:SHObound){reference-type="eqref"
reference="eq:SHObound"}, characteristic frequency as indicated, with
$\omega_{max}=2/k$.](Figures/boundsSHO.png){#fig:SHObounds
width="\\linewidth"}

Note that the same condition may be arrived at via *von Neumann*
analysis. Recall the DTFT of the $\delta_{tt}$ operator, as per
[\[eq:DTFTdtt\]](#eq:DTFTdtt){reference-type="eqref"
reference="eq:DTFTdtt"}. Transforming
[\[eq:Scheme1\]](#eq:Scheme1){reference-type="eqref"
reference="eq:Scheme1"} in the frequency domain accordingly, one gets

$$\begin{equation}
\left(-\frac{4}{k^2}\sin^2\left( \frac{\omega k}{2}\right) + \omega_0^2\right)\mathcal{X}\left\{ x^n\right\} = 0.
\end{equation}
$$
One must impose $0 \leq \sin^2\left( \frac{\omega k}{2}\right) \leq 1$,
which is possible if and only if
[\[eq:StabCondSHO\]](#eq:StabCondSHO){reference-type="eqref"
reference="eq:StabCondSHO"} holds. If this condition is violated, one
has that the frequency $\omega$ becomes pure imaginary, thus the sine
becomes a hyerbolic sine, with unbounded growth (instability).

#### Stability via energy analysis

The discussion of Sec.
[\[sec:EnAnGen\]](#sec:EnAnGen){reference-type="ref"
reference="sec:EnAnGen"} suggests that, if the model problem can be
shown to have a conserved total energy, with non-negative kinetic and
potential terms, then the solution can be bounded in terms of the
energy. It may be tempting to try to find a discrete version of
[\[eq:EnAnGen\]](#eq:EnAnGen){reference-type="eqref"
reference="eq:EnAnGen"}, valid in the discrete case. To that end,
[\[eq:Scheme1\]](#eq:Scheme1){reference-type="eqref"
reference="eq:Scheme1"} is multiplied by $m \delta_{t\cdot}x^n$, to get
$$\begin{equation}\label{eq:SHOen1}
    m \delta_{t\cdot}x^n (\delta_{tt}x^n +\omega_0^2 x^n) = 0\end{equation}$$ 

A couple of useful identities (that will be used throughout) are given here:
$$\begin{equation}\label{eq:IdsFD}
    \delta_{t\cdot}x^n \, \delta_{tt}x^n = \delta_{t+}\left( \frac{(\delta_{t-}x^n)^2}{2}\right), \,\, \delta_{t\cdot}x^n \, x^n = \delta_{t+}\left( \frac{x^n e_{t-}x^n}{2}\right).\end{equation}$$

These are proven by simple algebra. Using these identities in
[\[eq:SHOen1\]](#eq:SHOen1){reference-type="eqref"
reference="eq:SHOen1"}, one gets


$$\begin{equation}
\delta_{t+}\left( \frac{m}{2}{(\delta_{t-}x^n)^2} + \frac{K}{2} {x^n e_{t-}x^n} \right) = 0,
\end{equation}
$$

that is a discrete counterpart of
[\[eq:En2\]](#eq:En2){reference-type="eqref" reference="eq:En2"}.
Multiplication by $m$ was here used so to yield units of energy in the
expression within the brakets, though of course one may obtain conserved
energy per unit mass via multiplication by $\delta_{t\cdot}x^n$ alone.
It is easy, in the above, to recognise a discrete approximation to the
continuous energy balance [\[eq:En1\]](#eq:En1){reference-type="eqref"
reference="eq:En1"}. In this case, one may define an *interleaved* time
series $\mathfrak{h}^{n-1/2}$, corresponding to the discrete conserved
energy: 
$$\begin{equation}\label{eq:DiscEnSHO}
    \mathfrak{h}^{n-1/2} \triangleq \frac{m}{2}{(\delta_{t-}x^n)^2} + \frac{K}{2} {x^n e_{t-}x^n} = \mathfrak{h}^{1/2}.\end{equation}$$
In light of the discussion in the continuous case, one may of course use
energy conservation as a means to bound the growth of solutions over
time. The problem here, is that $\mathfrak{h}^{n-1/2}$ may *not* be
positive, since the potential energy is of indefinite sign. Instances
leading to negative energy overall are a manifestation of instability,
and must be avoided. It may be useful, then, to bound the potential term
in the energy expression. Using

$$\begin{equation}
x^n e_{t-}x^n  =  (\mu_{t-}x^n)^2 - \frac{k^2}{4}\left(\delta_{t-}x^n \right)^2,
\end{equation}
$$
the total energy is 
$$\begin{equation}\label{eq:ModEnSHO}
    \mathfrak{h}^{n-1/2} = \frac{m\left(\delta_{t-}x^n \right)^2}{2} \left(1 - \frac{\omega_0^2 k^2}{4}\right) + \frac{K (\mu_{t+}x^n)^2}{2} \geq \frac{m\left(\delta_{t-}x^n \right)^2}{2} \left(1 - \frac{\omega_0^2 k^2}{4}\right)\end{equation}$$
and thus the total energy will be non-negative if and only if
[\[eq:StabCondSHO\]](#eq:StabCondSHO){reference-type="eqref"
reference="eq:StabCondSHO"} is satisfied. Fig.
[1.3](#fig:EnConsSHO){reference-type="ref" reference="fig:EnConsSHO"}
shows the energy components and error for scheme
[\[eq:Scheme1\]](#eq:Scheme1){reference-type="eqref"
reference="eq:Scheme1"}. The energy components are given as per
[\[eq:DiscEnSHO\]](#eq:DiscEnSHO){reference-type="eqref"
reference="eq:DiscEnSHO"}, and it is remarked that the potential term
*does* become negative at times. It is the *overall* energy that is
positive. Of course, the equivalent expression in
[\[eq:ModEnSHO\]](#eq:ModEnSHO){reference-type="eqref"
reference="eq:ModEnSHO"} has modified expressions for the kinetic and
potential energies, that are always non-negative under stability
condition [\[eq:StabCondSHO\]](#eq:StabCondSHO){reference-type="eqref"
reference="eq:StabCondSHO"}.

![Energy behaviour of scheme
[\[eq:Scheme1\]](#eq:Scheme1){reference-type="eqref"
reference="eq:Scheme1"}. Left: energy components. Kinetic (dashed),
potential (dash-dotted), total (solid). Right: energy error
$\mathfrak{h}^{n-1/2}/\mathfrak{h}^{1/2}-1$. The oscillator is
initialised with $x_0=v_0=1$, and $\omega_0=100$
rad/s.](Figures/EnergyErr.png){#fig:EnConsSHO width="\\linewidth"}

#### Consistency and accuracy

We are now introducing the idea of *local truncation error* (LTE),
denoted here $\varepsilon^n$. Applying the finite difference scheme to
the true solution $x(t)$ yields a definiton of the LTE as
$$\begin{equation}\label{eq:LTEdef}
    \delta_{tt}x(t_n) + \omega_0^2 x(t_n) = \varepsilon^n.\end{equation}$$ 

Using
Taylor series arguments, as per
[\[eq:Errs\]](#eq:Errs){reference-type="eqref" reference="eq:Errs"}, one
gets

$$\begin{equation}
\left( \frac{d^2 x(t)}{dt} + \omega_0^2 x(t)\right)|_{t=t_n} + O(k^2) = \varepsilon^n,
\end{equation}
$$
and since $x(t)$ is the true solution, one recovers
$\varepsilon^n = O(k^2)$. The behaviour of the LTE as a function of $k$
describes the idea of *consistency*: a scheme is said to be consistent
if 
$$\begin{equation}
\lim_{k\rightarrow 0}\varepsilon^n = 0.
\end{equation}
$$ 

In practice, consistent
schemes are such that the local error becomes small as $k$ is decreased.
Usually, $\varepsilon = O(k^p)$, and one may conclude that the scheme is
$p^{th}$-order accurate. Of course, this is not entirely true, since the
question of accuracy is tightly bound to the ideas of stability and
convergence: higher-accurate schemes may *never* converge for a given
model problem. The idea of accuracy is only going to be meaningful when
a scheme is provably stable in some manner. As an example, consider a
fourth-order accurate difference operator discretising the second time
derivative: $$\begin{equation}\label{eq:FourthOrderdtt}
    \bar\delta_{tt} x(t) = \left(\frac{-e_{t+}^2 + 16 e_{t+}- 30 + 16e_{t+}- e_{t+}^2}{12k^2}\right) x(t) = \frac{d^2 x}{dt^2} + O(k^4).\end{equation}$$
Though technically "higher" accurate, this approximation is always
unstable, even for the simple problem of a free particle ($\phi = 0$ in
[\[eq:PhiF\]](#eq:PhiF){reference-type="eqref" reference="eq:PhiF"}).
Using the test solution $x^n = \hat x z^n$ for this test case, one has

$$\begin{equation}
\hat x z^n \left(-z^2 + 16 z - 30 + 16 z^{-1} - z^{-2}\right) = 0,
\end{equation}
$$
and it is easy to verify that there exist one (real) root
$z \approx 13.9282$ for which clearly $|z|>1$. The scheme is unstable,
and the higher accuracy of error of $\bar\delta_{tt}$ has no real
advantage.

#### Convergence

In turn, what we are really interested in is the evolution of the global
error $E^n$, which must remain bounded. For stable, consistent schemes,
the global error $E^n$ can be expected to maintain the same trend as the
local truncation error $\varepsilon^n$, though this claim would require
a formal proof. As an example, the output of scheme
[\[eq:Scheme1\]](#eq:Scheme1){reference-type="eqref"
reference="eq:Scheme1"} is compared against the exact solution given in
[\[eq:SHOexact\]](#eq:SHOexact){reference-type="eqref"
reference="eq:SHOexact"}. The initial conditions in the continuous
system are set as $x_0 =1$, $v_0 = 0$, giving $x(t) = \cos(\omega_0 t)$.
The initial conditions in the discrete scheme are given as
$x^0 = x_0 = 1$, $x^1 = \cos(\omega_0 k)$ (which discretises the exact
solution at the time $t=k$). Fig.
[1.4](#fig:SHOerrs){reference-type="ref" reference="fig:SHOerrs"} (a
double log plot) shows that indeed slopes of 2 are recovered.

![Global error of scheme
[\[eq:Scheme1\]](#eq:Scheme1){reference-type="eqref"
reference="eq:Scheme1"}. Here, $nk = 1$. Initial conditions are given as
$x_0 = 1$, $v_0 = 0$, giving $x(t) = \cos(\omega_0 t)$. The numerical
initial conditions are $x^0 = 1$, $x^1 = \cos(\omega_0 k)$. The three
lines correspond to $\omega_0 = 100$ rad/s (solid), $\omega_0 = 200$
rad/s (dashed), $\omega_0 = 300$ rad/s (dash-dotted). Sample rate
$f_s = 2000$ Hz.](Figures/OscError.pdf){#fig:SHOerrs
width="\\linewidth"}

The behaviour of the global error as a function of the time step $k$
encapsulates the idea of convergence. A scheme is convergent if
$$\begin{equation}\label{eq:convDef}
    \lim_{k\rightarrow 0} (x(t_n)-x^n) = \lim_{k\rightarrow 0} E^n = 0.\end{equation}$$
In practice, as the time step is decreased, the global error goes to
zero. We will come back to the idea of convergence later on. In general,
a stable and consistent method is also convergent (this should be proven
rigourously, and we will postpone this discussion until later sections).

#### Initialisation {#sec:Init}

In the previous subsection, the scheme was initialised exactly using
knowledge coming from the exact solution $x(t)$ as per
[\[eq:SHOexact\]](#eq:SHOexact){reference-type="eqref"
reference="eq:SHOexact"}. Of course, generally an exact solution is not
available, and schemes must be initialised in some other manner. Since
we usually know the initial position and velocity of the oscillator (in
continuous time) $x_0, v_0$, we would like to know how to use this
information to extract suitable initial values for the time series, i.e.
$x^0$, $x^1$.

![Global error of scheme
[\[eq:Scheme1\]](#eq:Scheme1){reference-type="eqref"
reference="eq:Scheme1"}. Initialisation with first-order accurate
(markers) and second-order accurate (lines) initial conditions. Here,
$nk = 1$. Initial conditions are given as $x_0 = 0.5$, $v_0 = 0.5$. The
numerical initial conditions are $x^0 = x_0$,
$x^1 = kv_0 + x_0 - 0.5 k^2 \omega_0^2 x_0$ (lines), $x^0 = x_0$,
$x^1 = kv_0 + x_0$ (markers). The three cases correspond to
$\omega_0 = 100$ rad/s, $\omega_0 = 200$ rad/s, $\omega_0 = 300$
rad/s.](Figures/OscErrorOrder.pdf){#fig:SHOerrsOrders
width="1\\linewidth"}

Obviously, one can set 
$$\begin{equation}
x^0=x_0
\end{equation}
$$ 

For $x^1$, one possible solution is
to use a simple forward difference to compute it from $v_0$ and $x^0$,
i.e.

$$\begin{equation}
\delta_{t+}x^0 = v_0 \,\,\, \rightarrow \,\,\, x^1 = x^0 + kv_0.
\end{equation}
$$
This approximation to the initial conditions is only *first-order
accurate.* This is easily proven via Taylor series arguments. One may be
tempted then to use the centered difference $\delta_{t\cdot}$, instead
of the forward difference $\delta_{t+}$, since it yields a second-order
accurate approximation to the time derivative. Of course, this is not
possible directly, since the stencil of $\delta_{t\cdot}$ is too large:
we do not know what $x^{-1}$ is (in fact, this value is not defined).
However, consider the following expression for the centered time
difference: $$\begin{equation}\label{eq:dtdInit}
    \delta_{t\cdot}= \delta_{t+}- \frac{k}{2} \delta_{tt}= \frac{d}{dt} + O(k^2)\end{equation}
    $$
Of course, the value of $\delta_{tt}$ is not substituted directly,
rather via [\[eq:Scheme1\]](#eq:Scheme1){reference-type="eqref"
reference="eq:Scheme1"}, i.e. $\delta_{tt}= -\omega_0^2$. Thus, a
second-order accurate intial condition can be given as

$$\begin{equation}
\delta_{t\cdot}x^0 = v_0 \,\,\, \rightarrow \,\,\, x^1 = x^0 + kv_0 - \frac{k^2}{2}\omega_0^2x^0.
\end{equation}
$$
Fig. [1.5](#fig:SHOerrsOrders){reference-type="ref"
reference="fig:SHOerrsOrders"} shows the error plots, computed against
the exact solution
[\[eq:SHOexact\]](#eq:SHOexact){reference-type="eqref"
reference="eq:SHOexact"}, displaying the expected trends in the limit of
high sample rate. Note that, for lower values of the sample rate, the
error of the first-order accurate scheme may in fact be lower than the
second-order accurate scheme, though in the limit of vanishing time step
the correct trends are recovered, and convergence is of course faster
for the second-order accurate schemes.

Higher order accurate approximations are of course possible. Here is a
list, obtained using Taylor-series arguments. $$\begin{aligned}
\label{eq:HigherOrderICsOscillator}
\delta_{t+}x^0 &= v_0 \quad \text{first order}\\
\left(\delta_{t+}-\frac{k}{2}\delta_{tt}\right) x^0 &= v_0 \quad \text{second order}\\
\left(\delta_{t+}-\frac{k}{2}\delta_{tt}- \frac{k^2}{6}\delta_{t+}\delta_{tt}\right) x^0 &= v_0 \quad \text{third order}\\
\left(\delta_{t+}-\frac{k}{2}\delta_{tt}- \frac{k^2}{6}\delta_{t+}\delta_{tt}- \frac{k^3}{24}\delta_{tt}^2 \right) x^0 &= v_0 \quad \text{fourth order}\end{aligned}$$
In the expressions above, substituting $(\delta_{tt})^p = (-\omega_0)^p$
gives a way to compute $x^1$, knowing $x^0$ and $v_0$.

### Frequency warping and modified equation techniques {#eq:ModEqTechniques}

The discussion about accuracy has so far dealt with the idea of
order-accuracy, i.e. how the global error $E^n$ behaves as the time step
$k$ is decreased. It was seen that finite difference scheme usually
behave in such a way that $|E^n| = O(k^p)$, where $p \geq  1$ is the
order of accuracy. This is, of course, one way of looking at how well a
scheme performs. For the oscillator, it may be useful to measure the
degree of accuracy in the frequency domain. Solution
[\[eq:sol_z\_SHO\]](#eq:sol_z_SHO){reference-type="eqref"
reference="eq:sol_z_SHO"} suggests that the solutions are oscillating
with natural frequency $$\begin{equation\}
\label{eq:ErrFreqs}
    \omega = \frac{1}{k}\arctan{\frac{\omega_0^2k^2\sqrt{\frac{4}{\omega_0^2 k^2}-1}}{(2-\omega_0^2 k^2)}} = \omega_0 \left(1 + \frac{\omega_0^2 k^2}{24} +  \frac{3 \omega_0^4 k^4}{640} + O(\omega_0^6k^6) \right)\end{equation}
$$
and hence the natural frequency computed by scheme
[\[eq:Scheme1\]](#eq:Scheme1){reference-type="eqref"
reference="eq:Scheme1"} is second-order accurate compared to the natural
frequency of the continuous system. See also left panel of Fig.
[1.6](#eq:FreqWarping){reference-type="ref" reference="eq:FreqWarping"}.
The error in frequency can be quite audible, as $\omega_0$ approaches
the limit of stability $\omega_{max}=2/k$. In practice, the cent
deviation from [\[eq:ErrFreqs\]](#eq:ErrFreqs){reference-type="eqref"
reference="eq:ErrFreqs"} is given by

$$\begin{equation}
1200 \log_{2}\frac{\omega}{\omega_0} = 1200\log_2\left(1 + \frac{\omega_0^2 k^2}{24} +  \frac{3 \omega_0^4 k^4}{640} + O(\omega_0^6k^6) \right),
\end{equation}
$$
and this is shown in the right panel of Fig.
[1.6](#eq:FreqWarping){reference-type="ref" reference="eq:FreqWarping"}.

![Frequency warping (left) and cent deviation (right) of scheme
[\[eq:Scheme1\]](#eq:Scheme1){reference-type="eqref"
reference="eq:Scheme1"}. In the left panel, three natural frequencies
are shown: 100 rad/s (solid), 200 rad/s (dahsed), 300 rad/s
(dash-dotted). The right panel shows cent deviation up to
$O(\omega_0^6k^6)$ (solid), and exact
(dashed).](Figures/FreqWarp.pdf){#eq:FreqWarping width="\\linewidth"}

The frequency warping effects are quite evident. It may be preferable,
then, to construct schemes with a higher accuracy. This, of course,
cannot be done by merely using difference operators with a larger
stencil, such as the one given in
[\[eq:FourthOrderdtt\]](#eq:FourthOrderdtt){reference-type="eqref"
reference="eq:FourthOrderdtt"}: this would result in unstable behaviour!
A different approach, known as *modified equation method*, can be
constructed starting from Taylor-series arguments. One has
$$\begin{equation}\label{eq:dttExpasion}
    \delta_{tt}= \sum_{l=1}^{\infty}\frac{2k^{2(l-1)}}{(2l)!}\frac{d^{2l}}{dt^{2l}}= \frac{d^2}{dt^2} + \frac{k^2}{12}\frac{d^4}{dt^4} + \frac{k^4}{360}\frac{d^6}{dt^6} + O(k^6)\end{equation}$$
Considering for the moment the expansion up to the term $l=2$, it is
natural to add a term $-\frac{k^2}{12}\delta_{tt}\delta_{tt}x^n$ on the
left-hand side of [\[eq:Scheme1\]](#eq:Scheme1){reference-type="eqref"
reference="eq:Scheme1"}, in order to cancel the $O(k^2)$ error. Then,
one may use $\delta_{tt}= -\omega_0^2$ twice, to get

$$\begin{equation}
\delta_{tt}x^n = \frac{1}{k^2}\left(-\omega_0^2k^2 + \frac{\omega_0^4 k^4}{12} \right) x^n,
\end{equation}
$$
and of course the local truncation error $\varepsilon^n$ is now
$O(k^4)$. Considering now the next term in the series, proportional to
$k^4$, and using $\delta_{tt}= -\omega_0^2$ three times, one gets

$$\begin{equation}
\delta_{tt}x^n = \frac{1}{k^2}\left(-\omega_0^2k^2 + \frac{\omega_0^4 k^4}{12} - \frac{\omega_0^6k^6}{320}\right) x^n,
\end{equation}
$$
and this approximation is $O(k^6)$. Of course, one may go on and add
more and more terms to the series. The expansion can be rewritten
conveniently as

$$\begin{equation}
\frac{1}{k^2}\left(-\omega_0^2k^2 + \frac{\omega_0^4 k^4}{12} - \frac{\omega_0^6 k^6}{320} + ...\right) = \frac{2}{k^2}\left(-1 + \underbrace{1 - \frac{\omega_0^2k^2}{2} + \frac{\omega_0^4k^4}{4!} - \frac{\omega_0^6k^6}{6!} + ...}_{\cos(\omega_0k)} \right),
\end{equation}
$$
and thus, adding infinite terms to the series results in
$$\begin{equation}\label{eq:SHOexactNum}
    \delta_{tt}x^n =  \frac{2}{k^2}\left(-1 + \cos(\omega_0 k) \right)x^n\end{equation}$$
Since the series expansion converges to a known function in this case,
we can say that scheme
[\[eq:SHOexactNum\]](#eq:SHOexactNum){reference-type="eqref"
reference="eq:SHOexactNum"} solves
[\[eq:SHOreal\]](#eq:SHOreal){reference-type="eqref"
reference="eq:SHOreal"} *exactly*. Of course, this is somewhat too
strong a statement: following the discussion in Sec.
[1.1.1.4](#sec:Init){reference-type="ref" reference="sec:Init"}, we know
that the scheme is only going to be as accurate as its initial
conditions in this case. However, scheme
[\[eq:SHOexactNum\]](#eq:SHOexactNum){reference-type="eqref"
reference="eq:SHOexactNum"} does not introduce errors in the frequency
domain.

## Loss {#sec:LossSHO}

![Time evolution and phase portraits of lightly damped oscillator
((a),(b)), and overdamped oscillator ((c),(d)). The natural frequency is
$\omega_0=100$ rad/s, and the decay times are $\tau_{60}=5$ s (lightly
damped), 0.05 s (overdamped). Initial conditions are $x_0=-0.01$ m,
$v_0=0.04$ m/s.](Figures/PhaseSpLoss.png){#fig:dampedOsc
width="0.9\\linewidth"}

We are now going to study the behaviour of the oscillator in the
presence of dissipative forces. These are always present in some form in
real systems. Terms proportional to the velocity are usually a good
starting point to model viscous damping. For the oscillator, a
modification of [\[eq:SHOreal\]](#eq:SHOreal){reference-type="eqref"
reference="eq:SHOreal"} can then be given as $$\begin{equation}\label{eq:SHOLoss}
    \frac{d^2 x}{dt^2} = -\omega_0^2 x - 2c \frac{dx}{dt}.\end{equation}$$ 

Here,
$c \geq 0$ is a loss parameter (assumed constant, and measured in
s$^{-1}$). We remark that
[\[eq:SHOLoss\]](#eq:SHOLoss){reference-type="eqref"
reference="eq:SHOLoss"} is still a linear, time invariant system, and we
may infer its stability properties via *ansatz*
[\[eq:ansatz\]](#eq:ansatz){reference-type="eqref"
reference="eq:ansatz"}. Thus, $$\begin{equation}\label{eq:LossAnsatz}
    \hat x e^{st}\left(s^2 + 2c s + \omega_0^2 \right) = 0, \,\, \text{implying}\,\, s_{\pm} = -{c} \pm \sqrt{{c}^2 - \omega_0^2}\end{equation}$$
The qualitative behaviour of the oscillator will depend on the value of
$c$ compared to $\omega_0$, i.e. whether the square root in the
expression for $s_\pm$ is real or pure imaginary. Two cases of interest
may be extracted as follows:

1.  $c<\omega_0$. In this case, the oscillator is only lightly damped,
    and $s_\pm = -c \pm j \sqrt{\omega_0^2-c^2}$. Remembering the
    definition of $s$ in [\[eq:LapT\]](#eq:LapT){reference-type="eqref"
    reference="eq:LapT"}, one may extract $\sigma = - c$,
    $\omega = \sqrt{\omega_0^2-c^2}$. The solution to
    [\[eq:SHOLoss\]](#eq:SHOLoss){reference-type="eqref"
    reference="eq:SHOLoss"} may then be written as
    $x(t) = A_+ e^{s_+ t} + A_- e^{s_- t} = e^{-c t}\left(A_+ e^{j \omega t} + A_- e^{-j\omega t} \right)$.
    As per usual, the complex constants $A_+,A_-$ are determined from
    the intial conditions. The solution is in this case the product of
    an oscillating solution, times an exponentially damped envelope. The
    frequency of vibration is
    $\omega = \omega_0\sqrt{1 - c^2/\omega_0^2}\approx \omega_0\left(1 - \frac{c^2}{2\omega_0^2} \right)$
    and is thus lower than the natural frequency of the undamped
    oscillator. In this case, motion is not strictly periodic, since the
    mass is never really going to extend as far at each oscillation
    because of the energy given away to losses. However, it still makes
    sense to speak of period of vibration, as $\tau = 2\pi / \omega.$

2.  $c>\omega_0$. In this case, $s_\pm$ are both real, and negative,
    i.e. $s_\pm < 0$. Thus, the solutions are still "physical" in that
    they die out exponentially as time increases. However, the mass does
    not oscillate. This case is sometimes referred to as *overdamped
    oscillator.*

Fig. [1.7](#fig:dampedOsc){reference-type="ref"
reference="fig:dampedOsc"} shows illustrative examples of such cases.
Note that, in phase space, the orbits spiral toward the centre, and
motion is strictly speaking not periodic. A third case (*critically
damped*) is obtained whenever $c=\omega_0$. This case serves as a
mathematical "boundary" between the overdamped and lightly damped cases,
and is never really realised in practice. A useful quantity to quantify
loss is the *decay time* $\tau_{60}$. This is defined as the time taken
by the oscillator to reduce its amplitude of vibration by 60 dB. In the
lightly damped case, the amplitude envelope is simply $e^{-ct}$. Thus,
the implicit definition of $\tau_{60}$ is obtained as $$\begin{equation}\label{eq:tau60}
    -60 = 20 \log_{10}e^{-c\tau_{60}}, \,\, \text{and upon inversion: } \,\, \tau_{60} = \frac{3}{c} \ln(10).\end{equation}$$
This shows that the loss constant $c$ is most easily interpreted as a
function of the decay time.

### Energy analysis

Qualitative results on stability can be obtained via energy analysis.
Multiplying [\[eq:SHOLoss\]](#eq:SHOLoss){reference-type="eqref"
reference="eq:SHOLoss"} by $m \frac{dx}{dt}$ and using the same
indentities as for [\[eq:En1\]](#eq:En1){reference-type="eqref"
reference="eq:En1"}, we get $$\begin{equation}\label{eq:EnBalLoss}
    \frac{d}{dt}\left( \frac{m}{2} \left(\frac{dx}{dt}\right)^2 + \frac{K x^2}{2}   \right) = -Q(t) \triangleq - 2mc \left( \frac{dx}{dt} \right)^2 \leq 0,\end{equation}$$
implying that the total energy is not increasing. Here, $Q(t) \geq 0$ is
the power dissipated by the oscillator. Thus, bounds
[\[eq:bnds\]](#eq:bnds){reference-type="eqref" reference="eq:bnds"} hold
here as well. It is somewhat harder to draw more quantitative results
here, without knowledge on the form of $x(t).$ However, it is easy to
draw trajectories in phase space: instead of a closed loop, the mass now
spirals toward the centre as a result of losses.

### A finite difference scheme

As a discretisation for
[\[eq:SHOLoss\]](#eq:SHOLoss){reference-type="eqref"
reference="eq:SHOLoss"} is obtained as $$\begin{equation}\label{eq:FDSHOLoss}
    \delta_{tt}x^n = -\omega_0^2 x^n -2c \delta_{t\cdot}x^n,\end{equation}$$ 

yielding
an update equation

$$\begin{equation}
\left(ck+1\right) x^{n+1} = (2-\omega_0^2 k^2) x^n + (ck - 1)x^n.
\end{equation}
$$
Using the definition of local truncation error $\varepsilon^n$, as per
[\[eq:LTEdef\]](#eq:LTEdef){reference-type="eqref"
reference="eq:LTEdef"}, one gets $\varepsilon^n = O(k^2)$, $\forall n$,
showing that the LTE is second-order accurate. Stability may be inferred
using either frequency-domain analysis, or energy methods. Application
of the former via the ansatz $x^n = \hat x z^n$ gives
$$\begin{equation}\label{eq:SolFdLoss}
    \hat x z^n \left((1+ck)z - (2-\omega_0^2 k^2) -(ck-1)z^{-1}\right), \,\, z_\pm = \frac{2-\omega_0^2k^2 \pm \sqrt{(2-\omega_0^2k^2)^2 - 4(1-ck)(1+ck)}}{2(1+ck)},\end{equation}$$
and, using the Schur-Cohn condition
[\[eq:SchurCohnStab\]](#eq:SchurCohnStab){reference-type="eqref"
reference="eq:SchurCohnStab"}, $|z_\pm|<1$ if and only if

$$
\begin{equation}
k < \frac{2}{\omega_0}
\end{equation}
$$ 

that is the same as
[\[eq:StabCondSHO\]](#eq:StabCondSHO){reference-type="eqref"
reference="eq:StabCondSHO"}. The same condition may be arrived at via
energy analysis. To that end, multiply
[\[eq:FDSHOLoss\]](#eq:FDSHOLoss){reference-type="eqref"
reference="eq:FDSHOLoss"} by $\delta_{t\cdot}x^n$,

$$\begin{equation}
\delta_{t\cdot}x^n \, \delta_{tt}x^n = - \delta_{t\cdot}x^n \, \omega_0^2 x^n - 2c (\delta_{t\cdot}x^n)^2.
\end{equation}
$$
Using identities [\[eq:IdsFD\]](#eq:IdsFD){reference-type="eqref"
reference="eq:IdsFD"}, and multiplying by the mass $m$ to restore units
of energy, one gets

$$\begin{equation}
\delta_{t+}\left( \frac{m}{2}{(\delta_{t-}x^n)^2} + \frac{K}{2} {x^n e_{t-}x^n} \right) = - 2mc (\delta_{t\cdot}x^n)^2 \leq 0,
\end{equation}
$$
which is a discrete counterpart of
[\[eq:EnBalLoss\]](#eq:EnBalLoss){reference-type="eqref"
reference="eq:EnBalLoss"}. Thus, the discrete energy is non-increasing,
and when the total energy is itself non-negative, boundedness of the
solution results. Thus, the stability condition
[\[eq:StabCondSHO\]](#eq:StabCondSHO){reference-type="eqref"
reference="eq:StabCondSHO"} is recovered in this case as well. Of
course, [\[eq:SHObound\]](#eq:SHObound){reference-type="eqref"
reference="eq:SHObound"} holds in this case too. The numerical decay
time may be established via knowledge of the solutions $z_\pm$ in
[\[eq:SolFdLoss\]](#eq:SolFdLoss){reference-type="eqref"
reference="eq:SolFdLoss"}. Assuming oscillating behaviour, $z_\pm$ are
complex conjugates with absolute value

$$\begin{equation}
|z_\pm| = \sqrt{\frac{1-ck}{1+ck}}.
\end{equation}
$$ 

Thus, the numerical decay time
index $n_{60}$ is given by

$$\begin{equation}
-60 = 20\log_{10}\left({\frac{1-ck}{1+ck}}\right)^{n_{60}/2},\,\,\,\, n_{60}k = \frac{6k \ln(10)}{\ln \frac{1+ck}{1-ck}}\approx \tau_{60} - c k^2 \ln(10),
\end{equation}
$$
showing that the numerical decay time (in seconds) is $O(k^2)$ compared
to the exact decay time $\tau_{60}$ defined in
[\[eq:tau60\]](#eq:tau60){reference-type="eqref" reference="eq:tau60"}.

### Higher-order schemes

Higher-order accurate schemes may of course be obtained in this case as
well. To that end, consider the definition of the LTE, as

$$\begin{equation}
\left(\delta_{tt}+2c\delta_{t\cdot}\right) x(t_n) = - \omega_0^2 x(t_n) + \varepsilon^n,
\end{equation}
$$
where $x(t_n)$ is the true solution. Expanding in a Taylor series, one
has

$$\begin{equation}
\left(\frac{d^2}{dt^2} + 2c \frac{d}{dt} \right)\left(1 + \frac{k^2}{6}\frac{d^2}{dt^2}\right)x(t_n) - \frac{k^2}{12}\frac{d^4}{dt^4}x(t_n) + O(k^4) = - \omega_0^2 x(t_n) + \varepsilon^n.
\end{equation}
$$
This suggests the use the following modified scheme, in order to cancel
the terms proportional to $k^2$:

$$\begin{equation}
\left(\delta_{tt}+ 2c \delta_{t\cdot}\right)\left(1 - \frac{k^2}{6}\delta_{tt}\right)x^n + \frac{k^2}{12}\delta_{tt}\delta_{tt}x^n = -\omega_0^2 x^n.
\end{equation}
$$
Since $\delta_{tt}+ 2c \delta_{t\cdot}= -\omega_0^2 + O(k^2)$, the
scheme above can be written as (to the order $O(k^4)$):

$$\begin{equation}
\left(\delta_{tt}+ 2c \delta_{t\cdot}\right)x^n - \frac{k^2}{6}(-\omega_0^2) \delta_{tt}x^n + \frac{k^2}{12}\delta_{tt}\delta_{tt}x^n = -\omega_0^2 x^n.
\end{equation}
$$
The hard bit left is to find a suitable approximation to
$\delta_{tt}\delta_{tt}$, involving at most a stencil of width 2. This
can be accomplised in the following way $$\begin{equation}\label{eq:dttsq}
    \delta_{tt}\delta_{tt}\approx \left(-\omega_0^2 e_{t-}-2c \delta_{t-}\right)\left(-\omega_0^2 e_{t+}-2c \delta_{t+}\right) = \omega_0^4 + 4c \omega_0^2 \delta_{t\cdot}+ 4c^2 \delta_{tt}\end{equation}$$
Putting it all together, one obtains a fourth-order accurate
approximation to [\[eq:SHOLoss\]](#eq:SHOLoss){reference-type="eqref"
reference="eq:SHOLoss"} as

$$\begin{equation}
\left(1 + \frac{k^2}{6}(\omega_0^2 + 2c^2) \right)\delta_{tt}x^n = -\omega_0^2\left(1 + \frac{\omega_0^2 k^2}{12} \right) x^n - 2c \left(1 + \frac{\omega_0^2k^2}{6} \right)\delta_{t\cdot}x^n.
\end{equation}
$$
Higher-order accurate schemes can may be obtained this way, i.e. finding
approximations to $\delta_{tt}^p$, involving only operators of width 2.
A sketch of the idea is given briefly here. From
[\[eq:dttsq\]](#eq:dttsq){reference-type="eqref" reference="eq:dttsq"},
one may construct $\delta_{tt}^3$ in the following way:
$$\begin{aligned}
    \delta_{tt}\delta_{tt}\delta_{tt}\approx \left( \omega_0^4 + 4c \omega_0^2 \delta_{t\cdot}+ 4c^2 \delta_{tt}\right)\delta_{tt}\approx                                              
    \left( \omega_0^4e_{t-}+ 4c \omega_0^2 \delta_{t-}+ 4c^2 (-\omega_0^2e_{t-}-2c\delta_{t-})\right)\delta_{tt}\approx                                          \\\left( \omega_0^4e_{t-}+ 4c \omega_0^2 \delta_{t-}+ 4c^2 (-\omega_0^2e_{t-}-2c\delta_{t-})\right)(-\omega_0^2e_{t+}-2c\delta_{t+})= \\
    (-\omega_0^6 + 4c^2 \omega_0^4) + \left(-6 c \omega_0^4 + 16 c^3 \omega_0^2 \right)\delta_{t\cdot}+ \left( -8 c^2 \omega_0^2 + 16 c^4\right)\delta_{tt}, \end{aligned}$$
showing that $\delta_{tt}^3$ can be approximated using a stencil of
width 2. One may of course use the modified equation technique described
above to any desired order. Luckily, the oscillator with loss also
posseses an exact solution, where "exact" is intended in the same way as
for [\[eq:SHOexact\]](#eq:SHOexact){reference-type="eqref"
reference="eq:SHOexact"} (i.e. exact up to the accuracy order of the
intial conditions). Considering again the continuous equation with loss,
[\[eq:SHOLoss\]](#eq:SHOLoss){reference-type="eqref"
reference="eq:SHOLoss"}, under the following transformation

$$\begin{equation}
X(t) = e^{ct} x(t),
\end{equation}
$$ 

one gets

$$\begin{equation}
\frac{d^2X}{dt^2} + \left(\omega_0^2-c^2 \right) X = 0.
\end{equation}
$$ 

Thus, the
exact scheme
[\[eq:SHOexactNum\]](#eq:SHOexactNum){reference-type="eqref"
reference="eq:SHOexactNum"} for the undamped oscillator can be applied
to the transformed variable $X$. When transformed back to $x$, this
gives $$\begin{equation}\label{eq:SHO2}
    % \left(\dtt + \frac{2e^{-c k}}{k^2} \left(e^{c k}-\cos\left(\sqrt{\omega_0^2-c^2} k\right)\right) + \frac{e_{t-}}{k^2} (e^{-2c k}-1)\right)x^n = 0.
    \delta_{tt}x^n = \left(-\frac{2}{k^2}\left(1 - \cos\left((\sqrt{\omega_0^2-c^2} \, k\right) \right)- \frac{e_{t+}(e^{ck}-1) + e_{t-}(e^{-ck}-1)}{k^2}\right)x^n\end{equation}$$
This scheme solves [\[eq:SHOLoss\]](#eq:SHOLoss){reference-type="eqref"
reference="eq:SHOLoss"} exactly, in particular, the frequency of
oscillation and the numerical decay time are exact.

## Forced Oscillations

The equation of the oscillator including loss and source terms is given
as $$\begin{equation}\label{eq:SHOForced}
    \frac{d^2 x}{dt^2} = -\omega_0^2 x - 2c \frac{dx}{dt} + f(t),\end{equation}$$
where $f(t)$ is a time-dependent force per unit mass. Energy analysis
leads here to the following energy balance

$$\begin{equation}
\frac{d}{dt}\left( \frac{m}{2} \left(\frac{dx}{dt}\right)^2 + \frac{K x^2}{2}   \right) = - Q(t) + P(t),
\end{equation}
$$
where $Q(t) = 2mc \left(\frac{dx}{dt}\right)^2$ is the dissipated power,
and where $P(t) = m\frac{dx}{dt}f(t)$ is the injected power.

In this case, the system is still linear, but it is not time invariant.
Solutions via transform techniques can be obtained, involving the use of
the one-sided Laplace transform
[\[eq:LapT\]](#eq:LapT){reference-type="eqref" reference="eq:LapT"} so
to incorporate the effects of the initial conditions and of the external
forcing. In this case, substitution of the simpler *ansatz*
[\[eq:ansatz\]](#eq:ansatz){reference-type="eqref"
reference="eq:ansatz"} is not possible, because of the presence of the
forcing term. Considering $t\geq 0$, the application of the one-sided
Laplace transform in
[\[eq:SHOForced\]](#eq:SHOForced){reference-type="eqref"
reference="eq:SHOForced"} gives

$$\begin{equation}
\int_0^{\infty} \left( \frac{d^2 x}{dt^2} + \omega_0^2 x + 2c \frac{dx}{dt} - f(t)\right)e^{-st} \, \text{d}t = 0
\end{equation}$$

Using integration by parts, one gets (remember that we defined
$x_0 = x(0), v_0 = dx(0)/dt$):

$$
\begin{equation}
\hat x(s)\left( s^2 + 2cs + \omega_0^2\right) = (s+2c)x_0 + v_0 + \hat{f}(s),
\end{equation}
$$
and thus

$$
\begin{equation}
\hat x(s) = \frac{(s+2c)x_0 + v_0 + \hat{f}(s)}{s^2 + 2cs + \omega_0^2}.
\end{equation}
$$
The expression above may be decomposed into the *transient* and *forced
response*. The transient is given by the contribution of the initial
conditions only, without external forcing, so:

$$\begin{equation}
\hat x(s) = \hat x_{tr}(s) + \hat x_{fr}(s) = \frac{(s+c)x_0 + (v_0 + cx_0)}{s^2 + 2cs + \omega_0^2} + \frac{\hat{f}(s)}{s^2 + 2cs + \omega_0^2}.
\end{equation}
$$
This shows that the contributions of the transient and of the forced
response are independent of each other: they add up in the final
response. The solution in the time domain is obtained upon inversion of
$\hat x(s).$ It is best to write the denominator of $\hat x(s)$ as
$(s+c)^2 + \left(\sqrt{\omega_0^2-c^2}\right)^2$, since this is the form
reported in [\[eq:LaplTtable\]](#eq:LaplTtable){reference-type="eqref"
reference="eq:LaplTtable"}. Using these (with
$a = \sqrt{\omega_0^2-c^2}$), the solution to the transient is

$$\begin{equation}
x_{tr}(t) = e^{-ct}\left(x_0 \cos\left(\sqrt{\omega_0^2-c^2}\,\,t\right) + \frac{v_0 + c x_0}{\sqrt{\omega_0^2-c^2}}\sin\left(\sqrt{\omega_0^2-c^2}\,\,t\right)\right),
\end{equation}
$$
which is of course the same as
[\[eq:LossAnsatz\]](#eq:LossAnsatz){reference-type="eqref"
reference="eq:LossAnsatz"} (we did not go through the substitution of
the initial conditions there, but the result is the same). For the
forced response, we may employ the *convolution* property of the Laplace
transform, which states that

$$\begin{equation}
\mathcal{L}^{-1}(\hat a(s)\hat b(s)) = \int_{0}^t a(t-u)b(u) \, \mathop{}\!\mathrm{d}u,
\end{equation}
$$
and thus, using $\hat a = 1/((s+c)^2 + (\sqrt{\omega_0^2-c^2})^2)$,
$\hat b = \hat f$, one gets $$\begin{equation}\label{eq:steadyStSHO}
    x_{fr}(t) = \int_0^t \frac{e^{-c(t-u)}}{\sqrt{\omega_0^2-c^2}}\sin\left(\sqrt{\omega_0^2-c^2}\,\,(t-u)\right) f(u) \, \mathop{}\!\mathrm{d}u.\end{equation}$$
This integral is not generally computable analytically. However, it also
encapsulates the idea that the forced response to *any* forcing can be
obtained as the convolution with the *impulse response*. To that end,
considering $f(t)=\delta(t-t_0)$, with delta being Dirac's delta here,
and denoting as per usual the initial time by $t_0$, one gets from
[\[eq:steadyStSHO\]](#eq:steadyStSHO){reference-type="eqref"
reference="eq:steadyStSHO"}: $$\begin{equation\}
\label{eq:GreenSHO}
    G(t|t_0) = \frac{e^{-c(t-t_0)}}{\sqrt{\omega_0^2-c^2}}\sin\left(\sqrt{\omega_0^2-c^2}\,\,(t-t_0)\right), \,\,\,\, t\geq t_0,$$\end{equation}

and $G$ is zero for $t<t_0$. The symbol $G(t|t_0)$ denotes the *Green's
function* of the harmonic oscillator, i.e. the response to a Dirac
impulse at $t=t_0$. In particular, for zero initial conditions, the
forced response is also the total response of the system.

![Time evolution of forced oscillator. Dashed-dotted line is response to
initial conditions, dahsed line is Green's function
[\[eq:GreenSHO\]](#eq:GreenSHO){reference-type="eqref"
reference="eq:GreenSHO"}, solid line is total response. The oscillator
is activated with a Dirac impulse at $t_0=0$. Initial conditions are
$x_0=-0.01$ m, $v_0=0.04$ m/s. The natural frequency of the oscillator
is $\omega_0 = 100$ rad/s, and the decay time is $\tau_{60} = 5$
s.](Figures/ForcedOsc.png){width="0.9\\linewidth"}

The stability of system
[\[eq:SHOForced\]](#eq:SHOForced){reference-type="eqref"
reference="eq:SHOForced"} may be adapted to include the effects of
external forcing. We are not going to bother as much here, and we will
assume that the solution will not "blow up" in a finite time if the
source remains bounded. In particular, if the source is itself bounded,
we remark that integral
[\[eq:steadyStSHO\]](#eq:steadyStSHO){reference-type="eqref"
reference="eq:steadyStSHO"} remains bounded.

### Response to harmonic input forcing

Consider now a harmonic input of the form $$\begin{equation}\label{eq:harmoForceLinear}
    f(t) = F e^{j\omega t},\end{equation}$$ 

with $F \in \mathbb{C}$ being a complex
forcing amplitude, and $\omega \in \mathbb{R}^+_0$ being the input
forcing radian frequency (not to be confused with the natural frequency
of the oscillator). We are only going to assume that the oscillator will
eventually fall into a steady-state here, since from the previous
discussion we know that the transient response will die out after
sufficient time has elapsed, and since the forcing is of harmonic type.
According to
[\[eq:steadyStSHO\]](#eq:steadyStSHO){reference-type="eqref"
reference="eq:steadyStSHO"}, one may compute the steady-state via the
convolution integral. However, in this case one may equivalently assume
that the steady state vibrates at the frequency of the input, thus

$$\begin{equation}
x(t) = X e^{j\omega t},
\end{equation}
$$ 

where we removed the index $st$ since we are
now assuming that the transient has completely died out, and thus we can
identify the whole solution $x(t)$ of
[\[eq:SHOForced\]](#eq:SHOForced){reference-type="eqref"
reference="eq:SHOForced"} with the steady-state. $X$ is here a complex
amplitude. Substituting into the equation of motion results in
$$\begin{equation}\label{eq:TransXF}
    \frac{X}{F} = \frac{1}{\omega_0^2-\omega^2+2jc\omega} =  \frac{\omega_0^2-\omega^2-2jc\omega}{(\omega_0^2-\omega^2)^2+4c^2\omega^2}.\end{equation}$$
Since $X,F$ are complex constants, the tangent of the phase angle is
obtained as the ratio between the imaginary and real parts, i.e.

$$\begin{equation}
\tan \left(\angle \frac{X}{F}\right) = -\frac{2c\omega}{\omega_0^2 - \omega^2}
\end{equation}
$$
One should pay attention to the sign of the denominator when inverting
the tangent function. Hence, the phase starts out at $0$ when
$\omega \approx 0$; then it reaches $-\pi/2$ when
$\omega \approx \omega_0$, and then approaches $-\pi$ as
$\omega \rightarrow \infty$. The absolute value of the transfer function
is obtained as $$\begin{equation\}
\label{eq:XFlinearOsc}
    \left|\frac{X}{F}\right|  = \left((\omega_0^2-\omega^2)^2+4c^2\omega^2\right)^{-1/2}\end{equation}
$$
which has a maximum at
$\omega = \omega_0 \left(1-2(c/\omega_0)^2 \right)^{1/2}$. The maximum
is

$$\begin{equation}
\left|\frac{X}{F}\right|(\omega = \omega_0 \left(1-2(c/\omega_0)^2 \right)^{1/2}) \approx \frac{1}{2c\omega_0},
\end{equation}
$$
where factors of the order $O(c^4)$ were disregarded. The *bandwidth* of
the transfer function is defined as the interval in frequency occurring
between the frequencies $\omega_+,\omega_-$ that are found at
$\frac{1}{\sqrt{2}}$ the maximum (these are the *half-power points*,
since power is the square of the absolute value). To obtain
$\omega_\pm$, one simply uses this defintion, hence:

$$\begin{equation}
\left((\omega_0^2-\omega^2)^2+4c^2\omega^2\right)^{-1/2} = \left(2\sqrt{2} c\omega_0\right)^{-1}.
\end{equation}
$$
Solving for $\omega$, and disregarding small terms, one gets

$$\begin{equation}
\omega_\pm \approx \omega_0\left(1 \pm {c} \right), \quad \rightarrow (\omega_+ - \omega_-) \approx 2c.
\end{equation}
$$
This shows that, for small damping values, one may recover the value of
the decay time from the frequency response, after measuring the bandwith
of the peak in the transfer function. As a cautionary note, the transfer
function $X/F$ was obtained here in the case of sinusoidal input
forcing, by "scanning" the frequency axis. The same transfer function
may however be obtained via the impulse response
[\[eq:GreenSHO\]](#eq:GreenSHO){reference-type="eqref"
reference="eq:GreenSHO"}. To show this, it is sufficient to compute a
Fourier transform and, again, we may resort to tables for this.
Considering the transforms
[\[eq:FouTtable\]](#eq:FouTtable){reference-type="eqref"
reference="eq:FouTtable"} (with $a = \sqrt{\omega_0^2 - c^2}$,
$x(t) = e^{-ct}$), one gets for the impulse response
[\[eq:GreenSHO\]](#eq:GreenSHO){reference-type="eqref"
reference="eq:GreenSHO"}

$$\begin{equation}
\mathcal{F}\left\{G(t|t_0=0)\right\}(\omega) = \frac{1}{\sqrt{2\pi}}\frac{1}{(\omega_0^2-\omega^2)+2jc\omega} = \frac{1}{\sqrt{2\pi}}\frac{\omega_0^2-\omega^2-2jc\omega}{(\omega_0^2-\omega^2)^2+4c^2\omega^2},
\end{equation}
$$
that is the same as [\[eq:TransXF\]](#eq:TransXF){reference-type="eqref"
reference="eq:TransXF"} (up to a constant factor). Thus, knowledge of
the impulse response is equivalent to "scanning" the frequency axis one
frequency at a time. This equivalence is often employed experimentally,
where the impulse response is obtained via deconvolution of appropriate
sine sweeps.

#### Mechanical impedance

Though [\[eq:TransXF\]](#eq:TransXF){reference-type="eqref"
reference="eq:TransXF"} expresses the general relationship between
output displacement and input forcing, it may be preferable to obtain
the transfer function between output velocity and input forcing. Thus,

$$\begin{equation}
\frac{dx(t)}{dt} = j\omega e^{j\omega t} \triangleq V e^{j\omega t}.
\end{equation}
$$
Using this in [\[eq:TransXF\]](#eq:TransXF){reference-type="eqref"
reference="eq:TransXF"} results in

$$\begin{equation}
\frac{V}{F} = \frac{j\omega(\omega_0^2 - \omega^2)-2c\omega^2}{(\omega_0^2-\omega^2)^2+4c^2\omega^2}.
\end{equation}
$$
The ratio $V/F$ (often denoted $Y$) is the *mechanical admittance.* The
inverse $F/V$ (often denoted $Z$) is the *mechanical impedance.* The
reason why it may be preferable to work with the impedance (or the
admittance), rather than $X/F$, is that velocity and force and
power-conjugated quantities: in an energy framework, knowledge of the
impedance/admittance allows to describe mechanical systems in an
energy-consistent manner, as we will see in due course. For the phase
angle, one has

$$\begin{equation}
\tan \left(\angle \frac{V}{F}\right) = \frac{\omega_0^2 - \omega^2}{(-2c\omega)}.
\end{equation}
$$
Thus, in this case the phase angle starts out at $\pi/2$ when
$\omega \approx 0$, then it goes to zero when $\omega \approx \omega_0$,
and then it approaches $-\pi/2$ as $\omega \rightarrow \infty$. The
absolute value is given by

$$\begin{equation}
\left| \frac{V}{F} \right| = \omega \left((\omega_0^2-\omega^2)^2+4c^2\omega^2\right)^{-1/2}.
\end{equation}
$$
Differentiating with respect to $\omega$, one gets a maximum at
$\omega = \omega_0$. The maximum is

$$\begin{equation}
\left|\frac{V}{F}\right|(\omega = \omega_0) = (2c)^{-1}.
\end{equation}
$$

![Absolute values and phase angles of the transfer functions for the
linear oscillator. The natural frequency of the oscillator is
$\omega_0 = 100$ rad/s, the decay times are set as
$\tau_{60}=\left\{1.0,1.2,1.4,1.6,1.8,2.0\right\}$
s](Figures/ImpedanceSHO.png){#fig:LinearTransFunctPlots
width="0.85\\linewidth"}

#### Finite difference schemes

![Error slopes of scheme
[\[eq:shoFDForced\]](#eq:shoFDForced){reference-type="eqref"
reference="eq:shoFDForced"}. The oscillator is activated with a Dirac
delta at $t_0=0$. Initial conditions are $x_0=-0.01$ m, $v_0=0.04$ m/s.
The natural frequencies of the oscillator are selected as
$\omega_0 = 100$ rad/s (solid), $\omega_0 = 200$ rad/s (dashed),
$\omega_0 = 300$ rad/s (dash-dotted), and the decay time is
$\tau_{60} = 5$ s.](Figures/SlopesForced.pdf){#fig:ErrSlopesForced
width="0.85\\linewidth"}

Finite difference schemes may be obtained in this case by any
discretisation of $f(t)$, of the desired accuracy. For second-order
accuracy,
$f^n = \left\{f(t_n),\mu_{t\cdot}f(t_n), \mu_{tt}f(t_n) \right\}$, are
all valid discretisations. Hence, a suitable second-order accurate
scheme is obtained as 

$$\begin{equation}\label{eq:shoFDForced}
    \delta_{tt}x^n = -\omega_0^2 x^n -2c\delta_{t\cdot}x^n + f^n
\end{equation}$$

Initialisation should be performed in such a way that second-order
accuracy is preserved. In this respect, one may set $x^0 = x_0$. Then,
[\[eq:dtdInit\]](#eq:dtdInit){reference-type="eqref"
reference="eq:dtdInit"} is used again, to get

$$\begin{equation}
(\delta_{t+}- \frac{k}{2}\delta_{tt})x^0 = v_0.
\end{equation}
$$ 

One has

$$\begin{equation}
\delta_{tt}x^0 \approx -\omega_0^2 x_0 - 2c \delta_{t+}x^0 + f^0.
\end{equation}
$$
Using these, one can extract $x^1$ as

$$\begin{equation}
x^1 = x^0 + \frac{k v_0 + \frac{k^2}{2}\left(-\omega_0^2 x^0 + f^0 \right)}{1+kc}
\end{equation}
$$
When approximating a Dirac delta at $t_0=0$, one should use $f^0 = 2/k$.
Error slopes for scheme
[\[eq:shoFDForced\]](#eq:shoFDForced){reference-type="eqref"
reference="eq:shoFDForced"} are shown in Fig.
[1.9](#fig:ErrSlopesForced){reference-type="ref"
reference="fig:ErrSlopesForced"}.
