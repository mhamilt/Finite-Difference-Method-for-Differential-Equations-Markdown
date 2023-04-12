---
layout: post
title: Lecture 5
---

-   [Numerical methods for the wave
    equation](#numerical-methods-for-the-wave-equation)
    -   [Spatial difference operators](#spatial-difference-operators)
        -   [Frequency-domain intepretation of spatial difference
            operators](#sec:FDtransformationsSpatial)
        -   [Matrix form of difference
            operators](#matrix-form-of-difference-operators)
        -   [Numerical Boundary Conditions](#sec:WEQNmatrices)
        -   [Eigenvalues and eigenvectors of second difference
            matrices](#eigenvalues-and-eigenvectors-of-second-difference-matrices)
    -   [An explicit finite difference
        scheme](#an-explicit-finite-difference-scheme)
        -   [Stability via Energy
            Analysis](#stability-via-energy-analysis)
        -   [Accuracy](#accuracy)
        -   [Initialisation. Global error.
            ](#initialisation.-global-error.)
    -   [Loss and forcing](#loss-and-forcing)
        -   [Finite difference schemes](#finite-difference-schemes)
        -   [Spreading and interpolation](#spreading-and-interpolation)

# Numerical methods for the wave equation

The previous chapter has introduced the formalism of the one-dimensional
wave equation. In this chapter we explore the numerical methods that can
be employed to solve it. While exact solutions are available in this
case, these do not generalise other systems. In fact, the exact
solutions serve here as a useful benchmark against which to test the
numerical schemes presented. We begin the discussion of the numerical
mehtods illustrating the use of finite differences.

## Spatial difference operators

Just like time can be discretised by means of a sample rate (or
equivalently, a time step), the spatial domain may also be discretised
using an appropriate *grid spacing*. Let such spacing be uniform along
the domain of interest, so that the continuous domain
$\mathcal{D} = \{ x : x \in [0,L]\}$ is mapped onto the discrete domain
$\mathfrak{d} = \{ mh : 0 \leq m \leq M = L/h, m \in {\mathbb N}\}$. In
the same way as, in Chapter
[\[chap:Intro\]](#chap:Intro){reference-type="ref"
reference="chap:Intro"}, we approximated the continuous time function
$x(t)$ via the discrete time series $x^n$, here we are going to
approximate the continuous function $y(x)$ via the *grid function*
$y_m$. The identity, forward and backward shift operators are then given
as
$${1}{y}_m = {y}_m, \quad e_{x+}{y}_m = {y}_{m+1}, \quad e_{x-}{y}_m = {y}_{m-1}.$$
From these, one may define the difference operators, all approximating
the first spatial derivative, as

::: subequations
$$\begin{aligned}
        \delta_{x+}& = \frac{(e_{x+}- 1)}{h} \approx \frac{d}{dx},     \\ 
        \delta_{x-}& = \frac{(1 - e_{x-})}{h}\approx \frac{d}{dx},      \\ 
        \delta_{x\cdot}& = \frac{(e_{x+}- e_{x-})}{2h}\approx \frac{d}{dx} . 
    \end{aligned}$$
:::

An approximation to the second time derivative is constructed from the
above as
$$\delta_{xx}= \delta_{x+}\delta_{x-}\approx \frac{d^2}{dx^2}.$$
Averaging operators (all approximating the identity) are

::: subequations
$$\begin{aligned}
        \mu_{x+} & = \frac{(e_{x+}+ 1)}{2} \approx 1,    \\ 
        \mu_{x-} & = \frac{(1 + e_{x-})}{2} \approx 1,    \\ 
        \mu_{x\cdot}& = \frac{(e_{x+}+ e_{x-})}{2} \approx 1. 
    \end{aligned}$$
:::

Taylor series arguments can be used to infer the order of the
approximation. The truncation errors of the difference operators may
again be found by applying them to the continuous function $y(x)$, and
expanding in a Taylor series. Hence

::: subequations
[\[eq:ErrsSpatial\]]{#eq:ErrsSpatial label="eq:ErrsSpatial"}
$$\begin{aligned}
\
     \frac{dy(x_m)}{dx} - \delta_{t+}y(x_m) &= O(h), \\
     \frac{dy(x_m)}{dx} - \delta_{t-}y(x_m) &= O(h), \\   
     \frac{dy(x_m)}{dx} - \delta_{t\cdot}y(x_m) &= O(h^2),\\   
     \frac{d^2y(x_m)}{dx^2} - \delta_{tt}y(x_m) &= O(h^2),\end{aligned}$$
:::

and so on.

### Frequency-domain intepretation of spatial difference operators {#sec:FDtransformationsSpatial}

In Section
[\[sec:FDtransformations\]](#sec:FDtransformations){reference-type="ref"
reference="sec:FDtransformations"}, the temporal difference operators
were applied to the kernel of the Fourier transform. A similar idea can
be employed for the spatial difference operators, where of course the
Fourier transform is now
$$\mathcal{X}\left\{ y_m \right\} = \sum_{m=-\infty}^\infty y_m e^{-jmh\gamma},$$
where $\gamma$ is the wavenumber. The spatial difference operators
transform in a manner analougous to the temporal ones. Take the second
difference: $$\label{eq:DTFTdss}
\mathcal{X}\left\{\delta_{xx}y_m\right\} = \frac{e^{j\gamma h}-2+e^{-j\gamma h}}{h^2}\mathcal{X}\left\{y_m\right\}=\frac{2}{h^2}\left(\cos(\gamma k)-1 \right)\mathcal{X}\left\{y_m\right\} = -\frac{4}{h^2} \sin^2\left(\frac{\gamma h}{2}\right)\mathcal{X}\left\{y_m\right\}$$
Analogously, one has

::: subequations
$$\begin{aligned}
\mathcal{X}\left\{\mu_{x\cdot} y_m\right\} &= \cos(\omega k) \mathcal{X}\left\{y_m \right\}, \\
\mathcal{X}\left\{\mu_{xx} y_m\right\} &= \frac{1}{2}\left(\cos(\omega k) + 1\right)\mathcal{X}\left\{y_m \right\} \\
\mathcal{X}\left\{\delta_{x\cdot}y_m\right\} &= \frac{j}{k}\sin(\omega k) \mathcal{X}\left\{y_m \right\}\end{aligned}$$
:::

These identities will prove useful in the study of the stability of the
wave equation, as will be seen shortly.

### Matrix form of difference operators

Since $y_m$ is interpreted as a grid function, it is natural to think of
the spatial difference operators as matrices acting on the vector
${\bf y} \in \mathbb{R}^{M+1}$. Consider the operator $\delta_{x+}$. For
interior points, one has
$\delta_{x+}y_m = \left({\bf D}^+ {\bf y}\right)_m = \frac{y_{m+1}-y_m}{h}$.
This is $${\bf D}^+ {\bf y} = \frac{1}{h}\begin{bmatrix}
& & {\bf \ddots}  & \ddots & & & \\ 
& & & {\bf -1} & 1 & & &\\
& & & & {\bf -1} & 1 & & \\
& & & & & \ddots & \ddots & \\
\end{bmatrix}
\begin{bmatrix}
\vdots \\
y_m \\
y_{m+1}\\
\vdots
\end{bmatrix}.$$ Here, the boldface indicates the elements on the main
diagonal. Analogously, one can define matrices for the operator
$\delta_{x-}y_m = \left({\bf D}^- {\bf y}\right)_m = \frac{y_{m}-y_{m-1}}{h}$
$${\bf D}^- {\bf y} = \frac{1}{h}\begin{bmatrix}
& & \ddots  & \ddots & & & \\ 
& & & -1 & {\bf 1} & & &\\
& & & & -1 & {\bf 1} & & \\
& & & & & \ddots & \ddots & \\
\end{bmatrix}
\begin{bmatrix}
\vdots \\
y_{m-1}\\
y_{m}\\
\vdots
\end{bmatrix}.$$ The centered derivative is given by
$\delta_{x\cdot}y_m = \left({\bf D} {\bf y}\right)_m = \frac{y_{m+1} - y_{m-1}}{2h}$,
or $${\bf D} {\bf y} = \frac{1}{2 h}\begin{bmatrix}
& & \ddots  & \ddots & \ddots & & & &\\ 
& & & -1 & {\bf 0} & 1 & & & &\\
& & & & -1 & {\bf 0 } & 1 &  & &\\
& & & & & -1 & {\bf 0 } & 1 &  & \\
& & & & & & \ddots & \ddots &  \ddots & \\
\end{bmatrix}
\begin{bmatrix}
\vdots \\
y_{m-1}\\
y_{m}\\
y_{m+1} \\
\vdots
\end{bmatrix}.$$ The second derivative is given by
$\delta_{xx}y_m = \left({\bf D}^2 {\bf y}\right)_m = \frac{y_{m+1} -2y_{m} + y_{m-1}}{h^2}$,
so that $${\bf D}^2 {\bf y} = \frac{1}{h^2}\begin{bmatrix}
& & \ddots  & \ddots & \ddots & & & &\\ 
& & & 1 & {\bf -2} & 1 & & & &\\
& & & & 1 & {\bf -2 } & 1 &  & &\\
& & & & & 1 & {\bf -2 } & 1 &  & \\
& & & & & & \ddots & \ddots &  \ddots & \\
\end{bmatrix}
\begin{bmatrix}
\vdots \\
y_{m-1}\\
y_{m}\\
y_{m+1} \\
\vdots
\end{bmatrix}.$$ The averaging operators can be defined analogously.
Denoting them by the letter ${\bf W}$, one has

::: subequations
$$\begin{aligned}
{\bf W}^+  &= \frac{1}{2}\begin{bmatrix}
& & {\bf \ddots}  & \ddots & & & \\ 
& & & {\bf 1} & 1 & & &\\
& & & & {\bf 1} & 1 & & \\
& & & & & \ddots & \ddots & \\
\end{bmatrix}, \\
{\bf W}^-  &= \frac{1}{2}\begin{bmatrix}
& & \ddots  & \ddots & & & \\ 
& & & 1 & {\bf 1} & & &\\
& & & & 1 & {\bf 1} & & \\
& & & & & \ddots & \ddots & \\
\end{bmatrix}, \\
{\bf W} &= \frac{1}{2 }\begin{bmatrix}
& & \ddots  & \ddots & \ddots & & & &\\ 
& & & 1 & {\bf 0} & 1 & & & &\\
& & & & 1 & {\bf 0 } & 1 &  & &\\
& & & & & 1 & {\bf 0 } & 1 &  & \\
& & & & & & \ddots & \ddots &  \ddots & \\
\end{bmatrix}.\end{aligned}$$
:::

### Numerical Boundary Conditions {#sec:WEQNmatrices}

While the forms of the matrices above hold for central points, the
points near the boundary deserve special treatment. One needs of course
to discretise the prescribed boundary conditions. We are considering now
the two types of boundary conditions encountered in Section
[\[sec:WaveEquationsBCs\]](#sec:WaveEquationsBCs){reference-type="ref"
reference="sec:WaveEquationsBCs"}: Dirichlet ($y=0$) and Neumann
$\left(\frac{dy}{dx}=0\right)$.

#### Fixed boundary conditions

For conditions of Dirichlet type, one can apply $y_0=y_M=0$. Clearly one
need not store the end points, since these are identically zero. Then,
one may define $${\bf D}^-_d  = \frac{1}{h}\begin{bmatrix}
{\bf 1} &  &   &   \\ 
-1& {\bf 1} & & &    \\
&-1 & {\bf 1} & &    \\
& & \ddots & \ddots &    \\
& &  & -1&   
\end{bmatrix},$$ and note that this matrix is *rectangular*, with
dimensions $M \times M-1$. The subscript $d$ is for *Dirichlet*. The
elements on the main diagonal were again written in boldface. The matrix
${\bf D}^{+}$ can be defined simply as
$${\bf D}^+_{d} = -({\bf D}^-_d)^\intercal$$ or
$${\bf D}^+_d  = \frac{1}{h}\begin{bmatrix}
{\bf -1} &  1 &   &   &\\ 
& {\bf -1} & 1 & &  &  \\
& & {\bf -1} & 1&   & \\
& & & \ddots &  \ddots & \\
& & &  & {\bf -1}&   1
\end{bmatrix},$$ and is a rectangular matrix of dimensions
$M-1 \times M$. The second difference operator is then given by
composing the first difference matrices, as $$\label{eq:SecondDiffMat}
{\bf D}^2_d = {\bf D}^+_d \,\, {\bf D}^-_d$$ Note that this composition
is obvisously *not* symmetric! That is, the product
${\bf D}^-_d \,\, {\bf D}^+_d$ is not the same as
[\[eq:SecondDiffMat\]](#eq:SecondDiffMat){reference-type="eqref"
reference="eq:SecondDiffMat"}, since the resulting matrix has different
dimensions. The correct composition is given by
[\[eq:SecondDiffMat\]](#eq:SecondDiffMat){reference-type="eqref"
reference="eq:SecondDiffMat"}, since the resulting matrix is
$M-1 \times M-1$. Explicitly $$\label{eq:D2Diri}
{\bf D}^2_d = \frac{1}{h^2}\begin{bmatrix}
{\bf -2}& 1&   &  &  &  \\ 
-1& {\bf -2} & 1 &  &  &\\
& \ddots & \ddots & \ddots &  &      \\
& &  1 & {\bf -2}& 1   \\
& &   & 1 & {\bf -2}   \\
\end{bmatrix}.$$

#### Free boundary conditions (Non-centred)

Consider now an approximation to the Neumann condition as
$\delta_{x+}y_0 = \delta_{x-}y_M=0$. This gives $y_0=y_1$ and
$y_M=y_{M-1}$. Again, one need not store the end points $y_0$, $y_M$,
since their value is the same as the inner points $y_1$, $y_{M-1}$.
Then, one can again define rectangular matrices as the forward and
backward difference operators. For the backward difference, take
$${\bf D}^-_{n1}  = \frac{1}{h}\begin{bmatrix}
{\bf 0} &  &   &   \\ 
-1& {\bf 1} & & &    \\
&-1 & {\bf 1} & &    \\
& & \ddots & \ddots &    \\
& &  & 0&   
\end{bmatrix},$$ This is a $M \times M-1$ matrix, whose first and last
rows are filled with zeros. Again, boldaface symbols are on the main
diagonal. The subscript ${n1}$ indicates *first-order Neumann*
conditions. Just like previously, the forward difference operator is
obtained as $${\bf D}^+_{n1} = -({\bf D}^-_{n1})^\intercal,$$ or
$${\bf D}^+_{n1} = \frac{1}{h}\begin{bmatrix}
{\bf 0} &  1 &   &   &\\ 
& {\bf -1} & 1 & &  &  \\
& & {\bf -1} & 1&   & \\
& & & \ddots &  \ddots & \\
& & &  & {\bf -1}&   0
\end{bmatrix},$$ and is a $M-1\times M$ matrix. The second difference
operator is given by $$\label{eq:SecondDiffMatNeumann}
{\bf D}^2_{n1} = {\bf D}^+_{n1} \,\, {\bf D}^-_{n1},$$ or
$$\label{eq:D2Neum1}
{\bf D}^2_{n1} = \frac{1}{h^2}\begin{bmatrix}
{\bf -1}& 1&   &  &  &  \\ 
-1& {\bf -2} & 1 &  &  &\\
& \ddots & \ddots & \ddots &  &      \\
& &  1 & {\bf -2}& 1   \\
& &   & 1 & {\bf -1}   \\
\end{bmatrix},$$ and is $M-1 \times M-1$.

#### Free boundary conditions (Centred)

A second-order accurate implementation of the Neumann condition results
from the application of $\delta_{x\cdot}y_0=\delta_{x\cdot}y_M=0$. In
this case, it results that $y_{-1}=y_1$, $y_{M-1}=y_{M+1}$, and hence,
application of the second-difference operator at the boundary points
results in $$\label{eq:D2Neum2}
{\bf D}^2_{n2} = \frac{1}{h^2}\begin{bmatrix}
{\bf -2}& 2&   &  &  &  \\ 
-1& {\bf -2} & 1 &  &  &\\
& \ddots & \ddots & \ddots &  &      \\
& &  1 & {\bf -2}& 1   \\
& &   & 2 & {\bf -2}   \\
\end{bmatrix},$$ and the matrix is $M+1 \times M+1$. Note that this
matrix is *not* symmetric, and it is therefore not possible to write it
as the product of a matrix times its transpose. Note, however, that the
matrix can be symmetrised by applying a linear transformation in the
form of a positive-definite diagonal matrix. To that end, consider the
$M+1\times M+1$ diagonal matrix
${\bf T} = \text{diag}([\frac{1}{2},1,1,...,1,1,\frac{1}{2}])$. Then
$$\label{eq:LinTransfn2}
{\bf T}\,{\bf D}^{2}_{n2} = \frac{1}{h^2}\begin{bmatrix}
{\bf -1}& 1&   &  &  &  \\ 
-1& {\bf -2} & 1 &  &  &\\
& \ddots & \ddots & \ddots &  &      \\
& &  1 & {\bf -2}& 1   \\
& &   & 1 & {\bf -1}   \\
\end{bmatrix},$$ that is,
[\[eq:D2Neum1\]](#eq:D2Neum1){reference-type="eqref"
reference="eq:D2Neum1"}, except the dimensions of the matrix are here
$M+1 \times M+1$.

### Eigenvalues and eigenvectors of second difference matrices

The analysis is Section
[1.1.1](#sec:FDtransformationsSpatial){reference-type="ref"
reference="sec:FDtransformationsSpatial"} revealed that an eigenfunction
of the $\delta_{xx}$ operator is $e^{-jm h \gamma}$, with eigenvalue
$-\frac{4}{h^2}\sin^2\left( \frac{\gamma h}{2}\right)$. Of course, a
second eigenfunction is given by $e^{jm h \gamma}$, so that the
eigenfunctions of the second difference matrix can be written as
$${\bf D}^2 \hat {\bf y} = -\frac{4}{h^2}\sin^2\left( \frac{\gamma h}{2}\right)\hat {\bf y},$$
where $\hat {\bf y} = A_+ \hat {\bf y}^+ + A_- \hat {\bf y}^-$ and where
$\hat { y}^+_m = e^{jmh\gamma}$, $\hat { y}^-_m = e^{-jmh\gamma}$. Here,
${\bf D}^2$ is either one of
[\[eq:D2Diri\]](#eq:D2Diri){reference-type="eqref"
reference="eq:D2Diri"},
[\[eq:D2Neum1\]](#eq:D2Neum1){reference-type="eqref"
reference="eq:D2Neum1"} or
[\[eq:D2Neum2\]](#eq:D2Neum2){reference-type="eqref"
reference="eq:D2Neum2"}. In practice, the quantised wavenumbers $\gamma$
and the corresponding eigenvectors are obtained after application of the
boundary conditions.

#### Eigenvalues of ${\bf D}^2_d$

Application of ${\hat{y}}_{m=0}={\hat{y}}_{m=M}=0$ results in the
following system of equations for $A_+$, $A_-$ $$\label{eq:eigenD2d}
\begin{bmatrix}
1 & 1 \\
e^{jMh\gamma} & e^{-jMh\gamma}
\end{bmatrix}
\begin{bmatrix}
A_+ \\
A_-
\end{bmatrix}=
\begin{bmatrix}
0 \\
0
\end{bmatrix},$$ which has non-trivial solutions if and only if the
determinant is set to zero, i.e. when $$\label{eq:waveNrD2d}
\sin ({Mh\gamma}) = 0  \implies  \gamma_p = \frac{p \pi}{M h}, \,\, p = 1,...,M-1.$$
The eigenvectors are obtained by solving either one of
[\[eq:eigenD2d\]](#eq:eigenD2d){reference-type="eqref"
reference="eq:eigenD2d"}. This gives $A_+=-A_-$, and thus
$$\label{eq:EigenVecsD2d}
{\hat y}^p_m = A_+ \sin \frac{mp\pi}{M },\,\, m = 1,..,M-1.$$ In the
above $\hat{ y}^p_m$ denotes the $m^{th}$ component of the $p^{th}$
eigenvector. The constant $A_+$ can be set to normalise the eigenvectors
conveniently, as shown below. The eigenvectors can be ordered columnwise
in a matrix,
$${\bf Y}_d = \begin{bmatrix}\hat {\bf y}^1 & \hat{\bf y}^2 & ... & \hat{\bf y}^{M-1}\end{bmatrix},$$
so that, from the spectral theorem, a decomposition of the
second-difference matrix is obtained as $$\label{eq:DecompD2d}
{\bf D}^{2}_d = -{\bf Y}_d \, {\bf \Lambda}_d \, {\bf Y}^\intercal_d,$$
where ${{\bf \Lambda}_d}$ is a diagonal matrix, whose diagonal elements
$\lambda_p$ are
$\frac{4}{h^2}\sin^2 \left( \frac{\gamma_p h}{2}\right)$. Here,
$$\label{eq:BndLambdaD2d}
0 < \lambda_p < 4/h^2.$$ Note that, by virtue of the bounds on the
eigenvalues, this matrix is clearly positive-definite. When the
eigenvectors are normalised using $A_+ = \sqrt{2/M}$, the matrix
${\bf Y}_d$ is orthonormal, and thus
${\bf Y}_d {\bf Y}_d^{\intercal}={\bf Y}_d^\intercal {\bf Y}_d = {\bf I}$.

#### Eigenvalues of ${\bf D}^2_{n1}$

Application of $\delta_{x+}\hat y_0 = \delta_{x-}\hat y_M = 0$ gives,
after a little algebra, the following system of equations
$$\label{eq:eigenD2n1}
\begin{bmatrix}
e^{jh\gamma} & -1 \\
e^{j(M-1)h\gamma} & -e^{-jMh\gamma}
\end{bmatrix}
\begin{bmatrix}
A_+ \\
A_-
\end{bmatrix}=
\begin{bmatrix}
0 \\
0
\end{bmatrix},$$ Setting the determinant to zero gives
$$\label{eq:EigenvaluesD2n1}
\sin ({(M-1)h\gamma}) = 0  \implies  \gamma_p = \frac{p \pi}{(M-1) h}, \,\, p = 0,...,M-2.$$
Note that the range for $p$ is now shifted, so to include the zero
eigenvalue. This must be included: looking at
[\[eq:D2Neum1\]](#eq:D2Neum1){reference-type="eqref"
reference="eq:D2Neum1"}, the vector $\hat {\bf y} = [1,1,...,1]$ is
necessarily an eigenvector, with eigenvalue equal to zero. The
eigenvectors are obtanained solving either one of
[\[eq:eigenD2n1\]](#eq:eigenD2n1){reference-type="eqref"
reference="eq:eigenD2n1"}. This gives
$$\hat y^p_m=A_+\cos\frac{p \pi \left( m-\frac{1}{2}\right)}{M-1}, \,\, m = 1,..,M-1.$$
As before, the eigenvectors can be arranged in columns,
$${\bf Y}_{n1} = \begin{bmatrix}\hat {\bf y}^1 & \hat{\bf y}^2 & ... & \hat{\bf y}^{M-1}\end{bmatrix},$$
so that, from the spectral theorem, a decomposition of the
second-difference matrix is obtained as $$\label{eq:DecompD2n1}
{\bf D}^{2}_{n1} = -{\bf Y}_{n1} \, {\bf \Lambda}_{n1} \, {\bf Y}^\intercal_{n1}.$$
where ${{\bf \Lambda}_{n1}}$ is a diagonal matrix, whose diagonal
elements $\lambda_p$ are equal to
$\frac{4}{h^2}\sin^2 \left( \frac{\gamma_p h}{2}\right)$. Here,
$$\label{eq:BndLambdaD2n1}
0 \leq \lambda_p < 4/h^2.$$ When the eigenvectors are normalised using
$A_+ = \sqrt{2/(M-1)}$, the matrix ${\bf Y}_{n1}$ is orthonormal, and
thus
${\bf Y}_{n1} {\bf Y}_{n1}^{\intercal}={\bf Y}_{n1}^\intercal {\bf Y}_{n1} = {\bf I}$.

#### Eigenvalues of ${\bf D}^2_{n2}$

Here, applying the boundary conditions
$\delta_{x\cdot}y_0 = \delta_{x\cdot}y_M = 0$ gives the following set of
equations $$\label{eq:eigenD2n2}
\begin{bmatrix}
1 & -1 \\
e^{jMh\gamma} & -e^{-jMh\gamma}
\end{bmatrix}
\begin{bmatrix}
A_+ \\
A_-
\end{bmatrix}=
\begin{bmatrix}
0 \\
0
\end{bmatrix},$$ and, upon setting the determinant to zero, the
quantised wavenumbers are such that $$\label{eq:EigenvaluesD2n2}
\sin ({Mh\gamma}) = 0  \implies  \gamma_p = \frac{p \pi}{M h}, \,\, p = 0,...,M.$$
These are the same as
[\[eq:waveNrD2d\]](#eq:waveNrD2d){reference-type="eqref"
reference="eq:waveNrD2d"}, except for the range of $p$, which again must
include the zero eigenvalue since the vector
$\hat {\bf y}=[1,1,..,1]^\intercal$ is necessarily an eigenvector. The
eigenvectors are then given by solving either one of
[\[eq:eigenD2n2\]](#eq:eigenD2n2){reference-type="eqref"
reference="eq:eigenD2n2"}, yielding
$$\hat y^p_m=A_+\cos\frac{m p \pi}{M}, \,\, m = 0,..,M.$$ The
eigenvectors can be ordered columnwise in a matrix,
$${\bf Y}_{n2} = \begin{bmatrix}\hat {\bf y}^1 & \hat{\bf y}^2 & ... & \hat{\bf y}^{M-1}\end{bmatrix},$$
so that, from the spectral theorem, a decomposition of the
second-difference matrix is obtained as $$\label{eq:DecompD2n2}
{\bf D}^{2}_{n2} = -{\bf Y}_{n2} \, {\bf \Lambda}_{n2} \, {\bf Y}^{-1}_{n2},$$
where ${{\bf \Lambda}_{n2}}$ is a diagonal matrix, whose diagonal
elements $\lambda_p$ are equal to
$\frac{4}{h^2}\sin^2 \left( \frac{\gamma_p h}{2}\right)$. Here,
$$\label{eq:BndLambdaD2n2}
0 \leq \lambda_p \leq 4/h^2.$$ Note that now, since ${\bf D}^2_{n2}$ is
*not* symmetric, the matrix ${\bf Y}_{n2}$ is not orthogonal.

## An explicit finite difference scheme

As a working finite difference scheme discretising
[\[eq:WE\]](#eq:WE){reference-type="eqref" reference="eq:WE"}, consider
$$\label{eq:WEFDexp}
\delta_{tt}{\bf y}^n = c^2 {\bf D}^2 {\bf y}^n.$$ Here, we are not
bothered with the particular form of the second-difference matrix. It
may be stem from application of either Neumann or Dirichlet boundary
conditions, so the subscript is left blank. Expanding out the second
time difference in [\[eq:WEFDexp\]](#eq:WEFDexp){reference-type="eqref"
reference="eq:WEFDexp"}, one gets
$${\bf y}^{n+1} = \left(2{\bf I}+c^2k^2{\bf D}^2 \right){\bf y}^n - {\bf y}^{n-1},$$
giving an explicit update. While simple, the properties of this scheme
are still unknown, particularly with respect to stability and accuracy.
Just like the choice of the system's parameters poses an upper bound on
the choice of the time step in the case of the harmonic oscillator, here
the choice of the time step and the grid spacing is limited by a
stability condition.

### Stability via Energy Analysis

An idea of stability from distributed systems is encapsulated in the
continuous energy balance given in Section
[\[sec:BoundWECnt\]](#sec:BoundWECnt){reference-type="ref"
reference="sec:BoundWECnt"}. A similar idea can be employed in the
discrete case. First, write
$${\bf D}^2 = -{\bf Y} \, {\bf \Lambda} \, {\bf Y}^{-1},$$ for a
positive-semidefinite diagonal matrix ${\bf \Lambda}$, which holds true
for the Dirichlet, first- and second-order Neumann conditions, as in
[\[eq:DecompD2d\]](#eq:DecompD2d){reference-type="eqref"
reference="eq:DecompD2d"},
[\[eq:DecompD2n1\]](#eq:DecompD2n1){reference-type="eqref"
reference="eq:DecompD2n1"},
[\[eq:DecompD2n2\]](#eq:DecompD2n2){reference-type="eqref"
reference="eq:DecompD2n2"}. Then, multiply
[\[eq:WEFDexp\]](#eq:WEFDexp){reference-type="eqref"
reference="eq:WEFDexp"} on the left by ${\bf Y}^{-1}$, to get
$$\delta_{tt}{\bf q}^n = -c^2 {\bf \Lambda} {\bf q}^n, \,\, \text{where }{\bf q}^n = {\bf Y}^{-1}{\bf y}^n.$$
Multiplying the equation above on the left by
$\delta_{t\cdot}{\bf q}^\intercal$, and using the usual identities, one
gets $$\delta_{t+}{\mathfrak h}^{n-1/2} = 0,$$ where
$$\label{eq:EnWaveEqnFD}
{\mathfrak h}^{n-1/2} = \frac{(\delta_{t-}{\bf q}^n)^\intercal \delta_{t-}{\bf q}^n}{2} + \frac{c^2(e_{t-}{\bf q}^n)^\intercal {\bf \Lambda} {\bf q}^n}{2}.$$
This expresses an energy balance in terms of the modified state variable
$\bf q$. While conserved, the discrete energy is not necessarily
positive, because the potential term is of indefinite sign. It may be
useful to re-write the discrete energy as a quadratic form for the
vector
${\bf p} = \frac{1}{\sqrt{2}}[({\bf q}^n)^\intercal, ({\bf q}^{n-1})^\intercal]^\intercal$.
This gives
$${\mathfrak h}^{n-1/2} = {\bf p}^\intercal \begin{bmatrix}\frac{\bf I}{k^2} & -\frac{\bf I}{k^2}+\frac{c^2}{2}{\bf \Lambda} \\ -\frac{\bf I}{k^2}+\frac{c^2}{2}{\bf \Lambda} & \frac{\bf I}{k^2} \end{bmatrix}{\bf p}.$$
Positive-definiteness may be established after inspection of the Schur
complement, that is, the matrix is positive definite if and only if its
Schur complement is:
$$\frac{\bf I}{k^2} - k^2 \left( -\frac{\bf I}{k^2}+\frac{c^2}{2}{\bf \Lambda}\right)^2 = -\frac{c^2k^2}{2}{\bf \Lambda}\left(-\frac{2{\bf I}}{k^2} + \frac{c^2}{2}{\bf \Lambda} \right) \geq 0.$$
Now, since we established that ${\bf \Lambda}$ is positive definite,
with largest eigenvalue $4/h^2$, the condition above results in
establishing when
$-\frac{2{\bf I}}{k^2} + \frac{c^2}{2}{\bf \Lambda}\leq 0$, that is:
$$\label{eq:CFL}
h\geq ck.$$ The condition above takes the name of
*Courant-Friedrichs-Lewy*, or CFL, condition. It is a stability
condition for scheme
[\[eq:WEFDexp\]](#eq:WEFDexp){reference-type="eqref"
reference="eq:WEFDexp"}. In practice, given a wave velocity $c$ and a
time step $k$, the grid spacing cannot be chosen to be smaller than
bound [\[eq:CFL\]](#eq:CFL){reference-type="eqref" reference="eq:CFL"}.

While [\[eq:EnWaveEqnFD\]](#eq:EnWaveEqnFD){reference-type="eqref"
reference="eq:EnWaveEqnFD"} expresses a discrete energy balance, it is
not immediate to see how this is related to the continuous energy
balance [\[eq:EnBalFull\]](#eq:EnBalFull){reference-type="eqref"
reference="eq:EnBalFull"}. This difficulty is only apparent: a natural
discretisation of
[\[eq:EnBalFull\]](#eq:EnBalFull){reference-type="eqref"
reference="eq:EnBalFull"} is already encoded in
[\[eq:EnWaveEqnFD\]](#eq:EnWaveEqnFD){reference-type="eqref"
reference="eq:EnWaveEqnFD"}, when the physical state ${\bf y}^n$ is
substituted back in. First, start with the Dirichlet and non-centred
Neumann conditions. For these, one has
${\bf Y}^\intercal = {\bf Y}^{-1}$, and thus
[\[eq:EnWaveEqnFD\]](#eq:EnWaveEqnFD){reference-type="eqref"
reference="eq:EnWaveEqnFD"} can be expressed as
$${\mathfrak h}^{n-1/2} = \frac{(\delta_{t-}{\bf y}^n)^\intercal \delta_{t-}{\bf y}^n}{2} + \frac{c^2(e_{t-}{\bf {\bf D}^-{\bf y}}^n)^\intercal \, {\bf D}^-{\bf y}^n}{2},$$
where the fact that ${\bf D}^2 = {\bf D}^+{\bf D}^-$ was used, as per
[\[eq:SecondDiffMat\]](#eq:SecondDiffMat){reference-type="eqref"
reference="eq:SecondDiffMat"} and
[\[eq:SecondDiffMatNeumann\]](#eq:SecondDiffMatNeumann){reference-type="eqref"
reference="eq:SecondDiffMatNeumann"}. For centered conditions, however,
the second-difference matrix is not symmetric, and ${\bf Y}$ is not
orthogonal, meaning that
[\[eq:EnWaveEqnFD\]](#eq:EnWaveEqnFD){reference-type="eqref"
reference="eq:EnWaveEqnFD"} cannot be directly converted to a discrete
counterpart of [\[eq:EnBalFull\]](#eq:EnBalFull){reference-type="eqref"
reference="eq:EnBalFull"}. To solve this issue, one may employ the
linear transformation
[\[eq:LinTransfn2\]](#eq:LinTransfn2){reference-type="eqref"
reference="eq:LinTransfn2"} on the wave equation
[\[eq:WEFDexp\]](#eq:WEFDexp){reference-type="eqref"
reference="eq:WEFDexp"}, to get
$$\delta_{tt}{\bf T}\,{\bf y}^n = c^2 {\bf T}\,{\bf D}^2 \,{\bf y}^n,$$
Remember that now the matrix ${\bf T}\,{\bf D}^2$ is symmetric, and can
be decomposed as ${\bf D}^2 = {\bf D}^+{\bf D}^-$. Multiplying the
equation above on the left by ${\bf y}^\intercal$ gives the energy
$${\mathfrak h}^{n-1/2} = \frac{(\delta_{t-}{\bf y}^n)^\intercal {\bf T} \delta_{t-}{\bf y}^n}{2} + \frac{c^2(e_{t-}{\bf {\bf D}^-{\bf y}}^n)^\intercal \, {\bf D}^-{\bf y}^n}{2},$$
and of course the stability analysis of this modified energy gives the
same CFL condition as [\[eq:CFL\]](#eq:CFL){reference-type="eqref"
reference="eq:CFL"}.

### Accuracy

The question of accuracy for scheme
[\[eq:WEFDexp\]](#eq:WEFDexp){reference-type="eqref"
reference="eq:WEFDexp"} is an interesting one. First, it may be useful
to derive the *numerical dispersion relation*. Remember that, for the
continuous system, the dispersion relation is given by
[\[eq:DispRelWECnt\]](#eq:DispRelWECnt){reference-type="eqref"
reference="eq:DispRelWECnt"}, in which the temporal frequencies are
proportional to the wavenumbers, with constant of proportionality given
by the wave speed $c$. Now, re-write scheme
[\[eq:WEFDexp\]](#eq:WEFDexp){reference-type="eqref"
reference="eq:WEFDexp"} in operator form, and transform in the frequency
domain: $$\label{eq:WEfdOperator}
\delta_{tt}y^n_m = c^2 \delta_{xx}y^n_m \implies -\frac{4}{k^2}\sin^2 \left( \frac{\omega k}{2}\right) = -\frac{4 c^2}{h^2}\sin^2 \left( \frac{\gamma h}{2}\right).{}$$
When the stability condition
[\[eq:CFL\]](#eq:CFL){reference-type="eqref" reference="eq:CFL"} is
satisfied with equality, the numerical dispersion relation above gives
$$\label{eq:NumDispRelWE}
\omega = c \gamma,$$ that is, the same as the continuous system! In
practice, the system is dispersionless. This is certainly too strong a
statement, since the dispersion relation above was obtained in the case
of an infinite grid, without boundaries. For all practical applications,
the grid is finite and the wavenumbers and frequencies are quantised. It
is the accuracy of these quantised frequencies and wavenumbers that must
be checked. From the previous section, we know that the continuous
boundary conditions can be discretised in a number of ways, leading to
different expressions for the quantised wavenumbers. First, consider the
numerical Dirichlet conditions. The wavenumbers are given by
[\[eq:waveNrD2d\]](#eq:waveNrD2d){reference-type="eqref"
reference="eq:waveNrD2d"}. When one chooses $h$ to be exactly equal to
$L/M$, then the numerical wavenumbers are the same as the continuous
ones, up to $m=M-1$, and hence the numerical frequencies are obtained
from [\[eq:NumDispRelWE\]](#eq:NumDispRelWE){reference-type="eqref"
reference="eq:NumDispRelWE"}, as
$$\omega_p = \frac{pc\pi}{Mh}, \,\, p = 1,...,M-1,$$ and these are the
same as the eigenfrequencies of the continuous system,
[\[eq:omeM\]](#eq:omeM){reference-type="eqref" reference="eq:omeM"}.
Note that, of course, the eigenvectors of the ${\bf D}^2_{d}$ matrix,
given in [\[eq:EigenVecsD2d\]](#eq:EigenVecsD2d){reference-type="eqref"
reference="eq:EigenVecsD2d"} are a sampled version of the eigenvectors
of the continuous system. So, for Dirchlet conditions, satisfying the
stability condition [\[eq:CFL\]](#eq:CFL){reference-type="eqref"
reference="eq:CFL"} and using the matrix ${\bf D}^2_{d}$ given in
[\[eq:D2Diri\]](#eq:D2Diri){reference-type="eqref"
reference="eq:D2Diri"} does indeed result in an exact scheme in the
frequency domain! In the time domain, things are slightly more
complicated, due to the discretisation of the initial conditions,
discussed below. This is nonetheless a powerful result, and indeed one
rarely found in other systems.

For Neumann conditions, the first-order accurate discretisation
${\bf D}^2_{n1}$ gives the quantised wavenumbers
[\[eq:EigenvaluesD2n1\]](#eq:EigenvaluesD2n1){reference-type="eqref"
reference="eq:EigenvaluesD2n1"}. When one chooses $h = L/M$, the
quantised wavenumbers are not exactly the same as the continuous one, so
here an approximation is introduced. However, the eigenvalues of the
matrix ${\bf D}^2_{n2}$, in
[\[eq:EigenvaluesD2n2\]](#eq:EigenvaluesD2n2){reference-type="eqref"
reference="eq:EigenvaluesD2n2"}, are again exact, and the eigenvectors
are sampled versions of the continuous eigenvectors.

The order of accuracy of scheme
[\[eq:WEFDexp\]](#eq:WEFDexp){reference-type="eqref"
reference="eq:WEFDexp"} may as well be approached in terms of the Taylor
series of its operators, given in
[\[eq:WEfdOperator\]](#eq:WEfdOperator){reference-type="eqref"
reference="eq:WEfdOperator"}. We want to get the local truncation error
of the scheme. To that end, assume now that $y(t,x)$ is in fact the true
solution of the continuous wave equation
[\[eq:WE\]](#eq:WE){reference-type="eqref" reference="eq:WE"}. Then,
apply the difference operators to $y(t,x)$. This gives the definition of
the local truncation error (LTE), as follows:
$$\tau^n_m \triangleq (\delta_{tt}- c^2 \delta_{xx})y(t=kn,x=mh).$$
Taylor-expanding the operators, one gets
$$\left(\frac{\partial^2}{\partial t^2} - c^2 \frac{\partial^2}{\partial x^2} + \frac{1}{12}\left(k^2 \frac{\partial^4}{\partial t^4} - c^2 h^2 \frac{\partial^4}{\partial x^4}\right) + \frac{1}{360}\left(k^4 \frac{\partial^6}{\partial t^6} - c^2 h^4 \frac{\partial^6}{\partial x^6}\right) + ... \right)y(t,x) = \tau^n_m.$$
When the stability condition
[\[eq:CFL\]](#eq:CFL){reference-type="eqref" reference="eq:CFL"} is
satisfied with equality, each term in the expansion contains the factor
$\frac{\partial^2}{\partial t^2} - c^2 \frac{\partial^2}{\partial x^2}$,
so that a factorisation of the above results as
$$\left(1 + \frac{k^2}{12}\left(\frac{\partial^2}{\partial t^2} + c^2 \frac{\partial^2}{\partial x^2}\right) + \frac{k^4}{360}\left( \frac{\partial^4}{\partial t^4} + c^4 \frac{\partial^4}{\partial x^4} +c^2 \frac{\partial^4}{\partial x^2 \partial t^2} \right) + ... \right)\left(\frac{\partial^2}{\partial t^2} - c^2 \frac{\partial^2}{\partial x^2}\right) y(t,x) = \tau^n_m.$$
But since $y(t,x)$ solves [\[eq:WE\]](#eq:WE){reference-type="eqref"
reference="eq:WE"}, then one gets $\tau^n_m=0$, that is, the LTE is
zero!

### Initialisation. Global error. 

The LTE computed above is exactly zero when $m$ is an inner point in the
domain, and $n \geq 2$ is some time step away from the initial time
steps. The global error of the finite difference scheme is defined as
$$\label{eq:ErrWaveEqn}
E^n_m = y(kn,mh)-y^n_m,$$ that is, the difference between the exact
solution computed at the time $t=kn$, $x=mh$, and the output of the
difference scheme at the corresponding time step and grid point.
Clearly, the global error includes the errors made in approximating the
boundary and initial conditions. The former, as seen, can be
approximated exactly via the matrices ${\bf D}^2_{d}$, ${\bf D}^2_{n2}$.
As for the latter, a suitable discretisation of the continuous initial
conditions [\[eq:ICs\]](#eq:ICs){reference-type="eqref"
reference="eq:ICs"} must be given, for some smooth initial displacement
$y_0(x)$ and velocity $v_0(x)$. As per usual, one may set
$$\label{eq:y0m}
y^0_m = y_0(x=mh),$$ that is, one may simply sample the continuous
initial shape to set the value of the grid function ${\bf y}^0$. The
grid function ${\bf y}^1$ is obtained approximating the initial
condition on the velocity. In analogy with
[\[eq:HigherOrderICsOscillator\]](#eq:HigherOrderICsOscillator){reference-type="eqref"
reference="eq:HigherOrderICsOscillator"}, obtained for the oscillator,
one may higher-order initial conditions may be obtained here by
expanding out the second time difference operator $\delta_{tt}$.

::: subequations
[\[eq:HigherOrderICsWE\]]{#eq:HigherOrderICsWE
label="eq:HigherOrderICsWE"} $$\begin{aligned}
\delta_{t+}y^0_m &= v_0(x=mh) \quad \text{first order} \label{eq:HigherOrderICsWE1}\\
\left(\delta_{t+}-\frac{k}{2}\delta_{tt}\right) y^0_m &= v_0(x=mh) \quad \text{second order} \label{eq:HigherOrderICsWE2}\\
\left(\delta_{t+}-\frac{k}{2}\delta_{tt}- \frac{k^2}{6}\delta_{t+}\delta_{tt}\right) y^0_m &= v_0(x=mh) \quad \text{third order}\label{eq:HigherOrderICsWE3}\\
\left(\delta_{t+}-\frac{k}{2}\delta_{tt}- \frac{k^2}{6}\delta_{t+}\delta_{tt}- \frac{k^3}{24}\delta_{tt}^2 \right) y^0_m &= v_0(x=mh) \quad \text{fourth order}\label{eq:HigherOrderICsWE4}\end{aligned}$$
:::

In the expressions above, substituting
$\delta_{tt}y^0_m = c^2 \delta_{xx}y^0_m$ gives a way to compute
$y^1_m$, knowing $y^0_m$ and $v_0$. While this idea is simple in theory,
things turn out to be a little more complicated than this. To figure out
what is going on, it may be useful in fact to plot the error curves as a
function of $k$, the time step. In order to do so, we are going to
compute the numerical error of the scheme against the exact solution.
The initial conditions are given as follows $$y_0(x)  = 
\left\{ 
\begin{array}{ll}
1 - \cos \left(\frac{2\pi x}{L/2}\right)& \text{if }0\leq x \leq L/2, \\
0 & \text{elsewhere,}
\end{array}
\right.$$ and $v_0(x)=0$. Under such intial conditions, an exact
solution is given by the D'Alembert soluition
[\[eq:Tr\]](#eq:Tr){reference-type="eqref" reference="eq:Tr"}. Boundary
conditions of fixed type are chosen, so that the wave reflects with a
change of sign at the boundary. In practice, a discrete counterpart of
[\[eq:Tr\]](#eq:Tr){reference-type="eqref" reference="eq:Tr"} is needed,
in the form of a digital waveguide. For that, divide the length of the
string in $M$ subintervals. Define $h=L/M$, $k = h/c$. Define the
left-travelling wave front as $l_{m}^n$, and the right-travelling wave
front as $r^n_m$. For $n=1$, one has
$$l^1_m = r^1_m = \frac{y^0_m}{2},\quad m=1,...,M-1.$$ where $y^0_m$ is
defined in [\[eq:y0m\]](#eq:y0m){reference-type="eqref"
reference="eq:y0m"}. Then, at the timestep $n\geq 2$, the travelling
wavefronts are updated as follows
$$l^n_m = l^{n-1}_{m+1}  \,\, (m < M-1), \quad
l^n_{M-1} = r^{n-2}_{M-1}, \quad
r^n_m = r^{n-1}_{m-1}, \,\, (m > 1), \quad
r^n_1 = l^{n-2}_{1}.$$ This discretises the wave equation exactly. Then
the ouput of the finite difference scheme, with initial conditions as
per
[\[eq:HigherOrderICsWE\]](#eq:HigherOrderICsWE){reference-type="eqref"
reference="eq:HigherOrderICsWE"}, is then compared against the output of
the exact solution, and the error
[\[eq:ErrWaveEqn\]](#eq:ErrWaveEqn){reference-type="eqref"
reference="eq:ErrWaveEqn"} is stored for various timesteps $k$. It is
convenient to use an even number $M$ of subintervals, and to check the
output of the wave equation at $m=M/2 +1$ (the central point of the
grid). The reference sample is `Ts=floor(te/k)`, where $t_e$ is the
duration in seconds of the simulation. The results are presented in Fig.
[1.1](#fig:ErrorCurvesWE){reference-type="ref"
reference="fig:ErrorCurvesWE"}.

![Error curves for scheme
[\[eq:WEFDexp\]](#eq:WEFDexp){reference-type="eqref"
reference="eq:WEFDexp"}, under various initial conditions. (a):
[\[eq:HigherOrderICsWE1\]](#eq:HigherOrderICsWE1){reference-type="eqref"
reference="eq:HigherOrderICsWE1"}; (b):
[\[eq:HigherOrderICsWE2\]](#eq:HigherOrderICsWE2){reference-type="eqref"
reference="eq:HigherOrderICsWE2"}; (c):
[\[eq:HigherOrderICsWE3\]](#eq:HigherOrderICsWE3){reference-type="eqref"
reference="eq:HigherOrderICsWE3"}; (d):
[\[eq:HigherOrderICsWE4\]](#eq:HigherOrderICsWE4){reference-type="eqref"
reference="eq:HigherOrderICsWE4"}. The error is computed as per
[\[eq:ErrWaveEqn\]](#eq:ErrWaveEqn){reference-type="eqref"
reference="eq:ErrWaveEqn"}, where $m = M/2 + 1$ (M is always chosen to
be even); and $n = \text{floor}(0.0165/k)$. Then, $h=L/M$, and
$k = h/c$, where $L=1$, and
$c=315$.](Figures/ErrorWEeqn.eps){#fig:ErrorCurvesWE
width="\\linewidth"}

The error plots are revealing: for
[\[eq:HigherOrderICsWE1\]](#eq:HigherOrderICsWE1){reference-type="eqref"
reference="eq:HigherOrderICsWE1"} and
[\[eq:HigherOrderICsWE3\]](#eq:HigherOrderICsWE3){reference-type="eqref"
reference="eq:HigherOrderICsWE3"}, the error curves attain the expected
trends, but for
[\[eq:HigherOrderICsWE2\]](#eq:HigherOrderICsWE2){reference-type="eqref"
reference="eq:HigherOrderICsWE2"} and
[\[eq:HigherOrderICsWE4\]](#eq:HigherOrderICsWE4){reference-type="eqref"
reference="eq:HigherOrderICsWE4"}, they do not! In fact,
[\[eq:HigherOrderICsWE2\]](#eq:HigherOrderICsWE2){reference-type="eqref"
reference="eq:HigherOrderICsWE2"} reveals that the error is down to
machine accuracy: in other words, the initialisation is *exact* in this
case. This is understood by a simple arguments: it is immediate to
verify that, under the condition that $h=ck$,
[\[eq:HigherOrderICsWE2\]](#eq:HigherOrderICsWE2){reference-type="eqref"
reference="eq:HigherOrderICsWE2"} gives
$y^1_m = \frac{y^0_{m+1}+y^0_{m-1}}{2}$, that is, the discrete
D'Alembert solution. This holds true for inner points, but, as we know
from the previous section, the boundary conditions are also discretised
exactly in this case via the matrix ${\bf D}^2_d$, and hence the scheme
overall is exact!

The interpretation of Fig.
[1.1](#fig:ErrorCurvesWE){reference-type="ref"
reference="fig:ErrorCurvesWE"}(d) is a little more obscure. Note that,
in order to implement the corresponding initial condition
[\[eq:HigherOrderICsWE4\]](#eq:HigherOrderICsWE4){reference-type="eqref"
reference="eq:HigherOrderICsWE4"}, one needs to apply the discrete
difference $\delta_{xx}$ twice or, equivalently, form a matrix
${\bf D}^4_d$. Here, the choice was to use the product of ${\bf D}^2_d$
with itself: ${\bf D}^4_d = {\bf D}^2_d \, {\bf D}^2_d$. This is the
problem: ${\bf D}^4_d$ is not of sufficiently high accuracy to increase
the accuracy of the scheme overall, which hence reamins third-order. A
higher-order discretisation of $\delta_{xx}^2$ needs to be implemented,
though this possibility will not be explored further here.

## Loss and forcing

The wave equation including loss and source terms can be given as
$$\frac{\partial^2 y}{\partial t^2} = c^2 \frac{\partial^2 y}{\partial x^2}-2 \sigma \frac{\partial y}{\partial t} + p(x,t).$$
Here, $\sigma \geq 0$ (measured in s$^{-1}$) is a viscous-loss, as in
[\[eq:WEloss\]](#eq:WEloss){reference-type="eqref"
reference="eq:WEloss"}, and $p$ (measured in m/s$^2$) is a source term.
For all practical purposes, one may set
$p(x,t) = \frac{\eta(x)}{\rho A} f(t)$, that is, the source is separable
in $x$ and $t$, with $f$ being a force signal measured in N. The energy
balance for this equation, an extension of
[\[eq:EnBalWE\]](#eq:EnBalWE){reference-type="eqref"
reference="eq:EnBalWE"}, reads $$\label{eq:EnBalanceWEQNloss}
\frac{d}{dt}\left( \int_{0}^{L} \frac{1}{2}\left( \frac{\partial y}{\partial t} \right)^2 dx + \int_{0}^{L} \frac{c^2}{2}\left( \frac{\partial y}{\partial x} \right)^2 dx \right) = -2\sigma \int_0^L \left(\frac{\partial y}{\partial t}\right)^2 \, dx + \frac{f(t)}{\rho A}\int_0^L \eta \frac{\partial y}{\partial t}\, dx.$$
Here, conservative boundary conditions are assumed, for instance of
fixed type, so that the corresponding boundary integrals vanish in the
energy balance.

### Finite difference schemes

An approximation using finite differences is obtained immediately after
defining the time-dependent vector ${\bf y}^n$, approximating the
solution $y(t,x)$ on the grid. Then, take
$$\delta_{tt}{\bf y}^n = c^2 {\bf D}^2_d {\bf y}^n - 2\sigma \delta_{t\cdot}{\bf y}^n + \frac{f^n}{\rho A}{\boldsymbol \eta}.$$
Here, the grid function ${\boldsymbol \eta}$ is a discrete approximation
to the continuous distribution $\eta(x)$, and $f^n$ is a time series
sampling the input $f(t)$ at the current time step. Regardless of the
particular discretisation for $\eta$, notice that the discrete energy
balance is obtained immediately after multiplication of the left by
$\delta_{t\cdot}{\bf y}^\intercal$, giving
$$\delta_{t+}\left(\frac{\delta_{t-}({\bf y}^n)^\intercal \, \delta_{t-}{\bf y}^n}{2} + \frac{c^2(e_{t-}{\bf {\bf D}^-{\bf y}}^n)^\intercal \, {\bf D}^-{\bf y}^n}{2} \right) = -2\sigma (\delta_{t\cdot}{\bf y}^n)^\intercal \delta_{t\cdot}{\bf y}^n +  \frac{f^n}{\rho A}(\delta_{t\cdot}{\bf y}^n)^\intercal {\boldsymbol \eta},$$
that is a discrete counterpart of
[\[eq:EnBalanceWEQNloss\]](#eq:EnBalanceWEQNloss){reference-type="eqref"
reference="eq:EnBalanceWEQNloss"}.

### Spreading and interpolation

In many practical applications, the distribution of the load is
concentrated in a small area around the loading point $x_p$, so that
$\eta(x) \approx \delta(x-x_p)$. In this case, $\boldsymbol \eta$ is
obtained using a spreading operator of sufficient accuracy. Spreading
and interpolation are closely-related concepts, in that the former
consists in yielding a grid function, starting from the knowledge of a
given continuous function; the former instead produces a continuous
function, starting from the sampled values of a grid function. There are
many different kinds of interpolants. One such popular ones is due to
Lagrange, and it uses a basis of polynomials. Suppose we want to
interpolate the value of a point $x_p$, such that
$$m_p = \text{floor}(x_p/h), \quad \alpha = x_p/h - m_p,$$ In practice
$0\leq m_p\leq M$ is the grid point to the left of $x_p$, and
$0\leq \alpha \leq h$ represent the remainder of the flooring operation.
Then, consider the following arrays ${\bf r}$ of coefficients, of length
$M+1$: $$\label{eq:LagrangeArrays}
\begin{bmatrix}
\vdots\\
r_{mp} \\
\vdots
\end{bmatrix} = \frac{1}{h}\begin{bmatrix}
\vdots\\
1 \\
\vdots
\end{bmatrix},\quad 
\begin{bmatrix}
\vdots\\
r_{mp-1} \\
r_{mp} \\
r_{mp+1} \\
\vdots
\end{bmatrix} = \frac{1}{h}\begin{bmatrix}
\vdots\\
\frac{1-\alpha}{2} \\
0 \\
\frac{1+\alpha}{2} \\
\vdots
\end{bmatrix},\quad
\begin{bmatrix}
\vdots\\
r_{mp-1} \\
r_{mp} \\
r_{mp+1} \\
\vdots
\end{bmatrix} = \frac{1}{h}\begin{bmatrix}
\vdots\\
\frac{\alpha(\alpha-1)}{2} \\
(1+\alpha)(1-\alpha) \\
\frac{\alpha(\alpha+1)}{2} \\
\vdots
\end{bmatrix},\quad
\begin{bmatrix}
\vdots\\
r_{mp-1} \\
r_{mp} \\
r_{mp+1} \\
r_{mp+2} \\
\vdots
\end{bmatrix} = \frac{1}{h}\begin{bmatrix}
\vdots\\
-\frac{\alpha(\alpha-1)(\alpha-2)}{6} \\
\frac{(\alpha+1)(\alpha-1)(\alpha-2)}{2} \\
-\frac{\alpha(\alpha+1)(\alpha-2)}{2} \\
\frac{\alpha(\alpha+1)(\alpha-1)}{6}\\
\vdots
\end{bmatrix}$$

![Errors of the discrete spreading operators
[\[eq:LagrangeArrays\]](#eq:LagrangeArrays){reference-type="eqref"
reference="eq:LagrangeArrays"}. Here, the error $E$, as per
[\[eq:ErrSpreading\]](#eq:ErrSpreading){reference-type="eqref"
reference="eq:ErrSpreading"}, is computed starting from the continuous
function $y(x) = \sin (\pi x)$, at $x_p=0.289$, and where the grid
function ${\bf y}$ contains the sampled values of $y$ at $x=mh$, where
$h=1/M$. The plots are obtained using a number of values for $M$. The
panels, from $(a)$ to $(d)$, include higher-order approximations to the
Delta function, as given in
[\[eq:LagrangeArrays\]](#eq:LagrangeArrays){reference-type="eqref"
reference="eq:LagrangeArrays"}.](Figures/ErrorsSpreading.eps){#fig:ErrsLagrange
width="\\linewidth"}

Lagrange interpolation is simply given by
$$y(x_p) = h{\bf r}^\intercal {\bf y},$$ where $\bf r$ is any one of the
arrays given in
[\[eq:LagrangeArrays\]](#eq:LagrangeArrays){reference-type="eqref"
reference="eq:LagrangeArrays"}. Of course, the arrays yields
coefficients such that the order of the approximation increases. The
error is defined simply as $$\label{eq:ErrSpreading}
E = y(x_p)-h{\bf r}^\intercal {\bf y}.$$ In particular,
[\[eq:LagrangeArrays\]](#eq:LagrangeArrays){reference-type="eqref"
reference="eq:LagrangeArrays"} gives approximations of order one to
four. Spreading, on the other hand, is somewhat the reversed operation.
Note that the equation above can be written in an equivalent form as
$$y(x_p) \triangleq \int_0^L y(x) \delta(x-x_p)\, dx = h{\bf r}^\intercal {\bf y},$$
in pratice yielding an approximation to the Dirac delta function as the
array $\bf r$. In Fig. [1.2](#fig:ErrsLagrange){reference-type="ref"
reference="fig:ErrsLagrange"}, the accuracy of the spreading operator is
checked in a numerical experiment, yielding the expected error trends.
