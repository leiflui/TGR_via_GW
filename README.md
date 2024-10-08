# Test of General Relativity via Gravitational-Wave Energy Deviations
The following repository is a compilation of codes used to test general relativity (GR) using the energy deviations from gravitational radiation. The majority of the calculations here are done in Python and compiled in a Jupyter notebook. Moreover, we are primarily using the PyCBC software package to generate the waveforms and modify the post-Newtonian (PN) coefficients. 

GitHub repository for PyCBC: https://github.com/gwastro/pycbc
## Energy Dissipated by Gravitational Waves of Coalescing binaries
Similar to electromagnetic waves, we can compute the energy of a gravitational wave (GW) using multipole expansion. To do so, we first compute the Isaacson stress-energy tensor
```math
    t_{\mu\nu}=-\frac{1}{8\pi}\left< R_{\mu\nu}^{(2)}-\frac{1}{2}\bar{g}_{\mu\nu} R^{(2)} \right>,
```
where $R_{\mu\nu}^{(2)}$ is the Ricci tensor to the expanded to second order in the perturbed metric $h_{\mu\nu}$, $`\bar{g}_{\mu\nu}`$ is the background metric, and the angle brackets $`\left< \dots \right>`$ denotes the average over several wavelengths. Basically, the Isaacson stress-energy tensor is the Einstein equation expanded to second order in $h_{\mu\nu}$. $R_{\mu\nu}^{(2)}$ usually involves many terms quadratic in the metric perturbation, however, we can simplify this expression by performing integration by parts and using the transverse-traceless(TT) gauge condition. This is basically where we choose $h_{0\alpha}=h^i_{i}=\nabla_jh^{ij}=0$. Thus, 
```math
    R_{\mu\nu}^{(2)}=\left<\partial_{\mu}h_{\alpha\beta}^{TT}\partial_{\nu}h^{\alpha\beta}_{TT}\right>.
```
Therefore the Isaacson stress-energy tensor can be written explicitly as 
```math
    t_{\mu\nu}=\frac{1}{32\pi}\left<\partial_{\mu} h_{\alpha\beta}^{\text{TT}}\partial_{\nu}h^{\alpha\beta}_{\text{TT}}\right>.
```
To compute the energy carried by a GW, we take the 00-component of the Isaacson stress-energy tensor and integrate over the volume $V$ 
```math
    \frac{\mathrm{d} E}{\mathrm{d} t}=\lim_{r\to\infty} \frac{1}{16\pi}\int_S  \left<\dot{h}_{+}^2+\dot{h}_{\times}^2\right>r^2\mathrm{d} \Omega  ,
```
where $h_{+}$ and $h_{\times}$ are the plus and cross polarizations of the GW, respectively, and the overdots represent the derivative with respect to coordinate time. 
## Aligned-Spin Binaries and Spin-Weighted Spherical Harmonics
To obtain an analytic expression for the integral over solid angle $\Omega$, we separate $h_{+}$ and $h_{\times}$ into time-dependent and angular parts. This can be done using spin-weighted spherical harmonics $`_{s}Y_{\ell m}`$. For outgoing GWs, we are concerned with the spin weight $s=-2$, as the graviton is a massless spin-2 particle and the negative sign corresponds to outgoing GWs. Hence, the $+$- and $\times$-polarizations can be written as
```math
        h_{+}-ih_{\times}=\sum_{lm}{_{-2}}Y_{lm}(\theta,\phi)h_{lm}(t)\approx{_{-2}}Y_{22}(\theta,\phi)h_{22}(t)+{_{-2}}Y_{2-2}(\theta,\phi)h_{2-2}(t),
```
where the second line has been truncated because, for aligned-spin binaries, the $h_{22}$ and $h_{2-2}$ modes are the leading order terms. Taking the time derivative of the equation above and multiplying this expression by its complex conjugate, we obtain the following expression
```math
    \dot{h}_{+}^2+\dot{h}_{\times}^2=|{_{-2}}Y_{22}|^2|\dot{h}_{2,2}|^2+{_{-2}}Y_{22}\;{_{-2}}Y_{2-2}^{*}\dot{h}_{22}\dot{h}^{*}_{2-2}+{{_{-2}}}Y_{22}^{*}\;{_{-2}}Y_{2-2}\dot{h}^{*}_{22}\dot{h}_{2-2}+|{_{-2}}Y_{2-2}|^2|\dot{h}_{2-2}|^2.
```
With this, the task at hand is to calculate $h_{22}$ and $h_{2-2}$. To do so, we first need to determine the spin-weighted spherical harmonics. Using the Wigner D matrix, ${_{-2}}Y_{22}$ and ${_{2}}Y_{22}$ in the $\theta$ and $\phi$ representation is as follows
```math
    {_{-2}}Y_{22}(\theta,\phi)=\sqrt{\frac{5}{64\pi}}(1+\cos\theta)^2e^{2i\phi},
```
```math
    {_{-2}}Y_{2-2}(\theta,\phi)=\sqrt{\frac{5}{64\pi}}(1-\cos\theta)^2e^{-2i\phi}.
```
Notice that $`{_{-2}}Y_{22}(0,0)={_{-2}}Y_{2-2}(0,\pi)=\frac{1}{2}\sqrt{\frac{5}{\pi}}`$ and $`{_{-2}}Y_{2-2}(0,0)={_{-2}}Y_{22}(0,\pi)=0`$ . Therefore, to solve for $h_{22}$ and $h_{2-2}$, we simply calculate $h_{+}-ih_{\times}$ at $\theta=\phi=0$, and $\phi=0$, $\theta=\pi$. This can be done using PyCBC and LALSimulation. Doing so we find that
```math
        h_{22}(t)=\sqrt{\frac{4\pi}{5}}[h_{+}(t,0,0)-ih_{\times}(t,0,0)],
```
and
```math
        h_{2-2}(t)=\sqrt{\frac{4\pi}{5}}[h_{+}(t,\pi,0)-ih_{\times}(t,\pi,0)].
```
Since we are considering non-precessing binaries, $\theta$ contains no time dependence. Therefore, to integrate Eq. \eqref{multipole} over solid angle, we simply compute the following integrals,
```math
    \int_S |{_{-2}}Y_{22}|^2\mathrm{d}\Omega=\int_S |{_{-2}}Y_{2-2}|^2\mathrm{d}\Omega=1,
```
and 
```math
\int_S {_{-2}}Y_{22}^{*}\;{_{-2}}Y_{2-2}\mathrm{d}\Omega= \int_S {_{-2}}Y_{22}\;{_{-2}}Y_{2-2}^{*}\mathrm{d}\Omega=\frac{1}{6}.
```
Using the above results, we can integrate $`\dot{h}_{+}^2+\dot{h}_{\times}^2`$ over solid angle to obtain the instantaneous power
```math
\frac{\mathrm{d} E}{\mathrm{d} t}=\lim_{r\to\infty}\frac{r^2}{16\pi}\left<|\dot{h}_{22}|^2+|\dot{h}_{2-2}|^2+\frac{1}{6}\left(\dot{h}^{*}_{22}\dot{h}_{2-2}+\dot{h}_{22}\dot{h}^{*}_{2-2}\right)\right>.
```
To avoid averaging over several wavelengths, we calculate the total energy by numerically integrating $\frac{\mathrm{d} E}{\mathrm{d} t}$ over the time array.
## Dephasing Coefficients
```math
    p_i\to (1+\delta p_i)p_i.
```
These fractional deviations are known as the dephasing coefficients. The phasing of IMRPhenomPv2 consists of three regimes. The first of which is the inspiral regime which is parameterized by PN coefficients [A. Buonanno et al., (2013).] $`\left\{\chi_0,\dots,\chi_7 \right\}`$ and $`\left\{\chi_{5l},\chi_{6l}\right\}`$. In this regime, there are also phenomenological parameters $`\left\{\sigma_0,\dots,\sigma_4\right\}`$ that contribute to the high effective PN order. This corrects for non-adiabaticity in the late inspiral phase and unknown high-order PN coefficients in the adiabatic regime. The second regime, is the intermediate regime, which is parameterized by the phenomenological coefficients $`\left\{\beta_0,\dots,\beta_3\right\}`$. Finally, there is the merger-ringdown regime which is parameterized by a combination of the phenomenological coefficients and the analytical black-hole perturbation theory parameters $`\left\{\alpha_0,\dots,\alpha_5\right\}`$ [J. Meidam et at., (2018).]. As one can see, if $\delta p_i=0$ this corresponds to a theory with no deviation from GR. 
## Truncating the Waveform in the Inspiral Regime
One potential obstacle is that the upper bound and lower bound for the integration of time was kind of chosen arbitrarily. Since we evaluate the energy over the same hypersurface, the total energy emitted by a GW should be the same, it is just that the energy is distributed over different time intervals. Therefore, if we integrate from $`-\infty < t < \infty`$, the $`\delta E/E_{GR}`$ should always be zero. This means we have to choose a fair and physically meaningful time interval. Otherwise, $\delta E/E_{GR}$ will obviously be different as the energy is distributed differently over time, but it will not be physically meaningful. Recall that $\delta\chi_i$ only changes the waveform in the inspiral phase, therefore we have to find a way to cut off the GW waveform before it reaches the merger phase. To do this, recall that the inspiral phase is characterized by $0.004\leq Mf\leq0.018$, where $M$ is the total mass and $f$ is the GW frequency. For this, we can constrain the upper bound of $t/M$ (Here, $t/M=0$ corresponds to the time of BBH coalescence.) by using the upper bound cutoff frequency $Mf=0.018$ by the following equation
```math
\frac{\mathrm{d} t}{\mathrm{d} f}=F(M,q,\delta\chi,f)
```
where $F(M,q,\delta\chi,f)$ is some function of $M$, $q$, $\delta\chi$, and $f$. There is no analytical form of $F(M,q,\delta\chi,f)$, so we need to compute the integral numerically. Similarly, we can set the lower bound to be $Mf=0.004$ and use the same formula to calculate the corresponding lower bound for $t/M$. This way, we can include most of the inspiral phase in our energy calculation. 

![Alt text](examples/relative_velocity.png)

## Interpolating the Gravitational-Wave Energy Deviation
As seen previously, each energy deviation is calculated through numerical integration and is quite computationally expensive to do this for many values of $\delta\chi_i$. Therefore, one can interpolate between points to obtain $`\delta E/E_{GR}`$ without large computational cost. We show a figure below where we vary $\delta\chi_4$, $M$, and $q$.
![Alt text](examples/4D_colormap_0p004=Mf=0p018_30by30by30_interpolant-1.png)
To make sure that the interpolant is giving numerically accurate results, we the interpolated results to that of the ones actually calculated by numerical integration.
![Alt text](examples/interpolant_check.png)
## Finding the Dephasing Coefficients that Yield the Largest Energy Deviations
Here, we use the interpolant shown previously to determine plot the contour and the vector field where changing the dephasing coefficient yields the largest energy deviation from GR. As an example, we plot the 3D contour along with the vector field of the energy deviation for varying $\delta \chi_4$ and $\delta \chi_{6}$.
![Alt text](examples/Varying_dchi4_and_dchi_6_corrected.png)
We show another example of a 3D vector field when varying the 3 dephasing coefficients. Here the arrows have been omitted to avoid clutter.
![Alt text](examples/4D_colormap_0p004=Mf=0p018_25by25_interpolant_varying_dchi0_dchi1_and_dchi7_q=1_M=60.png)
