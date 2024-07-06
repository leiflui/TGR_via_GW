# Test of General Relativity via Gravitational-Wave Energy Deviations
The following repository is a compilation of codes used to test general relativity (GR) using the energy deviations from gravitational radiation. The majority of the calculations here are done in Python and compiled in a Jupyter notebook. Moreover, we are primarily using the PyCBC software package to generate the waveforms and modify the post-Newtonian (PN) coefficients. 

GitHub repository for PyCBC: https://github.com/gwastro/pycbc
## Energy Dissipated by Gravitational Waves of Coalescing binaries
Similar to electromagnetic waves, we can compute the energy of a gravitational wave (GW) using multipole expansion. To do so, we first compute the Isaacson stress-energy tensor
```math
    t_{\mu\nu}=-\frac{1}{8\pi}\left< R_{\mu\nu}^{(2)}-\frac{1}{2}\bar{g}_{\mu\nu} R^{(2)} \right>,
```
where $R_{\mu\nu}^{(2)}$ is the Ricci tensor to the expanded to second order in the perturbed metric $h_{\mu\nu}$, $\bar{g}_{\mu\nu}$ is the background metric, and the angle brackets $\left<\dots\right>$ denotes the average over several wavelengths. Basically, the Isaacson stress-energy tensor is the Einstein equation expanded to second order in $h_{\mu\nu}$. $R_{\mu\nu}^{(2)}$ usually involves many terms quadratic in the metric perturbation, however, we can simplify this expression by performing integration by parts and using the transverse-traceless(TT) gauge condition. This is basically where we choose $h_{0\alpha}=h^i_{i}=\nabla_jh^{ij}=0$. Thus, 
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
To obtain an analytic expression for the integral over solid angle $\Omega$, we separate $h_{+}$ and $h_{\times}$ into time-dependent and angular parts. This can be done using spin-weighted spherical harmonics $_{s}Y_{\ell m}$. For outgoing GWs, we are concerned with the spin weight $s=-2$, as the graviton is a massless spin-2 particle and the negative sign corresponds to outgoing GWs. Hence, the $+$- and $\times$-polarizations can be written as
```math
        h_{+}-ih_{\times}=\sum_{lm}{_{-2}}Y_{lm}(\theta,\phi)h_{lm}(t)\approx{_{-2}}Y_{22}(\theta,\phi)h_{22}(t)+{_{-2}}Y_{2-2}(\theta,\phi)h_{2-2}(t),
```
where the second line has been truncated because, for aligned-spin binaries, the $h_{22}$ and $h_{2-2}$ modes are the leading order terms. Taking the time derivative of the equation above and multiplying this expression by its complex conjugate, we obtain the following expression
```math
    \dot{h}_{+}^2+\dot{h}_{\times}^2=|{_{-2}}Y_{22}|^2|\dot{h}_{2,2}|^2+{_{-2}}Y_{22}\;{_{-2}}Y_{2-2}^{*}\dot{h}_{22}\dot{h}^{*}_{2-2}+{{_{-2}}}Y_{22}^{*}\;{_{-2}}Y_{2-2}\dot{h}^{*}_{22}\dot{h}_{2-2}+|{_{-2}}Y_{2-2}|^2|\dot{h}_{2-2}|^2.
```
