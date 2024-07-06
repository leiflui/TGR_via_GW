# Test of General Relativity via Gravitational-Wave Energy Deviations
The following repository is a compilation of codes used to test general relativity (GR) using the energy deviations from gravitational radiation. The majority of the calculations here are done in Python and compiled in a Jupyter notebook. Moreover, we are primarily using the PyCBC software package to generate the waveforms and modify the post-Newtonian (PN) coefficients. 

GitHub repository for PyCBC: https://github.com/gwastro/pycbc
## Energy Dissipated by Gravitational Waves of Coalescing binaries
Similar to electromagnetic waves, we can compute the energy of a gravitational wave (GW) using multipole expansion. To do so, we first compute the Isaacson stress-energy tensor
```math
    t_{\mu\nu}=-\frac{1}{8\pi}\left<R_{\mu\nu}^{(2)}-\frac{1}{2}\bar{g}_{\mu\nu} R^{(2)}\right>,
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
    \frac{\dd E}{\dd t}=\lim_{r\to\infty} \frac{1}{16\pi}\int_S  \left<\dot{h}_{+}^2+\dot{h}_{\times}^2\right>r^2\dd \Omega  ,
```
where $h_{+}$ and $h_{\times}$ are the plus and cross polarizations of the GW, respectively, and the overdots represent the derivative with respect to coordinate time. 
