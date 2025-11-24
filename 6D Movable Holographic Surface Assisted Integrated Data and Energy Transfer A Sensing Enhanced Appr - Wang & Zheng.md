# 6D Movable Holographic Surface Assisted Integrated Data and Energy Transfer: A Sensing Enhanced Approach

Zhonglun Wang, Graduate Student Member, IEEE, Yizhe Zhao, Member, IEEE, Gangming Hu, Yali Zheng, Member, IEEE, and Kun Yang, Fellow, IEEE

Abstract—Reconfigurable holographic surface (RHS) enables cost-effective large-scale arrays with high spatial gain. However, its amplitude-controlled holographic beamforming suffers from directional fluctuations, making it difficult to fully exploit the spatial gain of RHS. Fortunately, the promising 6D movable antenna (6DMA) provides a potential solution to this problem. In this paper, we study a 6D movable holographic surface (6DMHS) integrated data and energy transfer (IDET) system, where a three-stage protocol is proposed, consisting of an uplink sensing stage, an orientation adjustment stage and a downlink transmission stage, to coordinate the 6DMHS and effectively serve the IDET receivers. Firstly, the holographic-based sensing technology is proposed and the sensing information of the IDET receivers is exploited. Secondly, by fixing the rotations with the sensing information, the orientation optimization problem is formulated for designing the holographic beamforming of the RHS and adjusting the translations of the 6DMHS. As a result, the directions with maximum beamforming gain are aligned with each IDET receiver. Thirdly, by fixing the orientation of the 6DMHS and the holographic beamforming, the equivalent wireless channel is obtained. The IDET performance optimization problem is formulated for obtaining the optimal digital beamforming, power splitting factor and energy harvesting (EH) power. Simulation results demonstrate that the proposed scheme is capable of improving the IDET performance compared to the benchmarks.

Index Terms—Reconfigurable holographic surface (RHS), 6D movable holographic surface (6DMHS), integrated data and energy transfer (IDET), holographic-based sensing.

# I. INTRODUCTION

# A. Backgrounds

In the 6G era, the number of the Internet of Things (IoT) devices is expected to increase exponentially to support the Internet of Everything [1]–[3]. Traditional power supply methods, such as batteries and wired grids, inevitably incur significant maintenance costs. To overcome this problem, the integrated data and energy transfer (IDET) technology has been proposed, enabling simultaneous wireless power supply and data transmission to the IoT devices. Moreover, the IDET receivers is anticipated to be widely deployed in the large-scale

Zhonglun Wang, Yizhe Zhao and Gangming Hu are with the School of Information and Communication Engineering, University of Electronic Science and Technology of China, Chengdu 611731, China, e-mail: 202021010524@std.uestc.edu.cn; yzzhao@uestc.edu.cn; 202521010916@std.uestc.edu.cn

Yali Zheng and Kun Yang are with the School of Computer Science and Electronic Engineering, University of Essex, Essex CO4 3SQ, U.K., e-mail: kunyang@essex.ac.uk.

IoT networks, thereby significantly enhancing the network efficiency while reducing the maintenance overhead [4].

One of the most effective ways to enhance IDET performance is to employ massive phased arrays and design beamforming schemes that improve wireless transmission efficiency [5]. However, the phased arrays exhibits certain inherent limitations. Specifically, the phased arrays require a lot of phase shifters, RF-chains and power amplifier, which may lead to a higher cost. Additionally, the complicated wiring and large physical size make it difficult to employ the massive phased array. Consequently, the cost and size may restrict the development of the traditional phased array [6].

Fortunately, the reconfigurable holographic surface (RHS), composed of the printed circuit board (PCB), DC control circuits, and diodes, has emerged as a promising transmit antenna technology for realizing massive MIMO and highly accurate beamforming [7]. In particular, the RHS employs an amplitude-controlled holographic beamformer, where the holographic principle is leveraged to generate directional beamforming for wireless transmission. However, due to the absence of the imaginary component in holographic beamforming, the beamforming gain of the RHS exhibits directional variations, i.e., the gains in different directions are not identical. This motivates us to rotate the RHS such that its maximum beamforming direction aligns with the IDET receiver, thereby fully exploiting its beamforming gain and improving IDET performance. Moreover, the emerging six-dimensional movable antenna (6DMA), containing of three rotation dimensions and three translation dimensions, offers a promising opportunity to realize this scheme [8]. Hence, by integrating the RHS with the 6DMA, we propose a novel six-dimensional movable holographic surface (6DMHS) to exhaustively exploit the beamforming gain of the RHS and enhance IDET performance.

# B. Related Works

Extensive research efforts have been devoted to enhancing IDET performance. The related studies on IDET and 6DMHS are summarized as follows.

Initially, the concept of the IDET is proposed in [9]. Then, many methods are adopts for improving the IDET performance. Bruno et al. [10] proposed the dedicated OFDM-based IDET waveform to compensate the nonlinear energy harvesting (EH) circuit, thereby improving the IDET performance.

Moreover, Kwon et al. [5] investigated the massive MIMO-assisted IDET transmitter and designing the beamforming to improve the IDET performance. Furthermore, to increase the coverage performance, Na et al. [11] studied the reconfigurable intelligent surface (RIS)-assisted IDET system and proved the existence of the RIS can improve the IDET performance effectively.

The aforementioned works mainly adopted the phased array-assisted transmitter. To overcome the costs and size drawbacks of the traditional phased array, the RHS-assisted IDET transmitter was proposed. The RHS was initially adopted for wireless communication. Deng eq al. [7] proposed the dedicated holographic beamforming for wireless communication, proving that the RHS can effectively improve the throughput performance. Then, for simplifying the design of the holographic beamforming, they [12] proposed the holographic-pattern division multiple access (HDMA) technology, which outperformed the traditional space-division multiple access (SDMA) in terms of throughput. Moreover, to further improve the throughput, Di et al. [13] adopted the RHS-assisted transmitter in the wideband OFDM system and designed the holographic to against the beamsquaint effect, proving that RHS outperformed the traditional phased array in the wideband system. For further improving the IDET performance, the RHS-assisted transmitter was also adopted in the IDET systems. Azarbahram et al. [14] jointly optimized the waveform and beamforming design in the RHS assisted wireless energy transfer (WET) system, demonstrating superior WET performance compared to the conventional phased arrays. Additionally, Huang et al. [15] optimized the holographic beamforming by considering the continuous-aperture RHS and practical electromagnetic-based wireless channels, thus attaining the performance limits achievable in the IDET systems.

Simultaneously, the aforementioned works adopted fixed-position antenna (FPA) to achieve IDET, which may not exploit the DoFs exhaustively. The movable antenna was proposed to address this problem. Firstly, the fluid antenna was proposed to exploit the higher spatial DoFs. Wong et al. [16] firstly proposed the fluid antenna-assisted communication systems, proving the effectiveness of this antenna. Then, Lin et al. [17] proposed the fluid antenna-assisted IDET systems, proving the improvement of the multiplexing gain for the fluid antenna. Note that the fluid antenna can only move in 2D space and the utilization of the DoFs is limited. To further improve the utilization of the DoFs, the 6DMA was proposed [18]. Shao et al. proposed the 6DMA-assisted system and designed the orientation of the 6DMA carefully, proving that the 6DMA was capable of achieving high DoFs and improving the IDET performance significantly. Moreover, by considering the finite discrete rotations and positions, they [8] designed the optimal orientation and power allocation scheme for improving the wireless communication performance. Furthermore, Wang et al. [19] adopted the 6DMA in the IDET system, demonstrating that the 6DMA was capable of improving the IDET performance effectively.

# C. Motivations and Contributions

However, the related works still have several drawbacks, which are summarized as follows.

- For RHS-assisted IDET systems, none of the existing studies considered the directional property of holographic beamforming caused by the absence of its imaginary component, nor did they investigate its impact on IDET performance.  
- For 6DMA-assisted IDET systems, none of the existing works considered equipping the system with large-scale antennas or RHS, which may limit further performance improvements.  
- For the holographic beamforming and the orientation optimization problems, the algorithm complexity scales with the square of the number of antennas, the number of rotations and the number of translations. When the antenna array or the RHS is sufficiently large, this complexity becomes prohibitively high, making the deployment of 6DMA in practice very challenging.

Motivated by these limitations, this paper proposes the 6DMHS-assisted IDET system. By considering the directional property of holographic beamforming, the RHS is rotated such that its maximum beamforming direction aligns with the IDET receiver to pursue higher performance. To address the prohibitive complexity of joint orientation and holographic beamforming design, a sensing-assisted scheme is further developed, while the sensing information is adopted to rotate all RHSs and design the holographic beamforming directly. As a result, the complexity may reduce significantly. Finally, a three-stage protocol is proposed to fully exploit the potential of 6DMHS in IDET systems. Our main contributions are summarized as follows.

- A 6DMHS-assisted IDET system is developed, in which the RHS is equipped with holographic-based sensing modules capable of performing sensing operations. The extracted sensing information is utilized to adjust the spatial orientation of the 6DMHS, thereby aligning the direction of the maximum beamforming gain of the RHS with the IDET receivers to improve the IDET performance.  
- To enable high-accuracy holographic sensing, an FFT-based detection method is adopted to extract the angular information of the receivers from the holographic image. Based on the sensing information, an alternating optimization algorithm is proposed to jointly optimize the spatial orientation of the 6DMHS and the holographic beamforming of the RHS. Subsequently, with the 6DMHS orientation fixed, the digital beamforming and power splitting factor are jointly optimized using fractional programming techniques.  
- Simulation results validate the superiority of the proposed holographic-based sensing design over conventional sensing method. The 6DMHS-assisted transmitter achieves significantly better IDET performance than systems relying on fixed-position RHS. Moreover, the proposed design outperforms the idealized perfect CSI-based scheme by

reducing pilot overhead while delivering improved overall performance.

# II. SYSTEM MODEL

In this section, we first introduce the model of the 6DMHS-assisted IDET system. Then, we introduce the wireless channel model of the system. Finally, we propose the WDT and WET model of the system.

# A. 6DMHS-Transmitter Model

As shown in Fig. 1(a), the 6DMHS-transmitter consists of a base station and  $B$  RHS. Each RHS is connected to the BS via a mechanical rod, along which a flexible cable is routed to transmit the IDET signals from the BS to the RHS. A translation motor and a rotary motor are installed at both ends of the mechanical rod to control its translation and rotation [18]. Moreover, each RHS is equipped with  $Q$  feeds and  $M = M_{x} \times M_{y}$  radiation elements, where  $M_{x}$  and  $M_{y}$  represent the number of the rows and columns of the RHS elements. The distance between two adjacent radiation elements is denoted as  $d$ . To fully utilize the beamforming gain of each RHS, we let each RHS serve one IDET receivers in a long frame. Therefore, there are  $K = B$  single-antenna IDET receivers in the system.

As shown in Fig. 1(b), the Rodrigues' rotation formula is adopted to formulate the rotation of the 6DMHS  $^1$  The  $b$ -th RHS is initially located on the  $xoy$  plane of the global coordinate system  $O - xyz$ , with its center at the origin  $O$ .  $\bar{\mathbf{u}}_b$  is the vector that represents the spatial information of the  $b$ -th RHS in the global coordinate system, such as its normal vector or the direction associated with the maximum beamforming gain. The vector  $\bar{\mathbf{n}}_b$  represents the outward normal vector of the  $b$ -th 6DMHS in the global coordinate system, which is aligned with the global  $z$ -axis. The coordinate of the  $(m_x,m_y)$ -th element for the  $b$ -th RHS in the global coordinate system is  $\bar{\mathbf{r}}_{m_x,m_y} = [m_xd,m_yd,0]^T$ . After translation and rotation operations, the center coordinate of the  $b$ -th RHS is  $q_b$ , while the normal vector of the  $b$ -th RHS is  $\mathbf{u}_b$ . The coordinate of the  $(m_x,m_y)$ -th element for the  $b$ -th RHS in the global coordinate system is given by

$$
\mathbf {r} _ {b, m _ {x}, m _ {y}} = \mathbf {q} _ {b} + \mathbf {R} _ {b} \bar {\mathbf {r}} _ {m _ {x}, m _ {y}}, \tag {1}
$$

where  $\mathbf{R}_b\in \mathbb{R}^{3\times 3}$  is the unitary rotation matrix and is expressed as

$$
\mathbf {R} _ {b} = \cos \alpha_ {b} \mathbf {I} _ {3 \times 3} + (1 - \cos \alpha_ {b}) \mathbf {v} _ {b} \mathbf {v} _ {b} ^ {T} + \sin \alpha_ {b} \mathbf {V} _ {b}. \tag {2}
$$

The matrix  $\mathbf{V}_b\in \mathbb{C}^{3\times 3}$  is given as

$$
\mathbf {V} _ {b} = \left[ \begin{array}{c c c} 0 & - \mathbf {v} _ {b} [ 3 ] & \mathbf {v} _ {b} [ 2 ] \\ \mathbf {v} _ {b} [ 3 ] & 0 & - \mathbf {v} [ 1 ] \\ \mathbf {v} _ {b} [ 2 ] & \mathbf {v} _ {b} [ 1 ] & 0 \end{array} \right]. \tag {3}
$$

<sup>1</sup>Note that the existing literature [18] adopts the Euler angles method, which consists three angles around three axes to describe the rotation.

Moreover,  $\mathbf{v}_b$  and  $\alpha_{b}$  represent the normalized rotation vector and the rotation angle of the  $b$ -th RHS, which can be expressed as

$$
\mathbf {v} _ {b} = \frac {\bar {\mathbf {u}} _ {b} \times \mathbf {u} _ {b}}{\| \bar {\mathbf {u}} _ {b} \times \mathbf {u} _ {b} \| _ {2}}, \tag {4}
$$

$$
\alpha_ {b} = \arccos  \frac {\left(\bar {\mathbf {u}} _ {b}\right) ^ {T} \mathbf {u} _ {b}}{\| \bar {\mathbf {u}} _ {b} \| _ {2} \| \mathbf {u} _ {b} \| _ {2}}, \tag {5}
$$

respectively. Additionally, in the local coordinate system  $O^{\prime} - x^{\prime}y^{\prime}z^{\prime}$ , the RHS lies on the  $x^{\prime}O^{\prime}y^{\prime}$  plane with its normal vector aligned with the  $z^{\prime}$ -axis, while the coordinate of the  $(m_x,m_y)$ -th element of the  $b$ -th RHS in this local coordinate system is  $\bar{\mathbf{r}}_{m_x,m_y}$ .

According to Eq. (2), we have

$$
\mathbf {r} _ {b, m _ {x}, m _ {y}} - \mathbf {q} _ {b} = m _ {x} d \mathbf {R} _ {b} \mathbf {e} _ {1} + m _ {y} d \mathbf {R} _ {b} \mathbf {e} _ {2}, \tag {6}
$$

where  $\mathbf{e}_1 = [1,0,0]^T$  and  $\mathbf{e}_2 = [0,1,0]^T$ .

Moreover, in order to ensure that the 6DMHS effectively radiates the IDET signal, three orientation constraints are introduced [20]:

- Rotation constraint for avoiding signal reflection:

$$
\mathbf {n} _ {b} ^ {T} \left(q _ {b _ {1}} - q _ {b _ {2}}\right) \leq 0, \forall b _ {1} \neq b _ {2}. \tag {7}
$$

- Rotation constraint for avoiding signal blockage:

$$
\mathbf {n} _ {b} ^ {T} \mathbf {q} _ {b} \geq 0, \forall b. \tag {8}
$$

- Minimum-distance constraint for avoiding collision:

$$
\left\| \mathbf {q} _ {b _ {1}} - \mathbf {q} _ {b _ {2}} \right\| _ {2} \geq d _ {\min }, \forall b _ {1} \neq b _ {2}. \tag {9}
$$

The vector  $\mathbf{n}_b$  denotes the outward normal of the  $b$ -th 6DMHS in the global coordinate system after translation and rotation operations.

# B. Channel Model

The wireless channel from the  $k$ -th IDET receiver to the  $b$ -th RHS is given as

$$
\mathbf {h} _ {k, b} = \sqrt {M} \sum_ {\iota} \Lambda_ {k, \iota , b} \eta_ {k, \iota , b} \mathbf {a} _ {b} \left(\theta_ {k, \iota}, \phi_ {k, \iota}\right), \tag {10}
$$

where  $\theta_{k,\iota}$  and  $\phi_{k,\iota}$  are the azimuth and elevation angles of the  $\iota$ -th scatterer for the  $k$ -th receiver, respectively.  $\eta_{k,\iota,b}$  is the channel complex channel gain from the  $k$ -th IDET receiver to the  $b$ -th RHS through the  $\iota$ -th scatterer, while  $\iota = 0$  represents the LoS link and  $\iota \geq 1$  represents the NLoS link.  $\eta_{k,\iota,b}$  can be further expressed as

$$
\eta_ {k, t, b} = \left\{ \begin{array}{l} \sqrt {\frac {K _ {\mathrm {R}}}{K _ {\mathrm {R}} + 1}} \sqrt {\Omega_ {k , b}}, \text {i f} \iota = 0, \\ \sqrt {\frac {1}{K _ {\mathrm {R}} + 1}} \sqrt {\Omega_ {k , b}} \varkappa_ {k, t, b}, \text {i f} \iota \geq 1, \end{array} \right. \tag {11}
$$

where  $K_{\mathbb{R}}$  is the Rician factor,  $\Omega_{k,b}$  is the path-loss from the  $k$ -th IDET receiver to the  $b$ -th RHS and  $\varkappa_{k,t,b} \sim CN(0,1)$  is the complex channel gain from the  $k$ -th IDET receiver to the  $b$ -th RHS through the  $\iota$ -th scatterer.  $\Lambda_{k,t,b}$  is the antenna gain from the  $k$ -th IDET receiver to the  $b$ -th RHS through the  $\iota$ -th

![](https://cdn-mineru.openxlab.org.cn/result/2025-11-24/7946ec08-51d9-4a1e-813e-6601cf747e28/d9c73d08d9364e3ec82f5801145efd8d44efa7f9ba98cd2e52d1ee682a5de52a.jpg)  
(a)

![](https://cdn-mineru.openxlab.org.cn/result/2025-11-24/7946ec08-51d9-4a1e-813e-6601cf747e28/6c8fe77844b387c0de6f2215ff9fbacfdfba795f0ba3f415584d4dae802f9e45.jpg)  
(b)  
Fig. 1. 6DMHS model: (a) 6DMHS-based transmitter; (b) Illustration of the coordinate system of the 6DMHS.

scatterer, which is given as [21]

$$
\Lambda_ {k, t, b} = \left\{ \begin{array}{l} 1, \text {i f} \mathbf {q} _ {b} ^ {T} \mathbf {f} \left(\theta_ {k, t}, \phi_ {k, t}\right) > 0, \\ 0, \text {e l s e .} \end{array} \right. \tag {12}
$$

where the direction vector  $\mathbf{f}(\theta, \phi)$  is defined as  $\mathbf{f}(\theta, \phi) = [f_1(\theta, \phi), f_2(\theta, \phi), f_3(\theta, \phi)]^T \in \mathbb{R}^{3 \times 1}$ , where we have  $f_1(\theta, \phi) = \cos \theta \cos \phi$ ,  $f_2(\theta, \phi) = \sin \theta \cos \phi$  and  $f_3(\theta, \phi) = \sin \phi$ . The steering vector  $\mathbf{a}_b(\theta_{k,\iota}, \phi_{k,\iota}) \in \mathbb{C}^{M \times 1}$  is defined as

$$
\begin{array}{l} \mathbf {a} _ {b} (\theta , \phi) = \sqrt {\frac {1}{M}} \left[ e ^ {j \frac {2 \pi}{\lambda} \mathbf {f} ^ {T} (\theta , \phi) \mathbf {r} _ {b, 0, 0}}, \dots , e ^ {j \frac {2 \pi}{\lambda} \mathbf {f} ^ {T} (\theta , \phi) \mathbf {r} _ {b, 0, M _ {y} - 1}}, \dots , \right. \\ \left. e ^ {j \frac {2 \pi}{\lambda} \mathbf {f} ^ {T} (\theta , \phi) \mathbf {r} _ {b, M _ {x} - 1, 0}}, \dots , e ^ {j \frac {2 \pi}{\lambda} \mathbf {f} ^ {T} (\theta , \phi) \mathbf {r} _ {b, M _ {x} - 1, M _ {y} - 1}} \right] ^ {T}. \tag {13} \\ \end{array}
$$

Define the wireless channel from the 6DMHS-assisted transmitter to the  $k$ -th IDET receiver as  $\mathbf{h}_k = \left[\mathbf{h}_{k,1}^T,\dots ,\mathbf{h}_{k,B}^T\right]^T\in$ $\mathbb{C}^{BM\times 1}$ .

# C. Downlink Transmission Signal Model

The downlink receive signal of the  $k$ -th IDET receiver is expressed as

$$
y _ {k} = \sum_ {b} \sum_ {k ^ {\prime}} \mathbf {h} _ {k, b} ^ {T} \operatorname {d i a g} \left(\boldsymbol {\Psi} _ {b}\right) \boldsymbol {\Theta} _ {b} \mathbf {X} _ {k ^ {\prime}, b} s _ {k ^ {\prime}} + z _ {k}, \tag {14}
$$

where  $\mathbf{X}_{k,b}$  represents the beamforming vector for the  $k$ -th IDET receiver at the  $b$ -th RHS,  $s_k$  is the IDET signal for the  $k$ -th receiver satisfying  $s_k \sim CN(0,1)$  and  $z_k$  is the antenna noise for the  $k$ -th receiver satisfying  $z_k \sim CN(0,\sigma_0^2)$ .  $\Theta_b$  represents the electromagnetic response of the  $b$ -th RHS and is denoted

as  $\Theta_{b} = [\Theta_{b,0},\dots ,\Theta_{b,Q - 1}]\in \mathbb{C}^{M\times Q}$ , where the  $(m_xM_x + m_y)$ -th entry of  $\Theta_{b,q}\in \mathbb{C}^{M\times 1}$  is given as

$$
\Theta_ {b, q, m _ {x} M _ {y} + m _ {y}} = \sqrt {\eta} \cdot e ^ {- \alpha \| \mathbf {r} _ {b, m _ {x}, m _ {y}} - \mathbf {r} _ {b, q} \| _ {2}} \cdot e ^ {- j \frac {2 \pi f _ {c} \varrho}{c} \| \mathbf {r} _ {b, m _ {x}, m _ {y}} - \mathbf {r} _ {b, q} \| _ {2}}, \tag {15}
$$

where  $\eta$  is the efficiency factor,  $\alpha$  is the real-valued propagation attenuation factor,  $\varrho$  is the refractive factor of RHS and  $\mathbf{r}_{b,q} \in \mathbb{R}^{3 \times 1}$  is the coordinate of the  $q$ -th feed at the  $b$ -th RHS, respectively. Generally, the real attenuation term  $e^{-\alpha \| \mathbf{r}_{b,m_x,m_y} - \mathbf{r}_{b,q} \|_2}$  is neglected, since  $\alpha$  is a very small value [22].  $\Psi_b \in \mathbb{C}^{M \times 1}$  is the RHS beamforming vector for the  $b$ -th RHS.

By employing the power splitter, the receive signal at the  $k$ -th IDET receiver is divided into two parts: one for information decoding and another for energy harvesting, which are expressed as

$$
y _ {\mathrm {I D}, k} =
$$

$$
\sqrt {1 - \rho_ {k}} \sum_ {b} \sum_ {k ^ {\prime}} \mathbf {h} _ {k, b} ^ {T} \operatorname {d i a g} \left(\boldsymbol {\Psi} _ {b}\right) \boldsymbol {\Theta} _ {b} \mathbf {X} _ {k ^ {\prime}, b} s _ {k ^ {\prime}} + \sqrt {1 - \rho_ {k}} z _ {k} + z _ {\text {c o v}, k}, \tag {16}
$$

$$
y _ {\mathrm {E H}, k} = \sqrt {\rho_ {k}} \sum_ {b} \sum_ {k ^ {\prime}} \mathbf {h} _ {k, b} ^ {T} \operatorname {d i a g} \left(\Psi_ {b}\right) \boldsymbol {\Theta} _ {b} \mathbf {X} _ {k ^ {\prime}, b} s _ {k ^ {\prime}} + \sqrt {\rho_ {k}} z _ {k}, \tag {17}
$$

where  $\rho_{k}$  is the power splitting factor for the  $k$ -th receiver and  $z_{\mathrm{cov},k}$  is the passband-to-baseband conversion noise satisfying  $z_{\mathrm{cov},k} \sim \mathcal{CN}(0, \sigma_{\mathrm{cov}}^2)$  [23].

The SINR of the  $k$ -th IDET receiver is expressed as

$$
\begin{array}{l} \gamma_ {k} = \\ \frac {(1 - \rho_ {k}) \mathbb {E} \left[ \left| \sum_ {b} \mathbf {h} _ {k , b} ^ {T} \operatorname {d i a g} (\boldsymbol {\Psi} _ {b}) \boldsymbol {\Theta} _ {b} \mathbf {X} _ {k , b} s _ {k} \right| ^ {2} \right]}{(1 - \rho_ {k}) \mathbb {E} \left[ \left| \sum_ {b} \sum_ {k ^ {\prime} \neq k} \mathbf {h} _ {k , b} ^ {T} \operatorname {d i a g} (\boldsymbol {\Psi} _ {b}) \boldsymbol {\Theta} _ {b} \mathbf {X} _ {k ^ {\prime} , b} s _ {k ^ {\prime}} \right| ^ {2} \right] + (1 - \rho_ {k}) \sigma_ {0} ^ {2} + \sigma_ {\mathrm {c o v}} ^ {2}} \\ = \frac {\left(1 - \rho_ {k}\right) \left| \sum_ {b} \mathbf {h} _ {k , b} ^ {T} \operatorname {d i a g} \left(\boldsymbol {\Psi} _ {b}\right) \boldsymbol {\Theta} _ {b} \mathbf {X} _ {k , b} \right| ^ {2}}{\left(1 - \rho_ {k}\right) \sum_ {k ^ {\prime} \neq k} \left| \sum_ {b} \mathbf {h} _ {k , b} ^ {T} \operatorname {d i a g} \left(\boldsymbol {\Psi} _ {b}\right) \boldsymbol {\Theta} _ {b} \mathbf {X} _ {k ^ {\prime} , b} \right| ^ {2} + \left(1 - \rho_ {k}\right) \sigma_ {0} ^ {2} + \sigma_ {\mathrm {c o v}} ^ {2}} \tag {18} \\ \end{array}
$$

Furthermore, the downlink throughput of the  $k$ -th IDET receiver is expressed as

$$
R _ {k} = \log_ {2} \left(1 + \gamma_ {k}\right), [ \mathrm {b i t / s / H z} ]. \tag {19}
$$

The EH power in terms of the RF signal for the  $k$ -th IDET receiver is expressed as

$$
\begin{array}{l} P _ {\mathrm {E H}, k} = \rho_ {k} \mathbb {E} \left[ \left| \sum_ {b} \sum_ {k ^ {\prime}} \mathbf {h} _ {k, b} ^ {T} \operatorname {d i a g} \left(\boldsymbol {\Psi} _ {b}\right) \boldsymbol {\Theta} _ {b} \mathbf {X} _ {k ^ {\prime}, b} s _ {k ^ {\prime}} \right| ^ {2} + \sigma_ {0} ^ {2} \right] \\ = \rho_ {k} \sum_ {k ^ {\prime}} \left| \sum_ {b} \mathbf {h} _ {k, b} ^ {T} \operatorname {d i a g} \left(\boldsymbol {\Psi} _ {b}\right) \boldsymbol {\Theta} _ {b} \mathbf {X} _ {k ^ {\prime}, b} \right| ^ {2} + \rho_ {k} \sigma_ {0} ^ {2}. \tag {20} \\ \end{array}
$$

After passing through the RF-to-DC conversion circuit, the harvested DC power for the  $k$ -th IDET receiver is expressed as

$$
\Gamma \left(P _ {\mathrm {E H}, k}\right) = \max  \left\{\frac {E _ {\mathrm {m}}}{e ^ {- \xi E _ {0} + \nu}} \left(\frac {1 + e ^ {- \xi E _ {0} + \nu}}{1 + e ^ {- \xi P _ {\mathrm {E H} , k} + \nu}} - 1\right), 0 \right\}, \tag {21}
$$

where  $\xi$  and  $\nu$  are the constant circuit parameters.  $E_0$  and  $E_{\mathrm{m}}$  represent the activation and saturation power of the EH circuit, respectively.

# III. DIRECTIONAL BEAMFORMING GAIN OF RHS AND PROTOCOL

# DESIGN

In this section, we first introduce the directional beamforming gain property of the RHS. Then, in view of this phe

![](https://cdn-mineru.openxlab.org.cn/result/2025-11-24/7946ec08-51d9-4a1e-813e-6601cf747e28/61b77afdc5359fd28d84beb1a5b1857b0fff8cf184ab098887c422278dd850e5.jpg)  
Fig. 2. An illustration of the directional beamforming gain of the RHS.

nomenon, we propose the practical sensing-assisted protocol to adequately utilize the 6DMHS for IDET.

# A. Directional Beamforming Gain of RHS

For the traditional phased array, the holographic beamformer is a complex vector and its phase varies over  $[- \pi, \pi]$ . However, for the RHS, the holographic beamformer  $\Psi_{b}$  is a real vector, meaning that the imaginary part of the beamformer vanishes. Hence, the RHS may not adjust the beamforming perfectly. We give an example as follows.

Denote the beamforming gain for the  $q$ -th feed of the  $b$ -th RHS towards the direction  $\mathbf{f}(\theta, \phi)$  IDET receiver as  $g_b(\theta, \phi)$ , which is expressed as Eq. (23) [13], [24], [25], where  $D_{b,q,m_x,m_y}$  and  $\Psi_b$  are given as

$$
\begin{array}{l} D _ {b, q, m _ {x}, m _ {y}} = \left\| \mathbf {r} _ {b, m _ {x}, m _ {y}} - \mathbf {r} _ {b, q} \right\| _ {2}. \\ \Psi_ {b} = \sum_ {q} \omega_ {b, q} \frac {\mathcal {R} \left\{\sqrt {\frac {M}{\eta}} \operatorname {d i a g} \left(\boldsymbol {\Theta} _ {b , q}\right) \mathbf {a} _ {b} (\theta , \phi) \right\} + 1}{2}. \tag {22} \\ \end{array}
$$

Observe from Eq. (23) that the beamforming gain consists of Parts I-III. The Part I in Eq. (23) contributes the mainly beamforming gain toward the direction  $(\theta, \phi)$ . The Part II in Eq. (23) is the twin holographic beamforming gain towards the direction  $(- \theta, -\phi)$  and it is generated due to the real-valued RHS beamforming. The Part III in Eq. (23) is the interference beamforming gain towards the other directions. The Part II and Part III in Eq. (23) may lead to a directional beamforming gain, i.e., the beamforming gain is different at different directions. A demonstration of directional beamforming gain is shown in Fig. 2, where the maximum beamforming gain is achieved at the direction  $\mathbf{f}(\theta, \phi) = [0.175, -0.275, 0.945]^T$ . Consequently, if an IDET receiver is always at the maximum beamforming gain direction of the RHS, a higher IDET performance may be achieved. In the next section, we adopt the 6DMA technique to adjust the orientation of the RHS for achieving the aforementioned method and enhancing the IDET performance. As a result, we give the definition of the directional beamforming gain of RHS.

Definition 1 (Directional Beamforming Gain of RHS): For an RHS with amplitude-controlled holographic beamform-

ing, the beamforming gain is different for different direction  $(\theta ,\phi)$ . Moreover, there exists the direction  $(\theta^{*},\phi^{*})$ , which satisfies

$$
g _ {b} \left(\theta^ {*}, \phi^ {*}\right) = \max  _ {\theta , \phi , \omega_ {b} \in S _ {1}, \mathbf {r} _ {b} ^ {\text {f e e d}} \in S _ {2}} g _ {b} (\theta , \phi), \tag {24}
$$

where  $\omega_{b} = [\omega_{b,0},\dots ,\omega_{b,Q - 1}]^{T}$  and  $\mathbf{r}_b^{\mathrm{feed}} = [\mathbf{r}_{b,0},\dots ,\mathbf{r}_{b,Q - 1}]$ . Recall that the beamforming gain for the  $q$  -th feed of the  $b$  -th RHS towards the  $k$  -th IDET receiver  $g_{b}(\theta^{*},\phi^{*})$  is defined as Eq. (23). The sets  $S_{1}$  and  $S_{2}$  are defined as

$$
\left. \boldsymbol {S} _ {1} = \left\{\left(\boldsymbol {\omega} _ {b}\right): \mathbf {1} ^ {T} \boldsymbol {\omega} _ {b} = 1 \right\}, \right. \tag {25}
$$

$$
\left. \boldsymbol {S} _ {2} = \left\{\left(\mathbf {r} _ {b} ^ {\text {f e e d}}\right): \left\| \mathbf {r} _ {b, i} - \mathbf {r} _ {b, j} \right\| _ {2} ^ {2} \geq \epsilon_ {r}, \forall i \neq j \right\}, \right. \tag {26}
$$

respectively.  $\epsilon_r$  is the minimal distance of the adjacent feeds.

According to Eq. (24), we can obtain the optimal feed positions and direction by adopting the Aquila Optimizer [26].

Remark 1: Once the RHS is fabricated, the feed positions are fixed. Therefore, in the subsequent sections, we adopt the optimal feed positions and keep them unchanged. Moreover, the direction  $(\theta ,\phi)$  obtained from the Aquila Optimizer is defined as the maximum beamforming gain direction of the RHS, and is adopted when designing the orientation of the 6DMHS in the following part.

# B. Sensing-Enhanced 6DMHS Protocol

Denote the number of discrete positions of 6DMHS as  $N$ . Denote the number of discrete rotations in each discrete position as  $L$ . To adjust the 6DMHS accurately,  $N$  and  $L$  are generally large enough. When adjusting the orientation of the 6DMHS, the algorithm complexity of the gradient-based algorithm for adopting the optimal scheme is  $O(M^2 N^2 L^2)$ . For RHS, the number of elements is sufficiently large, which leading to a very high algorithm complexity. This indicates that it is difficult to employ the 6DMHS in practice. To address this problem, we first acquire the angle information by high accurate holographic-based sensing method, and then align the RHS with the IDET receiver in the maximal beamforming gain direction for achieving higher WET performance by using the directional beamforming gain property of RHS. The sensing-enhanced 6DMHS protocol is illustrated in Fig. 3 and is summarized as follows.

- Stage I (Uplink Sensing stage): All the RHS of 6DMHS are positioned at their initial locations. The IDET receivers sequentially transmit the sensing signals to the 6DMHS-assisted transmitter using the time division multiple access (TDMA) protocol. The BS then estimates the angles of the IDET receivers using the holographic-based sensing method.  
- Stage II (Orientation adjustment): The 6DMHS-assisted transmitter jointly adjusts the orientation of the RHS and their holographic beamforming based on the acquired sensing information, ensuring that the direction of the maximum beamforming gain of each RHS aligns with the corresponding sensing angle.  
- Stage III (Downlink transmission): At the beginning of each IDET frame, the BS estimates the equivalent CSI

$$
\begin{array}{l} g _ {b} (\theta , \phi) = \left| \sum_ {q} \mathbf {a} _ {b} (\theta , \phi) ^ {T} \operatorname {d i a g} \left(\boldsymbol {\Psi} _ {b}\right) \boldsymbol {\Theta} _ {b, q} \right| \\ = \frac {1}{\sqrt {M}} \left| \sum_ {q, q ^ {\prime}} \omega_ {b, q ^ {\prime}} \sum_ {m _ {x}, m _ {y}} e ^ {- j \frac {2 \pi}{\lambda} \left(\mathbf {f} ^ {T} (\theta , \phi) \mathbf {r} _ {b, m _ {x}, m _ {y}} - \varrho D _ {b, q, m _ {x}, m _ {y}}\right)} \cdot \frac {\mathcal {R} \left\{e ^ {j \frac {2 \pi}{\lambda} \left(\mathbf {f} ^ {T} (\theta , \phi) \mathbf {r} _ {b , m _ {x} , m _ {y}} - \varrho D _ {b , q ^ {\prime} , m _ {x} , m _ {y}}\right)} + 1 \right\}}{2} \right| \\ = \frac {1}{\sqrt {M}} \left| \underbrace {\sum_ {q , q ^ {\prime}} \frac {\omega_ {b , q ^ {\prime}}}{4} \sum_ {m _ {x} , m _ {y}} e ^ {j \frac {2 \pi}{\lambda} \varrho (D _ {b , q , m _ {x} , m _ {y}} - D _ {b , q ^ {\prime} , m _ {x} , m _ {y}})}} _ {\text {P a r t I}} + \underbrace {\sum_ {q , q ^ {\prime}} \frac {\omega_ {b , q ^ {\prime}}}{4} \sum_ {m _ {x} , m _ {y}} e ^ {- j \frac {2 \pi}{\lambda} \left(\mathbf {f} ^ {T} (\theta , \phi) \mathbf {r} _ {b , m _ {x} , m _ {y}} - \varrho D _ {b , q , m _ {x} , m _ {y}}\right)} . \frac {e ^ {- j \frac {2 \pi}{\lambda} \left(\mathbf {f} ^ {T} (\theta , \phi) \mathbf {r} _ {b , m _ {x} , m _ {y}} - \varrho D _ {b , q ^ {\prime} , m _ {x} , m _ {y}}\right)} + 2}{4}} _ {\text {P a r t I I}} \right. \\ + \underbrace {\sum_ {q , q ^ {\prime}} \omega_ {b , q ^ {\prime}} \sum_ {m _ {x} , m _ {y}} e ^ {- j \frac {2 \pi}{\lambda} \left(\mathbf {f} ^ {T} (\theta , \phi) \mathbf {r} _ {b , m _ {x} , m _ {y}} - \varrho D _ {b , q , m _ {x} , m _ {y}}\right)} \cdot \frac {\mathcal {R} \left\{e ^ {j \frac {2 \pi}{\lambda} \left(\mathbf {f} ^ {T} (\theta , \phi) \mathbf {r} _ {b , m _ {x} , m _ {y}} - \varrho D _ {b , q ^ {\prime} , m _ {x} , m _ {y}}\right)} + 1 \right\}}{2}} _ {\text {P a r t I I I}} \Bigg |. \tag {23} \\ \end{array}
$$

![](https://cdn-mineru.openxlab.org.cn/result/2025-11-24/7946ec08-51d9-4a1e-813e-6601cf747e28/eaf01d69eb426e5b26e45249cf48187f38ec71b19672301e7921ffa16d7104a6.jpg)  
Fig. 3. Protocol of the 6DMHS-assisted IDET system.

![](https://cdn-mineru.openxlab.org.cn/result/2025-11-24/7946ec08-51d9-4a1e-813e-6601cf747e28/e5e2672ba71937bda00c1c4b61b6909189b75c524b8906d0cb4888bea3b20128.jpg)

Channel estimation

![](https://cdn-mineru.openxlab.org.cn/result/2025-11-24/7946ec08-51d9-4a1e-813e-6601cf747e28/ca9fa5c137f2f65321e3ceb1645164db10d0e85ecea46edd1ea4e5227e5592b3.jpg)

IDET transmission

of each receiver at each RF chain. This equivalent CSI is defined as the cascade of the wireless channel, the holographic beamforming and the electromagnetic response of the RHS. Based on the estimated equivalent CSI, the digital beamforming and the power splitting strategy for each IDET receiver are then optimized. As a result, the IDET signals are transmitted to the IDET receivers successfully.

# IV. HOLEGROPHIC-BASED SENSING METHOD AND ORIENTATION ADJUSTMENT FOR 6DMHS

In this section, the holographic-based sensing method for RHS is first proposed. Then, by utilizing the sensing information, we propose the alternative optimization algorithm for adjusting the orientation of 6DMHS.

# A. Holographic-based Sensing Method

For aligning the RHS with the IDET receiver accurately, a low-error angle information of the IDET receiver is required. However, the tradition sensing method may not acquire the sensing information accurately. The reasons are summarized as follows: 1) The loss of the imaginary part of the holographic beamformer; 2) In each feed, the receive sensing signal is the series superposition of the uplink sensing signal from each

![](https://cdn-mineru.openxlab.org.cn/result/2025-11-24/7946ec08-51d9-4a1e-813e-6601cf747e28/66f88fcb1346531e5945b2b7c3368f3bd120e2d5a57eff0bac2112c193aba791.jpg)  
(a)  
(b)  
Fig. 4. Illustration of the sensing RHS: (a) The architecture of the sensing RHS; (b) The circuit of the sensing element.

radiation element. Next, we propose the holographic-based sensing method for high accurate sensing.

The holography-based sensing model is illustrated in Fig. 4(a). The sensing elements are arranged around the RHS with the spacing of  $d_{\mathrm{S}}$ . The holographic-based sensing circuit of each sensing element is illustrated as Fig. 4(b). When the uplink sensing signal is fed into the sensing circuit, it is superimposed with the reference signal. Then, this composite signal passes through the diode-based filter circuit, which retains only the direct current (DC) component. Next, the power meter records the DC power and uploads it to the BS. The collection of all recorded power values from all the power meters constitutes the holographic image, which encapsulates the CSI of all the IDET receivers. The generation of the reference signals is summarized as follows: 1) Firstly, the reference source generates the signal with constant amplitude and phase; 2) Secondly, this signals pass through the random phase shifters, resulting in reference signals with randomized phases. In the following, we describe how the holographic-based sensing is achieved by using the holographic image.

Remark 2: The sensing elements consist of the low-cost diode, reference source, low-precision random phase shifter and power meter, while they do not require the high-cost RF-chains and ADC modules [27]. This indicates that the cost of

the sensing elements is low and it is practical to equip them in the RHS for achieving higher sensing performance.

The wireless channel from the  $k$ -th IDET receiver to the sensing elements of the  $b$ -th RHS in the local coordinate system is given as

$$
\begin{array}{l} \mathbf {H} _ {k, b} ^ {\mathrm {S}, \mathrm {L}} = \\ \left[ \begin{array}{c c c c c} h _ {k, b, 0, 0} ^ {\mathrm {S}, \mathrm {L}} & h _ {k, b, 0, 1} ^ {\mathrm {S}, \mathrm {L}} & \dots & h _ {k, b, 0, N _ {y} - 2} ^ {\mathrm {S}, \mathrm {L}} & h _ {k, b, 0, N _ {y} - 1} ^ {\mathrm {S}, \mathrm {L}} \\ h _ {k, b, 1, 0} ^ {\mathrm {S}, \mathrm {L}} & 0 & \dots & 0 & h _ {k, b, 1, N _ {y} - 1} ^ {\mathrm {S}, \mathrm {L}} \\ \vdots & \vdots & \ddots & \vdots & \vdots \\ h _ {k, b, N _ {x} - 2, 0} ^ {\mathrm {S}, \mathrm {L}} & 0 & \dots & 0 & h _ {k, b, N _ {x} - 2, N _ {y} - 1} ^ {\mathrm {S}, \mathrm {L}} \\ h _ {k, b, N _ {x} - 1, 0} ^ {\mathrm {S}, \mathrm {L}} & h _ {k, b, N _ {x} - 1, 1} ^ {\mathrm {S}, \mathrm {L}} & \dots & h _ {k, b, N _ {x} - 1, N _ {y} - 2} ^ {\mathrm {S}, \mathrm {L}} & h _ {k, b, N _ {x} - 1, N _ {y} - 1} ^ {\mathrm {S}, \mathrm {L}} \end{array} \right], \tag {27} \\ \end{array}
$$

where  $h_{k,b,n_x,n_y}^{\mathrm{S,L}}$  is expressed as

$$
h _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S}, \mathrm {L}} = \left\{ \begin{array}{c} \sum_ {\iota} \Lambda_ {k, \iota , b} \eta_ {k, \iota , b} e ^ {j \frac {2 \pi}{\lambda} \mathbf {f} ^ {T} \left(\theta_ {k, \iota , b} ^ {\mathrm {L}}, \phi_ {k, \iota , b} ^ {\mathrm {L}}\right) \mathbf {r} _ {b, n _ {x}, n _ {y}} ^ {\mathrm {S}, \mathrm {L}}}, \text {i f} n _ {x} \cdot n _ {y} = 0 \text {o r} \\ n _ {x} = N _ {x} - 1 \text {o r} n _ {y} = N _ {y} - 1, \\ 0, \text {e l s e}, \end{array} \right. \tag {28}
$$

where  $(\theta_{k,\iota,b}^{\mathrm{L}}, \phi_{k,\iota,b}^{\mathrm{L}})$  denote the angles of the  $\iota$ -th scatterer associated with the  $k$ -th IDET receiver in the local coordinate system of the  $b$ -th RHS. The vector  $\mathbf{r}_{b,n_x,n_y}^{\mathrm{S,L}} \in \mathbb{R}^{3 \times 1}$  represents the coordinates of the  $(n_x, n_y)$ -th sensing element at the  $b$ -th RHS in the local coordinate system. The uplink receive signal from the  $k$ -th IDET receiver at the  $(n_x, n_y)$ -th sensing element of the  $b$ -th RHS is given as

$$
y _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S}} = \sqrt {P ^ {\mathrm {S}}} h _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S , L}} s _ {k} ^ {\mathrm {S}} + z _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S}}, \tag {29}
$$

where  $P^{\mathrm{S}}$  is the transmit power of the uplink sensing signal,  $s_k^{\mathrm{S}} = 1$  is the sensing signal and  $z_{k,b,n_x,n_y}^{\mathrm{S}}\sim CN(0,\sigma_s^2)$  is the antenna noise.

After being superimposed with the reference signal in the holographic-based sensing circuit, the output value of the  $(n_x,n_y)$ -th power meter at the  $b$ -th RHS is expressed as

$$
\begin{array}{l} P _ {k, b, n _ {x}, n _ {y}} ^ {S} \\ = \left| y _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S}} + s _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {r e f}} \right| ^ {2} \\ = \left| y _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S}} \right| ^ {2} + \left| s _ {k, b, n _ {x}, n _ {y}} ^ {\text {r e f}} \right| ^ {2} + y _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S} *} s _ {k, b, n _ {x}, n _ {y}} ^ {\text {r e f}} + y _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S}} s _ {k, b, n _ {x}, n _ {y}} ^ {\text {r e f} *} \\ \stackrel {\text {(a)}} {\approx} \left| s _ {k, b, n _ {x}, n _ {y}} ^ {\text {r e f}} \right| ^ {2} + 2 \mathcal {R} \left\{y _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S} *} s _ {k, b, n _ {x}, n _ {y}} ^ {\text {r e f}} \right\} \tag {30} \\ \end{array}
$$

where  $s_{k,b,n_x,n_y}^{\mathrm{ref}}$  represents the reference signal of the  $(n_x, n_y)$ -th sensing element at the  $b$ -th RHS. We have  $s_{k,b,n_x,n_y}^{\mathrm{ref}} = \sqrt{A} e^{j2\pi \chi_{k,b,n_x,n_y}}$ , where  $A$  is a constant and  $\chi_{k,b,n_x,n_y} \sim \mathcal{U}(-\sigma_1, \sigma_1)$ . In addition, (a) holds due to  $\left| y_{k,b,n_x,n_y}^{\mathrm{S}} \right| \ll \left| s_{k,b,n_x,n_y}^{\mathrm{ref}} \right|$ . Note that the reference signal is known at the BS. By subtracting the reference signal from Eq. (30), the holographic

Note that the power of the reference wave can be artificially adjusted, whereas the receive signal power is inherently low due to the limited uplink transmit power and the severe exponential path loss [27]. As a result, the reference wave power is significantly higher than that of the uplink signal.

image is given as

$$
\tilde {P} _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S}} = 2 \mathcal {R} \left\{y _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S} *} s _ {k, b, n _ {x}, n _ {y}} ^ {\text {r e f}} \right\}. \tag {31}
$$

Then, by exciting the holographic image with the reference signal, the angle information is extracted. The excited holographic image is expressed as

$$
\begin{array}{l} \mathcal {H} _ {k, b, n _ {x}, n _ {y}} = \tilde {P} _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S}} s _ {k, b, n _ {x}, n _ {y}} ^ {\text {r e f}} \\ = h _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S}, \mathrm {L}} + h _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S}, \mathrm {L} *} \left(s _ {k, b, n _ {x}, n _ {y}} ^ {\text {r e f}}\right) ^ {2} + z _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S}} + z _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S} *} \left(s _ {k, b, n _ {x}, n _ {y}} ^ {\text {r e f}}\right) ^ {2}. \tag {32} \\ \end{array}
$$

Theorem 1: When the number of the sensing elements is sufficient large, i.e.,  $N_x, N_y \to \infty$ , we have

$$
\begin{array}{l} \tilde {\mathcal {H}} _ {k, b} = \frac {1}{N} \sum_ {n _ {x}, n _ {y}} e ^ {- j \frac {2 \pi}{\lambda} \mathbf {f} ^ {T} \left(\theta_ {k, b} ^ {\mathrm {L , e}}, \phi_ {k, b} ^ {\mathrm {L , e}}\right) \mathbf {r} _ {b, n _ {x}, n _ {y}} ^ {\mathrm {S , L}}} \mathcal {H} _ {k, b, n _ {x}, n _ {y}} \\ = \frac {1}{N} \sum_ {n _ {x}, n _ {y}} e ^ {- j \frac {2 \pi}{\lambda} \mathbf {f} ^ {T} \left(\theta_ {k, b} ^ {\mathrm {L , c}}, \varphi_ {k, b} ^ {\mathrm {L , c}}\right) \mathbf {r} _ {b, n _ {x}, n _ {y}} ^ {\mathrm {S , L}}} \left[ h _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S , L}} + h _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S , L *}} \left(s _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {R E F}}\right) ^ {2} \right. \\ \left. + z _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S}} + z _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S} ^ {*}} \left(s _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {R E F}}\right) ^ {2} \right] \\ \end{array}
$$

$$
\stackrel {N \rightarrow \infty} {=} \left\{\begin{array}{l l}\Lambda_ {k, t, b} \eta_ {k, t, b},&\text {i f} f _ {1} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) = f _ {1} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right) \text {a n d}\\&f _ {2} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) = f _ {2} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right),\\0,&\text {e l s e ,}\end{array}\right. \tag {33}
$$

where  $N = 2N_{x} + 2N_{y} - 4$ .  $(\theta_{k,b}^{\mathrm{L,e}},\phi_{k,b}^{\mathrm{L,e}})$  denote the estimated angles of the  $k$  th IDET receiver with respect to the  $b$  -th RHS.

Proof 1: Please see Appendix A.

According to Theorem 1, the direction corresponding to the maximum beamforming gain for the  $k$ -th IDET receiver is given as

$$
\left(\theta_ {k, b ^ {*}} ^ {\mathrm {L}, \mathrm {e} ^ {*}}, \phi_ {k, b ^ {*}} ^ {\mathrm {L}, \mathrm {e} ^ {*}}, b ^ {*}\right) = \arg \max  _ {b, \theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}} \tilde {\mathcal {H}} _ {k, b}. \tag {34}
$$

By adopting the FFT-based detection method, Eq. (34) can be efficiently solved with low computational complexity.

Theorem 2: Densely arranged antennas do not improve the angular resolution of the RHS. Therefore, a half-wavelength spacing represents the optimal element configuration when using the FFT-based detection method.

Proof 2: Please see Appendix B.

According to Theorem 2, it is not necessary to install the sensing module densely. Instead, we can place the sensing elements with a half-wavelength spacing, which helps reduce costs and simplifies the wiring complexity.

Denote the matrix of the holographic image as

$$
\mathcal {H} _ {k, b} =
$$

$$
\left[ \begin{array}{c c c c c} \mathcal {H} _ {k, b, 0, 0} & \mathcal {H} _ {k, b, 0, 1} & \dots & \mathcal {H} _ {k, b, 0, N _ {y} - 2} & \mathcal {H} _ {k, b, 0, N _ {y} - 1} \\ \mathcal {H} _ {k, b, 1, 0} & 0 & \dots & 0 & \mathcal {H} _ {k, b, 1, N _ {y} - 1} \\ \vdots & \vdots & \ddots & \vdots & \vdots \\ \mathcal {H} _ {k, b, N _ {x} - 2, 0} & 0 & \dots & 0 & \mathcal {H} _ {k, b, N _ {x} - 2, N _ {y} - 1} \\ \mathcal {H} _ {k, b, N _ {x} - 1, 0} & \mathcal {H} _ {k, b, N _ {x} - 1, 1} & \dots & \mathcal {H} _ {k, b, N _ {x} - 1, N _ {y} - 2} & \mathcal {H} _ {k, b, N _ {x} - 1, N _ {y} - 1} \end{array} \right]. \tag {35}
$$

Then, by adopting the 2D-FFT, we have

$$
\left(n _ {x} ^ {*}, n _ {y} ^ {*}, b ^ {*}\right) = \arg \max  _ {n _ {x}, n _ {y}, b} \left| \mathbf {F} _ {N _ {x} \times N _ {x}} ^ {T} \mathcal {H} _ {k, b} \mathbf {F} _ {N _ {y} \times N _ {y}} \right|. \tag {36}
$$

According to (36), the estimated direction vector in the local coordinate system is given as

$$
\begin{array}{l} \mathbf {f} (\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}) = \left[ \frac {(2 n _ {y} ^ {*} - N _ {y} + 1) d _ {\mathrm {S}}}{2 N \lambda}, \frac {(2 n _ {x} ^ {*} - N _ {x} + 1) d _ {\mathrm {S}}}{2 N \lambda}, \right. \\ \left. \sqrt {1 - \left(\frac {\left(2 n _ {y} ^ {*} - N _ {y} + 1\right) d _ {\mathrm {S}}}{2 N \lambda}\right) ^ {2} - \left(\frac {\left(2 n _ {x} ^ {*} - N _ {x} + 1\right) d _ {\mathrm {S}}}{2 N \lambda}\right) ^ {2}} \right] ^ {T}. \tag {37} \\ \end{array}
$$

Moreover, the estimated direction vector in the global coordinate system is given as

$$
\mathbf {f} \left(\theta_ {k} ^ {\mathrm {e}}, \phi_ {k} ^ {\mathrm {e}}\right) = \bar {\mathbf {R}} _ {b ^ {*}} \mathbf {f} \left(\theta_ {k, b ^ {*}} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b ^ {*}} ^ {\mathrm {L}, \mathrm {e}}\right), \tag {38}
$$

where  $\bar{\mathbf{R}}_{b^*}$  is the initial rotation matrix for the  $b^{*}$ -th RHS. Finally, we utilize the sensing information to adjust the orientation of the 6DMHS for improving the IDET performance.

# B. Orientation Adjustment for 6DMHS

In our proposed protocol, the  $b$ -th RHS is adopted to serve the  $b$ -th IDET receiver. The parameters  $\bar{\mathbf{u}}_b$  and  $\mathbf{u}_b$  are given as

$$
\bar {\mathbf {u}} _ {b} = \left[ \cos \theta^ {*} \cos \phi^ {*}, \sin \theta^ {*} \cos \phi^ {*}, \sin \phi^ {*} \right] ^ {T}, \tag {39}
$$

$$
\mathbf {u} _ {b} = \left[ \cos \theta_ {k = b} ^ {\mathrm {e}} \cos \phi_ {k = b} ^ {\mathrm {e}}, \sin \theta_ {k = b} ^ {\mathrm {e}} \cos \phi_ {k = b} ^ {\mathrm {e}}, \sin \phi_ {k = b} ^ {\mathrm {e}} \right] ^ {T}. \tag {40}
$$

According to Eqs. (2)-(5), (39), and (40), the unitary rotation matrix of the  $b$ -th RHS, denoted by  $\tilde{\mathbf{R}}_b$ , is obtained. In the following,  $\tilde{\mathbf{R}}_b$  is employed to adjust the rotation of the  $b$ -th RHS. The translation position of the  $b$ -th RHS  $\mathbf{q}_b$  is given as

$$
\mathbf {q} _ {b} = \mathbf {s} _ {b} ^ {T} \mathbf {q}, \tag {41}
$$

where  $\mathbf{q}$  is the translation positions and is denoted as  $\mathbf{q} = \left[\mathbf{q}^{(0)},\dots ,\mathbf{q}^{(M_1 - 1)}\right]\in \mathbb{R}^{M_1\times 3}$ .  $\mathbf{q}^{(m_1)}\in \mathbb{R}^{3\times 1}$  is the coordinate of the  $m_{1}$ -th translation position. The selection vector  $\mathbf{s}_b\in \mathbb{R}^{M_1\times 1}$  is defined as

$$
s _ {b} \left[ m _ {1} \right] = \left\{ \begin{array}{l l} 1, & \text {i f} m _ {1} = m _ {1, b}, \\ 0, & \text {e l s e .} \end{array} \right. \tag {42}
$$

Here, the  $m_{1,b}$ -th position is selected for the  $b$ -th RHS. By substituting the obtained  $\tilde{\mathbf{R}}_b$  and  $\mathbf{q}_b$  into Eq. (13), the steering vector  $\mathbf{a}_b(\theta_{k = b}^{\mathrm{e}},\phi_{k = b}^{\mathrm{e}})$  is obtained. Then, the holographic beamformer of the  $b$ -th RHS is expressed as

$$
\bar {\Psi} _ {b} = \sum_ {q} \bar {\omega} _ {b, q} \frac {\mathcal {R} \left\{\sqrt {\frac {M}{\eta}} \operatorname {d i a g} \left(\boldsymbol {\Theta} _ {b , q}\right) \mathbf {a} _ {b} \left(\theta_ {k = b} ^ {\mathrm {e}}, \phi_ {k = b} ^ {\mathrm {e}}\right) \right\} + 1}{2}. \tag {43}
$$

The power gain of the  $k$ -th RHS at the estimated direction  $(\theta_k^{\mathrm{e}},\phi_k^{\mathrm{e}})$  is given as

$$
\dot {g} _ {k} \left(\theta_ {k} ^ {\mathrm {e}}, \phi_ {k} ^ {\mathrm {e}}\right) = \sum_ {k ^ {\prime}} \left| \sum_ {b} \mathbf {a} _ {b} \left(\theta_ {k} ^ {\mathrm {e}}, \phi_ {k} ^ {\mathrm {e}}\right) \operatorname {d i a g} \left(\bar {\boldsymbol {\Psi}} _ {b}\right) \boldsymbol {\Theta} _ {b} \mathbf {X} _ {k ^ {\prime}, b} \right| ^ {2}. \tag {44}
$$

The orientation optimization problem is formulated as

$$
\left(\mathrm {P} 1\right) \max  _ {\mathbf {s} _ {b}, \bar {\omega} _ {b, q}, \mathbf {X} _ {k, b}} \min  _ {k} \dot {g} _ {k} \left(\theta_ {k} ^ {\mathrm {e}}, \phi_ {k} ^ {\mathrm {e}}\right), \tag {45}
$$

$$
s. t. \quad \mathbf {s} _ {b} [ i ] \in \{0, 1 \}, \forall i, \tag {46a}
$$

$$
\mathbf {1} ^ {T} \mathbf {s} _ {b} = 1, \forall b, \tag {46b}
$$

$$
\sum_ {q} \bar {\omega} _ {b, q} = 1, \forall b, \tag {46c}
$$

$$
\sum_ {b} \| \mathbf {X} _ {b} \| _ {2} ^ {2} \leq P _ {\mathrm {t x}}, \forall b. \tag {46d}
$$

$$
(7) - (9).
$$

The constraint (46a) ensures that  $\mathbf{s}_b$  is a binary vector, while the constraint (46b) enforces that exactly one element of the vector  $\mathbf{s}_b$  is equal to 1, with all other elements being 0. The constraint (46c) ensures that the holographic beamformer is not greater than 1 and the constraint (46d) ensures that the transmit power of the 6DMHS-assisted transmitter is not greater than  $P_{\mathrm{tx}}$ . The constraints (7)-(9) refer to the rotation constraint for avoiding signal reflection, rotation constraint for avoiding signal blockage and minimum-distance constraint for avoiding collision, respectively. The optimization problem (P1) is nonconvex due to the nonconvex objective function and the binary constraint (46a). To address this problem, we adopt the alternating optimization algorithm to optimize  $\mathbf{s}_b$  and adopt the FP-based algorithm to optimize  $\bar{\omega}_{b,q}$  and  $\mathbf{X}_{k,b}$ .

Firstly, by fixing  $\bar{\omega}_{b,q}$  and  $\mathbf{X}_{k,b}$ , the optimization problem (P1) can be reformulated as

$$
\left(\mathrm {P} 2\right) \max  _ {\mathbf {s} _ {b}} \min  _ {k} \dot {g} _ {k} \left(\theta_ {k} ^ {\mathrm {e}}, \phi_ {k} ^ {\mathrm {e}}\right), \tag {46}
$$

$$
\begin{array}{l} \text {s . t .} \quad (4 6 a), (4 6 b), (7) - (9). \end{array}
$$

The optimization problem (P2) can be solved by using the alternating optimization algorithm [22], which is summarized as Algorithm 1. The complexity for solving Algorithm 1 is  $\mathcal{O}(M_1BKQ(M + K))$ .

Algorithm 1 Alternating optimization algorithm for solving (P2).

Input: Initial values of  $\mathbf{s}_b$ ,  $\bar{\omega}_{b,a}$ ,  $\mathbf{X}_{k,b}$ ;  
Output: The selection vector  $\mathbf{s}_b$ ;  
1: Update  $t \gets 1$ ,  $\overline{g}^{\mathrm{new}} \gets +\infty$ ,  $\overline{g}^{\mathrm{old}} \gets -\infty$ ,  $\tilde{g}^{\mathrm{old}} \gets -\infty$ ,  $\hat{\mathbf{s}}_b \gets \mathbf{s}_b$ ;  
2: while  $\left| \overline{g}^{\mathrm{new}} - \overline{g}^{\mathrm{old}} \right| \geq \epsilon \, \mathrm{d}\omega$  
3: Update  $\overline{g}^{\mathrm{old}}\gets \overline{g}^{\mathrm{new}}$  
4: for  $m_{1} = 0$  to  $M_1 - 1$  do  
5: for  $b = 0$  to  $B - 1$  do  
6: Update  $\mathbf{s}_{b'} \gets \hat{\mathbf{s}}_{b'}$ ,  $\forall b' = 0, \dots, B - 1$ ;  
7: Update  $\mathbf{s}_b\gets \mathbf{0}$ $s_b[m_1]\gets 1$ $\tilde{g}_{k}^{\mathrm{new}}\left(\theta_{k}^{\mathrm{e}},\phi_{k}^{\mathrm{e}}\right)$  
8: if  $\tilde{g}^{\mathrm{new}} > \tilde{g}^{\mathrm{old}}$  and the constraints (7)-(9) are satisfied then  
9: Update  $\tilde{g}^{\mathrm{old}}\gets \tilde{g}^{\mathrm{new}},\hat{\mathbf{s}}_b\gets \mathbf{s}_b;$  
10: end if  
11: end for  
12: Update  $\overline{g}^{\mathrm{new}}\gets \widetilde{g}^{\mathrm{new}}$  
13: end for  
14: end while  
15: return  $\{\mathbf{s}_b,\min_k\dot{g}_k(\theta_k^e,\phi_k^e)\}$

Secondly, by fixing  $\mathbf{s}_b$ , (P1) can be reformulated as

$$
\begin{array}{l} \left(\mathrm {P} 3\right) \max  _ {\bar {\omega} _ {b, q}, \mathbf {X} _ {k, b}} \min  _ {k} \dot {g} _ {k} \left(\theta_ {k} ^ {\mathrm {e}}, \phi_ {k} ^ {\mathrm {e}}\right), \tag {47} \\ \begin{array}{l l} \text {s . t .} & (4 6 c), (4 6 d). \end{array} \\ \end{array}
$$

By adopting the FP method, the objective function can be converted to be convex and (P3) can be solved effectively. Denote the vector  $\Xi_k \in \mathbb{C}^{K \times 1}$  as

$$
\begin{array}{l} \boldsymbol {\Xi} _ {k} = \left[ \sum_ {b} \mathbf {a} _ {b} \left(\theta_ {k} ^ {\mathrm {e}}, \phi_ {k} ^ {\mathrm {e}}\right) \operatorname {d i a g} \left(\bar {\boldsymbol {\Psi}} _ {b}\right) \boldsymbol {\Theta} _ {b} \mathbf {X} _ {0, b}, \dots , \right. \\ \left. \sum_ {b} \mathbf {a} _ {b} \left(\theta_ {k} ^ {\mathrm {e}}, \phi_ {k} ^ {\mathrm {e}}\right) \operatorname {d i a g} \left(\bar {\boldsymbol {\Psi}} _ {b}\right) \boldsymbol {\Theta} _ {b} \mathbf {X} _ {K - 1, b} \right] ^ {T}. \tag {48} \\ \end{array}
$$

where the vector  $\mathbf{a}_b(\theta, \phi)$  is given as

$$
\begin{array}{l} \mathbf {a} _ {b} (\theta , \phi) = \sqrt {\frac {1}{M}}. \\ \left[ e ^ {j \frac {2 \pi}{\lambda} \mathbf {f} ^ {T} (\theta , \phi) \left(\mathbf {s} _ {b} ^ {T} \mathbf {q} + \mathbf {R} _ {b} \bar {\mathbf {r}} _ {0, 0}\right)}, \dots , e ^ {j \frac {2 \pi}{\lambda} \mathbf {f} ^ {T} (\theta , \phi) \left(\mathbf {s} _ {b} ^ {T} \mathbf {q} + \mathbf {R} _ {b} \bar {\mathbf {r}} _ {0, M y - 1}\right)}, \dots , \right. \\ \left. \left. e ^ {j \frac {2 \pi}{\lambda} \mathbf {f} ^ {T} (\theta , \phi) \left(\mathbf {s} _ {b} ^ {T} \mathbf {q} + \mathbf {R} _ {b} \bar {\mathbf {r}} _ {M _ {x} - 1, 0}\right)}, \dots , e ^ {j \frac {2 \pi}{\lambda} \mathbf {f} ^ {T} (\theta , \phi) \left(\mathbf {s} _ {b} ^ {T} \mathbf {q} + \mathbf {R} _ {b} \bar {\mathbf {r}} _ {M _ {x} - 1, M _ {y} - 1}\right)} \right] ^ {T}, \right. \tag {49} \\ \end{array}
$$

Then, by introducing the auxiliary variable  $\zeta_k \in \mathbb{C}^{K \times 1}$ , the objective function of (P3) can be reformulated as

$$
\ddot {g} _ {k} \left(\theta_ {k} ^ {\mathrm {e}}, \phi_ {k} ^ {\mathrm {e}}\right) = 2 \Re \left\{\zeta_ {k} ^ {H} \Xi_ {k} \right\} - \zeta_ {k} ^ {H} \zeta_ {k}. \tag {50}
$$

The problem (P3) can be reformulated as

$$
\begin{array}{l} \left. \left(\mathrm {P} 4\right) \max  _ {\bar {\omega} _ {b, q}, \mathbf {X} _ {k, b}, \boldsymbol {\zeta} _ {k}} \min  _ {k} \ddot {g} _ {k} \left(\theta_ {k} ^ {\mathrm {e}}, \phi_ {k} ^ {\mathrm {e}}\right), \right. \tag {51} \\ \begin{array}{l l} \text {s . t .} & (4 6 c), (4 6 d). \end{array} \\ \end{array}
$$

Next, we solve (P4) by alternatively optimizing  $\bar{\omega}_{b,q}$ ,  $\mathbf{X}_{k,b}$  and  $\zeta_{k}$ .

By fixing  $\bar{\omega}_{b,q}$  and  $\zeta_{k}$ , (P4) can be reformulated as

$$
\begin{array}{l} \left(\mathrm {P} 4 - 1\right) \max  _ {\mathbf {X} _ {k, b}} \min  _ {k} \ddot {g} _ {k} \left(\theta_ {k} ^ {\mathrm {e}}, \phi_ {k} ^ {\mathrm {e}}\right), \tag {52} \\ \begin{array}{l} \text {s . t .} \\ \text {(4 6 d) .} \end{array} \\ \end{array}
$$

By fixing  $\mathbf{X}_{k,b}$  and  $\zeta_{k}$ , (P4) can be reformulated as

$$
\begin{array}{l} \left(\mathrm {P} 4 - 2\right) \max  _ {\bar {\omega} _ {b, q}} \min  _ {k} \ddot {g} _ {k} \left(\theta_ {k} ^ {\mathrm {e}}, \phi_ {k} ^ {\mathrm {e}}\right), \tag {53} \\ \begin{array}{l l} \text {s . t .} & (4 6 c). \end{array} \\ \end{array}
$$

The problems (P4-1) and (P4-2) are convex and can be solved by adopting the interior point method [28].

By fixing  $\bar{\omega}_{b,q}$  and  $\mathbf{X}_{k,b}$ , (P4) can be reformulated as

$$
\left(\mathrm {P} 4 - 3\right) \max  _ {\zeta_ {k}} \min  _ {k} \ddot {g} _ {k} \left(\theta_ {k} ^ {\mathrm {e}}, \phi_ {k} ^ {\mathrm {e}}\right). \tag {54}
$$

By adopting the method of derivation, the optimal  $\zeta_{k}$  is given as

$$
\zeta_ {k} = \Xi_ {k}. \tag {55}
$$

By alternatively optimizing  $\bar{\omega}_{b,q}$ ,  $\mathbf{X}_{k,b}$  and  $\Xi_{k}$ , the optimal solution to (P4) can be obtained, which also yields the optimal value of (P3). Then, by alternatively solving (P2) and (P3), the optimal solution of (P1) is obtained. The algorithm for solving (P1) is summarized in Algorithm 2. The complexity of Algorithm 2 is  $O\big(T_{\mathrm{ite},1}\big(M_1BKQ(M + K) + T_{\mathrm{ite},2}B^{3.5}K^{3.5}Q^{3.5})\big)$ , where  $T_{\mathrm{ite},1}$  and  $T_{\mathrm{ite},2}$  represent the number of the inner iterations and the number of the outer iterations, respectively.

# Algorithm 2 Orientation optimization algorithm.

Input: Initial values of  $\mathbf{s}_b$ ,  $\bar{\omega}_{b,q}$ ,  $\mathbf{X}_{b,k}$ ,  $\zeta_k$ ;

Output: The selection vector  $\mathbf{s}_b$  and the weight factor  $\bar{\omega}_{b,q}$

$$
\begin{array}{l} 1: \text {U p d a t e} \overline {{g}} ^ {\text {n e w}} \leftarrow + \infty , \overline {{g}} ^ {\text {o l d}} \leftarrow - \infty , \tilde {g} ^ {\text {n e w}} \leftarrow + \infty , \tilde {g} ^ {\text {o l d}} \leftarrow - \infty ; \\ 2: \text {w h i l e} \left| \overline {{g}} ^ {\text {n e w}} - \overline {{g}} ^ {\text {o l d}} \right| \geq \epsilon \mathbf {d o} \\ 3: \quad \text {U p d a t e} \bar {g} ^ {\text {o l d}} \leftarrow \bar {g} ^ {\text {n e w}}; \\ 4: \quad \text {U p d a t e} \mathbf {s} _ {b} \text {b y a d o p t i n g A l g o r m} 1; \\ 5: \quad \text {w h i l e} \left| \tilde {g} ^ {\text {n e w}} - \tilde {g} ^ {\text {o l d}} \right| \geq \epsilon \mathbf {d o} \\ 6: \quad \text {U p d a t e} \tilde {g} ^ {\text {o l d}} \leftarrow \tilde {g} ^ {\text {n e w}}; \\ 7: \quad \text {U p d a t e} \mathbf {X} _ {k, b} \text {b y s o l v i n g (P 4 - 1)}; \\ 8: \quad \text {U p d a t e} \omega_ {b, q} \text {b y s o l v i n g (P 4 - 2)}; \\ 9: \quad \text {U p d a t e} \zeta_ {k} \text {b y E q . (5 5)}; \\ 1 0: \quad \text {U p d a t e} \tilde {g} ^ {\text {n e w}} \leftarrow \min  _ {k} \ddot {g} _ {k} \left(\theta_ {k} ^ {\mathrm {e}}, \phi_ {k} ^ {\mathrm {e}}\right); \\ 1 1: \quad \text {e n d w h i l e} \\ 1 2: \quad \text {U p d a t e} \bar {g} ^ {\text {n e w}} \leftarrow \min  _ {k} \dot {g} _ {k} \left(\theta_ {k} ^ {\mathrm {e}}, \phi_ {k} ^ {\mathrm {e}}\right); \\ 1 3: \text {e n d w h i l e} \\ 1 4: \text {r e t u r n} \left\{\mathbf {s} _ {b}, \bar {\omega} _ {b, q}, \mathbf {X} _ {b, k}, \zeta_ {k}, \min  _ {k} \dot {g} _ {k} \left(\theta_ {k} ^ {\mathrm {e}}, \phi_ {k} ^ {\mathrm {e}}\right) \right\}. \\ \end{array}
$$

# V. JOINT BEAMFORMING AND POWER SPLITTER DESIGN

In this section, we jointly optimize the digital beamforming and the power splitter based on the equivalent CSI. During the orientation adjustment stage, the optimal orientation of the 6DMHS and the holographic beamforming of the RHS are determined. Subsequently, in the downlink transmission stage, the equivalent CSI at each feed is estimated and exploited for IDET services. According to (10), the equivalent channel  $\bar{\mathbf{h}}_{k,b} \in \mathbb{C}^{1 \times Q}$  is expressed as

$$
\bar {\mathbf {h}} _ {k, b} = \sqrt {M} \sum_ {\iota} \Lambda_ {k, t, b} \eta_ {k, t, b} \mathbf {a} _ {b} ^ {T} \left(\theta_ {k, t, b}, \phi_ {k, t, b}\right) \operatorname {d i a g} \left(\Psi_ {b}\right) \boldsymbol {\Theta} _ {b}. \tag {56}
$$

In the sequel, the equivalent CSI is estimated within each short IDET frame, after which the downlink digital beamforming and power splitter are jointly designed to maximize the IDET performance.

The minimal EH power maximization problem is formulated as

$$
\left(\mathrm {P} 5\right) \max  _ {\mathbf {X} _ {k, b}, \rho_ {k}} \min  _ {k} \Gamma \left(P _ {\mathrm {E H}, k}\right), \tag {57}
$$

$$
\text {s . t .} \quad R _ {k} \geq R _ {0}, \tag {57a}
$$

$$
\sum_ {b, k} \left\| \mathbf {X} _ {k, b} \right\| _ {2} ^ {2} \leq P _ {\mathrm {t x}}. \tag {57b}
$$

$$
0 \leq \rho_ {k} \leq 1, \tag {57c}
$$

(P5) is a nonconvex problem, due to the objective function and the constraint (57b) are nonconvex.

Firstly, by introducing the auxiliary variable  $\vartheta_{k}\in \mathbb{C}^{1}$ ,  $\gamma_{k}$  can be reformulated as

$$
\begin{array}{l} \dot {\gamma} _ {k} = 2 \Re \left\{\sqrt {1 - \rho_ {k}} \vartheta_ {k} ^ {H} \sum_ {b} \overline {{\mathbf {h}}} _ {k, b} \mathbf {X} _ {k, b} \right\} \\ \left. - \vartheta_ {k} ^ {H} \vartheta_ {k} \left\{\left(1 - \rho_ {k}\right) \sum_ {k ^ {\prime} \neq k} \left| \sum_ {b} \overline {{\mathbf {h}}} _ {k, b} \mathbf {X} _ {k ^ {\prime}, b} \right| ^ {2} + \left(1 - \rho_ {k}\right) \sigma_ {0} ^ {2} + \sigma_ {\mathrm {c o v}} ^ {2} \right\}. \right. \tag {58} \\ \end{array}
$$

The downlink throughput  $R_{k}$  can be reformulated as

$$
\dot {R} _ {k} = \log_ {2} (1 + \dot {\gamma} _ {k}), [ \mathrm {b i t / s / H z} ]. \tag {59}
$$

Secondly, by introducing the variable  $\pmb{\varsigma}_k\in \mathbb{C}^{K\times 1}$ ,  $P_{\mathrm{EH},k}$  is reformulated as

$$
\dot {P} _ {\mathrm {E H}, k} = 2 \Re \left\{\boldsymbol {\varsigma} _ {k} ^ {H} \hat {\vec {\Xi}} _ {k} \right\} - \boldsymbol {\varsigma} _ {k} ^ {H} \boldsymbol {\varsigma} _ {k} + \rho_ {k} \sigma_ {0} ^ {2}, \tag {60}
$$

where  $\dot{\Xi}_k = \sqrt{\rho_k}\left[\sum_b\overline{\mathbf{h}}_{k,b}\mathbf{X}_{b,0},\dots ,\sum_b\overline{\mathbf{h}}_{k,b}\mathbf{X}_{b,K - 1}\right]^T\in \mathbb{C}^{K\times 1}$

Then, (P5) can be reformulated as

$$
\left(\mathrm {P} 6\right) \max  _ {\mathbf {X} _ {k, b}, \rho_ {k}, \vartheta_ {k}, \varsigma_ {k}} P _ {0}, \tag {61}
$$

$$
s. t. \quad \Gamma (\dot {P} _ {\mathrm {E H}, k}) \geq P _ {0} \tag {61a}
$$

$$
\dot {R} _ {k} \geq R _ {0}, \tag {61b}
$$

$$
(5 7 b) - (5 7 c).
$$

The constraint (61a) is nonconvex. By using the inverse function of the function  $\Gamma (\dot{P}_{\mathrm{EH},k})$ , the constraint (61a) can be reformulated as

$$
\dot {P} _ {\mathrm {E H}, k} \geq \Gamma^ {- 1} \left(P _ {0}\right), \tag {62}
$$

where  $\Gamma^{-1}(\cdot)$  is the inverse function of  $\Gamma (\cdot)$ . By letting  $\dot{P}_0 = \Gamma^{-1}(P_0)$ , (P6) can be reformulated as

$$
\left(\mathrm {P} 7\right) \max  _ {\mathbf {X} _ {k, b}, \rho_ {k}, \theta_ {k}, \varsigma_ {k}} \dot {P} _ {0}, \tag {63}
$$

$$
\text {s . t .} \quad \dot {P} _ {\mathrm {E H}, k} \geq \dot {P} _ {0} \tag {63a}
$$

$$
(6 1 b), (5 7 b) - (5 7 c).
$$

Next, by alternatively optimizing  $\mathbf{X}_{k,b},\rho_k,\vartheta_k$  and  $\varsigma_{k}$ , (P7) can be solved effectively.

By fixing  $\rho_{k},\vartheta_{k}$  and  $\varsigma_{k}$  (P7) is reformulated as

$$
\left(\mathrm {P} 7 - 1\right) \max  _ {\mathbf {X} _ {k, b}} \dot {P} _ {0}, \tag {64}
$$

$$
\begin{array}{l} \text {s . t .} \\ \hline (6 3 \mathrm {a}), (6 1 \mathrm {b}), (5 7 \mathrm {b}). \end{array}
$$

By fixing  $\mathbf{X}_{k,b},\vartheta_k$  and  $\varsigma_{k}$  (P7) can be reformulated as

$$
\left(\mathrm {P} 7 - 2\right) \max  _ {\rho_ {k}} \dot {P} _ {0}, \tag {65}
$$

$$
\begin{array}{l} \text {s . t .} \\ \hline (6 3 \mathrm {a}), (6 1 \mathrm {b}), (5 7 \mathrm {c}). \end{array}
$$

The problems (P7-1) and (P7-2) are convex, which can be solved by using the interior point method. By fixing  $\mathbf{X}_{k,b},\rho_{k}$  and  $\varsigma_{k}$ , (P7) can be reformulated as

$$
\left(\mathrm {P} 7 - 3\right) \max  _ {\vartheta_ {k}} \dot {P} _ {0}, \tag {66}
$$

$$
\begin{array}{l} \text {s . t .} \\ \hline \end{array} (6 1 \mathrm {b}).
$$

By adopting the method of derivation, the optimal  $\vartheta_{k}^{*}$  is given as

$$
\vartheta_ {k} ^ {*} = \frac {\sqrt {1 - \rho_ {k}} \sum_ {b} \overline {{\mathbf {h}}} _ {k , b} \mathbf {X} _ {k , b}}{(1 - \rho_ {k}) \sum_ {k ^ {\prime} \neq k} \left| \sum_ {b} \overline {{\mathbf {h}}} _ {k , b} \mathbf {X} _ {k ^ {\prime} , b} \right| ^ {2} + (1 - \rho_ {k}) \sigma_ {0} ^ {2} + \sigma_ {\mathrm {c o v}} ^ {2}}. \tag {67}
$$

By fixing  $\mathbf{X}_{k,b},\rho_k$  and  $\vartheta_{k}$  (P7) can be reformulated as

$$
\left(\mathrm {P} 7 - 4\right) \max  _ {\varsigma_ {k}} \dot {P} _ {0}, \tag {68}
$$

$$
\begin{array}{l l} \text {s . t .} & (6 3 \mathrm {a}). \end{array}
$$

Similarly, by adopting the method of derivation, the optimal  $S_{k}^{*}$  is given as

$$
\mathbf {s} _ {k} ^ {*} = \dot {\Xi} _ {k}. \tag {69}
$$

By alternatively solving (P7-1)-(P7-4), the optimal solutions of (P7) are obtained, which are also the optimal solutions of (P5) and (P6). The FP-based algorithm for solving (P5) is summarized in Algorithm 3. The complexity of the Algorithm

3 is  $O\big(T_{\mathrm{ite},3}\big(B^{3.5}Q^{3.5}K^{3.5}\big)\big)$ , where  $T_{\mathrm{ite},3}$  represents the number of the iterations.

Algorithm 3 FP-based optimization algorithm for solving (P2).

Input: Initial value of  $\mathbf{X}_{k,b}$ ,  $\rho_{k}$ ,  $\vartheta_{k}$  and  $\pmb{\varsigma}_{k}$ ;

Output: The digital beamforming vector  $\mathbf{X}_{k,b}$  and the power splitting factor  $\rho_{k}$ ;

$$
\begin{array}{l} 1: \quad \text {U p d a t e} \overline {{P}} _ {0} ^ {\text {n e w}} \leftarrow + \infty , \overline {{P}} _ {0} ^ {\text {o l d}} \leftarrow - \infty ; \\ 2: \text {w h i l e} \left| \overline {{P}} _ {0} ^ {\text {n e w}} - \overline {{P}} _ {0} ^ {\text {o l d}} \right| \geq \epsilon \text {d o} \\ 3: \quad \text {U p d a t e} \overline {{P}} _ {0} ^ {\text {o l d}} \leftarrow \overline {{P}} _ {0} ^ {\text {n e w}}; \\ 4: \quad \text {U p d a t e} \mathbf {X} _ {k, b} ^ {\prime} \text {b y s o l v i n g (P 7 - 1)}; \\ 5: \quad \text {U p d a t e} \rho_ {k} \text {b y s o l v i n g (P 7 - 2)}; \\ 6: \quad \text {U p d a t e} \vartheta_ {k} \text {b y E q . (6 7)}; \\ 7: \quad \text {U p d a t e} \underline {{\underline {{s}}}} _ {\text {h e w}} \text {b y E q .} (6 9); \\ 8: \quad \text {U p d a t e} P _ {0} ^ {\text {m e n}} \leftarrow \min  _ {k} \Gamma \left(P _ {\mathrm {E H}, k}\right); \\ 9: \text {e n d w h i l e} \\ 1 0: \text {r e t u r n} \left\{\mathbf {X} _ {k, b}, \rho_ {k}, \vartheta_ {k}, \boldsymbol {\varsigma} _ {k}, \min  _ {k} \Gamma \left(P _ {\mathrm {E H}, k}\right) \right\}. \\ \end{array}
$$

# VI. SIMULATION RESULTS

In this section, numerical results are presented to validate the performance of the proposed 6DMHS-assisted IDET system. The 3GPP channel model is employed, wherein the channel between the 6DMHS-assisted transmitter and the IDET receivers comprises one LoS path and three NLoS paths. The parameters are set as follows:  $f_{c} = 30 \mathrm{GHz}$ ,  $c = 3 \times 10^{8}$  m/s,  $M = 32 \times 32$ ,  $Q = 1$ ,  $M_{1} = 50$ ,  $P_{t} = 40$  dBm,  $K_{\mathrm{R}} = 10$ ,  $\sigma_{0}^{2} = -100$  dBm and  $\sigma_{\mathrm{cov}}^{2} = -50$  dBm.

Some benchmarks are listed as follows:

- FPA: In this scheme, the positions and rotations of all the RHS are fixed. The RHS is uniformly distributed on the surface of a sphere with a radius of  $1\mathrm{m}$ . Furthermore, the orientation of each RHS is aligned with the normal vector pointing from the BS to the center of the corresponding RHS.  
- 6DMHS with rotation only: In this scheme, the positions of RHS is fixed and each RHS is uniformly placed on the surface of a sphere with a radius of  $1\mathrm{m}$ . Additionally, the orientation of each RHS is aligned with the sensing direction of its corresponding IDET receiver.  
- 6DMHS with translation only: In this scheme, the rotations of RHS is fixed along the normal vector direction from the BS to the center of each RHS. The discrete center positions of the RHS are located on the surface of a sphere with a radius of  $1\mathrm{m}$  and are optimized to maximize the beamforming gain.  
- Least square (LS)-based sensing [29]: This scheme utilizes the least square-based sensing method to estimate the angles of the IDET receivers, while the positions and rotations of the 6DMHS-assisted transmitter are optimized using Algorithm 2 proposed in this paper.

Fig. 5 illustrates the variation of the EH power with respect to the downlink transmit power. As the transmit power increases, the EH power gradually improves. Meanwhile, the proposed scheme is compared with the 6DMHS with rotation only, 6DMHS with translation only and FPA benchmarks. The proposed scheme achieves enhanced spatial degrees of freedom by jointly optimizing the rotations and translations of

![](https://cdn-mineru.openxlab.org.cn/result/2025-11-24/7946ec08-51d9-4a1e-813e-6601cf747e28/134890d0b8d2424683bba81ed0ccad86397764fc62fe36f2acd3d4ac46a68b98.jpg)  
Fig. 5. EH power performance versus the transmit power under different transmitter configurations.

the 6DMHS, thereby attaining the best IDET performance. The 6DMHS with rotation only scheme outperforms the 6DMHS with translation only scheme, as the rotation operation enables alignment of the direction of maximum beamforming gain with the corresponding IDET receiver, resulting in improved IDET performance. Furthermore, both the 6DMHS with rotation only and 6DMHS with translation only schemes outperform the FPA scheme, since the added mobility in these schemes provide greater flexibility, leading to superior IDET performance.

Fig. 6 illustrates the trade-off between the downlink throughput threshold and the EH power. Firstly, we define the root-mean-square error of the sensing angles as  $\mathrm{RMSE} = \frac{1}{T_{\mathrm{s}}}\sum_{t_{\mathrm{s}} = 0}^{T_{\mathrm{s}} - 1}\left(\| \mathbf{f}(\theta_k^{\mathrm{e}},\phi_k^{\mathrm{e}}) - \mathbf{f}(\theta_k^{\mathrm{t}},\phi_k^{\mathrm{t}})\| _2^2\right)$ , where  $\theta_{k}^{\mathrm{t}}$  and  $\phi_k^{\mathrm{t}}$  are the true azimuth and elevation angles of the  $k$ -th IDET receiver or its scatterer with the maximal channel power gain. As the throughput threshold increases, the EH power gradually decreases. This is because a larger portion of the received RF power for the IDET receiver is allocated for information decoding, thereby reducing the input power of the energy harvester. Additionally, the impact of angular sensing accuracy on the IDET performance is analyzed. When there is no sensing error, i.e.,  $\mathrm{RMSE} = 0$ , the optimal IDET performance is achieved. As the sensing error increases from  $\mathrm{RMSE} = 0$  to  $\mathrm{RMSE} = 0.17$ , the performance steadily degrades. This result highlights the critical role of uplink sensing accuracy in downlink throughput, indicating that the higher sensing precision leads to an improved IDET performance.

Fig. 7 illustrates the impact of the throughput threshold on the IDET performance under different sensing elements spacing and different sensing scheme. It can be observed that when adopts the proposed sensing method and  $d_{\mathrm{S}} \leq \frac{\lambda}{2}$ , the IDET performance is approximately equivalent under different spacing of sensing elements. Moreover, the proposed holographic-based sensing method consistently outperforms the LS-based sensing method in terms of IDET performance, indicating that our proposed sensing method has a better

![](https://cdn-mineru.openxlab.org.cn/result/2025-11-24/7946ec08-51d9-4a1e-813e-6601cf747e28/030f70dcf1b181bb95a0ace820ff2082706fdb50c070c92e27a75784ddb08933.jpg)  
Fig. 6. EH power performance versus the throughput threshold under different sensing errors.

![](https://cdn-mineru.openxlab.org.cn/result/2025-11-24/7946ec08-51d9-4a1e-813e-6601cf747e28/08a91a076e066d3c14c596f3f9c8170fadad3f55481736ce3d336599c7661d3c.jpg)  
Fig. 7. EH power performance versus the adjustment factor under sensing methods.

sensing performance. However, when  $d_{\mathrm{S}} > \frac{\lambda}{2}$ , the IDET performance declines sharply. This degradation is attributed to the deterioration in sensing accuracy, as detailed below. Fig. 8 presents the sensing accuracy under different sensing element spacing. When  $d_{\mathrm{S}} \leq \frac{\lambda}{2}$ , the proposed holographic-based sensing method achieves at least a  $95\%$  probability that the RMSE remains below 0.1. In contrast, when  $d_{\mathrm{S}} > \frac{\lambda}{2}$ , the sensing accuracy of the proposed method deteriorates significantly. This observation supports Lemma 2, which states that placing sensing elements at half-wavelength intervals is sufficient to achieve near-optimal sensing performance.

Additionally, Fig. 8 shows that the LS-based sensing method performs worse than the proposed holographic-based sensing method. This is because the LS-based sensing method derives sensing information from the signals of the feeds, which aggregates the responses of all RHS elements. The signal aliasing resulting from this aggregation degrades angu

![](https://cdn-mineru.openxlab.org.cn/result/2025-11-24/7946ec08-51d9-4a1e-813e-6601cf747e28/de2a048fcd18725d5c33281639bec8813c25db6f06a7f7c4e294b89c8c5435e6.jpg)  
Fig. 8. CDF of the sensing error under different  $\kappa$  and different sensing methods.

![](https://cdn-mineru.openxlab.org.cn/result/2025-11-24/7946ec08-51d9-4a1e-813e-6601cf747e28/a812e51c1111828affe6209fac1bdb90b852c837fd13909f056fc82ed025b9c2.jpg)  
Fig. 9. EH power versus the transmit power under different angle alignment and feed positioning methods.

lar estimation accuracy. Moreover, the lacks of the imaginary part of the holographic beamformer may hurt the sensing performance of the LS-based sensing method. Consequently, the IDET performance adopting the LS-based sensing method reduces. In contrast, the proposed sensing method acquires sensing data in parallel from individual sensing elements, thereby avoiding inter-element aliasing and achieving superior sensing accuracy. As a result, the proposed holographic-based sensing method enables a higher IDET performance.

Fig. 9 compares the IDET performance under different feed position configurations and alignment strategies. It is observed that aligning the maximum beamforming gain direction of the RHS with the IDET receivers yields better IDET performance than simply aligning the normal vector of the RHS with the receiver direction. This confirms that directing the RHS's peak beamforming gain toward the receivers can significantly enhance the IDET performance, thereby validating the ef

![](https://cdn-mineru.openxlab.org.cn/result/2025-11-24/7946ec08-51d9-4a1e-813e-6601cf747e28/314434fdfbb79c9bcdb1c2ba8e85ca1e3249719582c38ba8bc88e327919bba12.jpg)  
Fig. 10. The EH power performance versus the number of the RHS element under different number of feeds.

![](https://cdn-mineru.openxlab.org.cn/result/2025-11-24/7946ec08-51d9-4a1e-813e-6601cf747e28/24139e46202b3eb1b2261dc9294aee93f895f40a49714b23f872b2b28b034c58.jpg)  
Fig. 11. EH power performance versus the rice factor under different protocols.

fectiveness of the proposed alignment strategy. Moreover, optimizing the feed positions further improves the IDET performance. This is because the optimized feed positions enable the RHS elements to achieve coherent signal combination, thereby increasing the radiated power in the direction of the main beam. As a result, the IDET performance is substantially enhanced. These results highlight the importance of feed position optimization during the RHS design and preparation process.

Fig. 10 illustrates the impact of the number of the RHS elements and feeds on the IDET performance. It is observed that as the number of elements increases, the EH power increases monotonically. This is because more elements contribute to a higher spatial gain, thereby enhancing the IDET performance. Moreover, increasing the number of feeds from 1 to 4 leads to a gradual improvement in IDET performance, due to the greater degrees of freedom in digital beamforming.

Assume that the CSI remains unchanged within a coherence time block of  $T_{c} = 140$  [30]. For our proposed scheme, the pilot overhead is given by  $T_{s,1} = \alpha KS\log_2(BQ)$ , and remaining time allocated for IDET transmission is  $T_{c} - T_{s,1}$ , where  $\alpha = 0.32$  is a constant parameter [31]. In contrast, for the perfect CSI scheme, the pilot overhead is  $T_{s,2} = \alpha KS\log_2(BM)$  [31] and the transmission time is  $T_{c} - T_{s,2}$ . Since  $Q \ll M$ , it follows that  $T_{s,1} \ll T_{s,2}$ . Fig. 11 illustrates the impact of the Rician factor on IDET performance. As the Rician factor increases from  $-10$  dB to  $-2$  dB, the IDET performance of our proposed scheme decreases. This is because when the Rician factor is low, the channel power gain is mainly concentrated on a scatterer having higher channel power gain, while the LoS link and the NLoS links generated by the other scatterers can be ignored. However, when we increase the Rician factor, the channel power gain of the LoS link and the NLoS links increase and cannot be ignored, leading to a worse IDET performance. Similarly, as we increase the Rician factor from  $-2$  dB to  $10$  dB, the IDET performance increases monotonously, since the channel power gain is gradually concentrated on the LoS link. Additionally, we compare our proposed scheme with the LoS-Only scheme. It is observed that when the Rician factor is greater than 2 dB, the LoS-Only scheme has the similar IDET performance with our proposed scheme. This is because both the two schemes design the holographic beamformer by adopting the LoS link. However, when the Rician factor is lower than 2 dB, our proposed scheme outperforms the LoS-Only scheme. This is because our proposed scheme designs the holographic beamformer by utilizing the NLoS link with the maximal channel power gain, which is higher than the channel power gain of the LoS link. Moreover, the IDET performance for the Perfect CSI scheme is almost unchanged when we increase the Rician factor from  $-10$  dB to  $10$  dB. This is indicates that designing the holographic beamformer with the perfect CSI can obtain the more robust IDET performance. Considering the overhead of the channel estimation, the IDET performance of the Perfect CSI scheme inferior to our proposed scheme when the Rician factor is greater than 2 dB or lower than  $-6$  dB. This is because the channel power gain is concentrated on the LoS link or the NLoS link with the maximal channel power gain and designing the holographic beamforming with the perfect CSI may not increase the IDET performance significantly, while the overwhelming overhead for the Perfect CSI scheme may decrease its IDET performance. However, when the Rician factor is in range  $-6$  dB to 2 dB, designing the holographic beamforming with the perfect CSI can significantly increase the IDET performance. This is because the multi-path interference is sever for the IDET receiver, and precisely designing the holographic beamforming with the perfect CSI may restrain the interference and improve the IDET performance effectively.

# VII. CONCLUSION

This paper investigated a sensing-enhanced 6DMHS-assisted IDET system. A holographic-based sensing method was proposed to acquire the angular information of IDET re

ceivers. Leveraging the obtained sensing information, the holographic beamforming of the RHS together with the rotations and positions of the 6DMHS were jointly optimized to align the maximum beamforming direction of the RHS with the IDET receiver. Subsequently, by fixing the holographic beamforming and the 6DMHS orientation, the downlink equivalent CSI was obtained, while the digital beamforming and power splitting factor were further optimized. As a result, the IDET performance was significantly improved. Simulation results show that: 1) The proposed parallel-form holographic sensing method achieves superior sensing performance over the conventional series-form method; 2) Aligning the maximum beamforming direction of the RHS with the IDET receiver yields notable IDET performance gains; 3) The proposed three-stage protocol outperforms the traditional perfect-CSI benchmark in the presence of an absolutely dominant path in the wireless channel, owing to its reduced pilot overhead.

# APPENDIX A PROOF OF THEOREM 1

Firstly, the limitation of  $\lim_{N\to \infty}\frac{1}{N}\sum_{n_x,n_y}e^{j\frac{2\pi}{A}\mathbf{f}^T (\theta_{k,b}^{L,e},\phi_{k,b}^{L,e})\mathbf{r}_{b,n_x,n_y}^{S,L}}h_{k,b,n_x,n_y}^{S,L}$  is expressed as Eq. (74), where  $\chi_1,\chi_2,\chi_3$  and  $\chi_4$  are given as

$$
\begin{array}{l} \mathcal {X} _ {1} = \frac {- N _ {x} + 1}{2} d _ {\mathrm {S}} \left(- f _ {1} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {1} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right)\right) \\ + \frac {- N _ {y} + 1}{2} d _ {\mathrm {S}} \left(- f _ {2} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {2} \left(\theta_ {k, \iota , b} ^ {\mathrm {L}}, \phi_ {k, \iota , b} ^ {\mathrm {L}}\right)\right), \tag {70} \\ \end{array}
$$

$$
\begin{array}{l} \mathcal {X} _ {2} = \frac {- N _ {x} + 1}{2} d _ {\mathrm {S}} \left(- f _ {1} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {1} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right)\right) \\ + \frac {N _ {\mathrm {y}} - 1}{2} d _ {\mathrm {S}} \left(- f _ {2} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {2} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right)\right), \tag {71} \\ \end{array}
$$

$$
\begin{array}{l} \mathcal {X} _ {3} = \frac {N _ {x} - 1}{2} d _ {\mathrm {S}} \left(- f _ {1} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {1} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right)\right) \\ + \frac {- N _ {y} + 1}{2} d _ {\mathrm {S}} \left(- f _ {2} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {2} \left(\theta_ {k, \iota , b} ^ {\mathrm {L}}, \phi_ {k, \iota , b} ^ {\mathrm {L}}\right)\right), \tag {72} \\ \end{array}
$$

$$
\begin{array}{l} \mathcal {X} _ {4} = \frac {N _ {x} - 1}{2} d _ {\mathrm {S}} \left(- f _ {1} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {1} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right)\right) \\ + \frac {N _ {\mathrm {y}} - 1}{2} d _ {\mathrm {S}} \left(- f _ {2} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {2} \left(\theta_ {k, \iota , b} ^ {\mathrm {L}}, \phi_ {k, \iota , b} ^ {\mathrm {L}}\right)\right). \tag {73} \\ \end{array}
$$

Then, the limitation of  $\lim_{N\to \infty}\frac{1}{N}\sum_{n_x,n_y}e^{j\frac{2\pi}{\lambda}\mathbf{r}^T (\theta_{k,b}^{\mathrm{L,e}},\theta_{k,b}^{\mathrm{L,e}})\mathbf{r}_{b,n_x,n_y}^{\mathrm{S,L}}}h_{k,b,n_x,n_y}^{\mathrm{S,L*}}(s_{k,b,n_x,n_y}^{\mathrm{ref}})^2$  is expressed as Eq. (79), where  $\mathcal{Y}_1,\mathcal{Y}_2,\mathcal{Y}_3$  and  $\mathcal{Y}_4$  are

$$
\begin{array}{l} \frac {1}{N} \sum_ {n _ {x}, n _ {y}} e ^ {- j \frac {2 \pi}{\lambda} \mathbf {f} ^ {T} (\theta_ {k, b} ^ {\mathrm {L , e}}, \phi_ {k, b} ^ {\mathrm {L , e}}) \mathbf {r} _ {b, n _ {x}, n _ {y}} ^ {\mathrm {S , L}}} h _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S , L}} \\ = \frac {1}{N} \sum_ {\iota} \Lambda_ {k, t, b} \eta_ {k, \iota , b} \sum_ {n _ {x}} e ^ {j \frac {2 \pi}{\lambda} \left[ n _ {x} d _ {S} \left(- f _ {1} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {1} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right)\right) + \mathcal {X} _ {1} \right]} + \frac {1}{N} \sum_ {\iota} \Lambda_ {k, t, b} \eta_ {k, \iota , b} \sum_ {n _ {x}} e ^ {j \frac {2 \pi}{\lambda} \left[ n _ {x} d _ {S} \left(- f _ {1} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {c}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {c}}\right) + f _ {1} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right)\right) + \mathcal {X} _ {2} \right]} \\ + \frac {1}{N} \sum_ {\iota} \Lambda_ {k, \iota , b} \eta_ {k, \iota , b} \sum_ {n _ {\mathrm {y}}} e ^ {j \frac {2 \pi}{\lambda} \left[ X _ {3} + n _ {\mathrm {y}} d _ {\mathrm {S}} \left(- f _ {2} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {2} \left(\theta_ {k, \iota , b} ^ {\mathrm {L}}, \phi_ {k, \iota , b} ^ {\mathrm {L}}\right)\right) \right]} + \frac {1}{N} \sum_ {\iota} \Lambda_ {k, \iota , b} \eta_ {k, \iota , b} \sum_ {n _ {\mathrm {y}}} e ^ {j \frac {2 \pi}{\lambda} \left[ X _ {4} + n _ {\mathrm {y}} d _ {\mathrm {S}} \left(- f _ {2} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {2} \left(\theta_ {k, \iota , b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right)\right) \right]} \\ - \frac {1}{N} \sum_ {\iota} \Lambda_ {k, \iota , b} \eta_ {k, \iota , b} e ^ {j \frac {2 \pi}{\lambda} X _ {1}} - \frac {1}{N} \sum_ {\iota} \Lambda_ {k, \iota , b} \eta_ {k, \iota , b} e ^ {j \frac {2 \pi}{\lambda} X _ {2}} - \frac {1}{N} \sum_ {\iota} \Lambda_ {k, \iota , b} \eta_ {k, \iota , b} e ^ {j \frac {2 \pi}{\lambda} X _ {3}} - \frac {1}{N} \sum_ {\iota} \Lambda_ {k, \iota , b} \eta_ {k, \iota , b} e ^ {j \frac {2 \pi}{\lambda} X _ {4}} \\ = \frac {1}{N} \sum_ {\iota} \Lambda_ {k, t, b} \eta_ {k, t, b} e ^ {j \frac {\pi}{\lambda} (\mathcal {X} _ {1} + \mathcal {X} _ {2})} e ^ {j \frac {\pi}{\lambda} (N _ {x} - 1) d _ {\mathrm {S}} \left(- f _ {1} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}} + f _ {1} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right)\right) \right.} \frac {\sin \frac {\pi}{\lambda} N _ {x} d _ {\mathrm {S}} \left(- f _ {1} \left(\theta_ {k , b} ^ {\mathrm {L} , \mathrm {e}}, \phi_ {k , b} ^ {\mathrm {L} , \mathrm {e}}\right) + f _ {1} \left(\theta_ {k , t , b} ^ {\mathrm {L}}, \phi_ {k , t , b} ^ {\mathrm {L}}\right)\right)}{\sin \frac {\pi}{\lambda} d _ {\mathrm {S}} \left(- f _ {1} \left(\theta_ {k , b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k , b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {1} \left(\theta_ {k , t , b} ^ {\mathrm {L}}, \phi_ {k , t , b} ^ {\mathrm {L}}\right)\right)} \\ + \frac {1}{N} \sum_ {\iota} \Lambda_ {k, t, b} \eta_ {k, t, b} e ^ {j \frac {2 \pi}{\lambda} (\mathcal {X} _ {3} + \mathcal {X} _ {4})} e ^ {j \frac {\pi}{\lambda} (N _ {y} - 1) d _ {\mathrm {S}} \left(- f _ {2} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}} + f _ {2} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right)\right)\right.} \frac {\sin \frac {\pi}{\lambda} N _ {y} d _ {\mathrm {S}} \left(- f _ {2} \left(\theta_ {k , b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k , b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {2} \left(\theta_ {k , t , b} ^ {\mathrm {L}}, \phi_ {k , t , b} ^ {\mathrm {L}}\right)\right)}{\sin \frac {\pi}{\lambda} d _ {\mathrm {S}} \left( \right.- f _ {2} \left(\theta_ {k , b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k , b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {2} \left(\theta_ {k , t , b} ^ {\mathrm {L}}, \phi_ {k , t , b} ^ {\mathrm {L}})\right)} \\ - \frac {1}{N} \sum_ {\iota} \Lambda_ {k, \iota , b} \eta_ {k, \iota , b} e ^ {j \frac {2 \pi}{\lambda} X _ {1}} - \frac {1}{N} \sum_ {\iota} \Lambda_ {k, \iota , b} \eta_ {k, \iota , b} e ^ {j \frac {2 \pi}{\lambda} X _ {2}} - \frac {1}{N} \sum_ {\iota} \Lambda_ {k, \iota , b} \eta_ {k, \iota , b} e ^ {j \frac {2 \pi}{\lambda} X _ {3}} - \frac {1}{N} \sum_ {\iota} \Lambda_ {k, \iota , b} \eta_ {k, \iota , b} e ^ {j \frac {2 \pi}{\lambda} X _ {4}} \\ \stackrel {N _ {x}, N _ {y} \rightarrow \infty} {=} \left\{\begin{array}{c}\Lambda_ {k, t, b} \eta_ {k, t, b}, \text {i f} f _ {1} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) = f _ {1} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right), f _ {2} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) = f _ {2} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right),\\0, \text {e l s e .}\end{array}\right. \tag {74} \\ \end{array}
$$

given as

$$
\begin{array}{l} \mathcal {Y} _ {1} = \frac {- N _ {x} + 1}{2} d _ {\mathrm {S}} \left(f _ {1} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {1} \left(\theta_ {k, \iota , b} ^ {\mathrm {L}}, \phi_ {k, \iota , b} ^ {\mathrm {L}}\right)\right) \\ + \frac {- N _ {\mathrm {y}} + 1}{2} d _ {\mathrm {S}} \left(f _ {2} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {2} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right)\right), \tag {75} \\ \end{array}
$$

$$
\begin{array}{l} \mathcal {Y} _ {2} = \frac {- N _ {x} + 1}{2} d _ {\mathrm {S}} \left(f _ {1} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {1} \left(\theta_ {k, \iota , b} ^ {\mathrm {L}}, \phi_ {k, \iota , b} ^ {\mathrm {L}}\right)\right) \\ + \frac {N _ {y} - 1}{2} d _ {\mathrm {S}} \left(f _ {2} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {2} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right)\right), \tag {76} \\ \end{array}
$$

$$
\begin{array}{l} \mathcal {Y} _ {3} = \frac {N _ {x} - 1}{2} d _ {\mathrm {S}} \left(f _ {1} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {1} \left(\theta_ {k, \iota , b} ^ {\mathrm {L}}, \phi_ {k, \iota , b} ^ {\mathrm {L}}\right)\right) \\ + \frac {- N _ {y} + 1}{2} d _ {\mathrm {S}} \left(f _ {2} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {2} \left(\theta_ {k, \iota , b} ^ {\mathrm {L}}, \phi_ {k, \iota , b} ^ {\mathrm {L}}\right)\right), \tag {77} \\ \end{array}
$$

$$
\begin{array}{l} \mathcal {Y} _ {4} = \frac {N _ {x} - 1}{2} d _ {\mathrm {S}} \left(f _ {1} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {1} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right)\right) \\ + \frac {N _ {y} - 1}{2} d _ {\mathrm {S}} \left(f _ {2} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {2} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right)\right). \tag {78} \\ \end{array}
$$

Similar to Eq. (74) and Eq. (79), we have

$$
\begin{array}{l} \lim  _ {N \rightarrow \infty} \frac {1}{N} \sum_ {n _ {x}, n _ {y}} e ^ {- j \frac {2 \pi}{\lambda} \mathbf {f} ^ {T} \left(\theta_ {k, b} ^ {\mathrm {L , e}}, \phi_ {k, b} ^ {\mathrm {L , e}}\right) \mathbf {r} _ {b, n _ {x}, n _ {y}} ^ {\mathrm {S , L}}} z _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S}} = 0, (80) \\ \lim  _ {N \rightarrow \infty} \frac {1}{N} \sum_ {n _ {x}, n _ {y}} e ^ {- j \frac {\pi}{\lambda} \mathbf {f} ^ {T} \left(\theta_ {k, b} ^ {\mathrm {L , e}}, \phi_ {k, b} ^ {\mathrm {L , e}}\right) \mathbf {r} _ {b, n _ {x}, n _ {y}} ^ {\mathrm {S , L}}} z _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S} *} \left(s _ {k, b, n _ {x}, n _ {y}} ^ {\text {r e f}}\right) ^ {2} = 0. (81) \\ \end{array}
$$

In conclusion, according to Eq. (74), Eq. (79), Eq. (80) and Eq. (81), we have

$$
\tilde {\mathcal {H}} _ {k, b} = \left\{ \begin{array}{l l} \Lambda_ {k, t, b} \eta_ {k, t, b}, & \text {i f} f _ {1} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) = f _ {1} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right) \text {a n d} \\ & f _ {2} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) = f _ {2} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right), \\ 0, & \text {e l s e .} \end{array} \right. \tag {82}
$$

The proof of Theorem 1 is completed.

APPENDIX B PROOF OF LEMMA 2

Denote  $(\theta_0^{\mathrm{L}},\phi_0^{\mathrm{L}})$  and  $(\theta^{\mathrm{L,e}},\phi^{\mathrm{L,e}})$  as the true angles and the estimated angles in the local coordinate system, which satisfies  $-\pi < \theta_0^{\mathrm{L}},\theta^{\mathrm{L,e}}\leq \pi$  and  $-\frac{\pi}{2} <  \phi_0^{\mathrm{L}},\phi^{\mathrm{L,e}}\leq \frac{\pi}{2}$ . The ambiguity function is given as

$$
\begin{array}{l} \left| \frac {1}{N} \sum_ {n _ {x}, n _ {y}} e ^ {- j \frac {2 \pi}{\lambda} \mathbf {f} ^ {T} \left(\theta^ {\mathrm {L}, \mathrm {c}}, \phi^ {\mathrm {L}, \mathrm {c}}\right) \mathbf {r} _ {b, n _ {x}, n _ {y}} ^ {\mathrm {S}, \mathrm {L}}} e ^ {j \frac {2 \pi}{\lambda} \mathbf {f} ^ {T} \left(\theta_ {0} ^ {\mathrm {L}}, \phi_ {0} ^ {\mathrm {L}}\right) \mathbf {r} _ {b, n _ {x}, n _ {y}} ^ {\mathrm {S}, \mathrm {L}}} \right| \\ = \left| \frac {1}{N} \sum_ {n _ {x}} e ^ {j \frac {2 \pi}{\lambda} \left[ n _ {x} d _ {S} \left(- f _ {1} \left(\theta^ {\mathrm {L}, \mathrm {e}}, \phi^ {\mathrm {L}, \mathrm {e}}\right) + f _ {1} \left(\theta_ {0} ^ {\mathrm {L}}, \phi_ {0} ^ {\mathrm {L}}\right)\right) + \mathcal {Z} _ {1} + \mathcal {Z} _ {2} \right]} \right. \\ + \frac {1}{N} \sum_ {n _ {y}} e ^ {j \frac {2 \pi}{\lambda} \left[ n _ {y} d _ {S} \left(- f _ {2} \left(\theta^ {\mathrm {L}, \mathrm {e}}, \phi^ {\mathrm {L}, \mathrm {e}}\right) + f _ {2} \left(\theta_ {0} ^ {\mathrm {L}}, \phi_ {0} ^ {\mathrm {L}}\right)\right) + \mathcal {Z} _ {3} + \mathcal {Z} _ {4} \right]} \\ \left. - \frac {1}{N} e ^ {j \frac {2 \pi}{\lambda}} \mathcal {Z} _ {1} - \frac {1}{N} e ^ {j \frac {2 \pi}{\lambda}} \mathcal {Z} _ {2} - \frac {1}{N} e ^ {j \frac {2 \pi}{\lambda}} \mathcal {Z} _ {3} - \frac {1}{N} e ^ {j \frac {2 \pi}{\lambda}} \mathcal {Z} _ {4} \right|. \tag {83} \\ \end{array}
$$

$$
\begin{array}{l} \frac {1}{N} \sum_ {n _ {x}, n _ {y}} e ^ {- j \frac {2 \pi}{\lambda} \mathbf {f} ^ {T} (\theta_ {k, b} ^ {\mathrm {L , e}}, \phi_ {k, b} ^ {\mathrm {L , e}}) \mathbf {r} _ {b, n _ {x}, n _ {y}} ^ {\mathrm {S , L}}} h _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S , L *}} \left(s _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {r e f}}\right) ^ {2} \\ = \frac {1}{N} \sum_ {n _ {x}, n _ {y}} e ^ {j \frac {2 \pi}{\lambda} \mathbf {f} ^ {T} \left(\theta_ {k, b} ^ {\mathrm {L}. e}, \phi_ {k, b} ^ {\mathrm {L}. e}\right) \mathbf {r} _ {b, n _ {x}, n _ {y}} ^ {\mathrm {S}, \mathrm {L}}} h _ {k, b, n _ {x}, n _ {y}} ^ {\mathrm {S}, \mathrm {L} *} e ^ {j 2 \pi \chi_ {k, b, n _ {x}, n _ {y}}} \\ = \frac {1}{N} \sum_ {\iota} \Lambda_ {k, \iota , b} \eta_ {k, t, b} \sum_ {n _ {x}} e ^ {- j \frac {2 \pi}{\lambda} \left[ n _ {x} d _ {S} \left(f _ {1} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {1} \left(\theta_ {k, x, b} ^ {\mathrm {L}}, \phi_ {k, x, b} ^ {\mathrm {L}}\right)\right) + \mathcal {Y} _ {1} \right] e ^ {j 2 \pi \chi_ {k, b, n _ {x}, 0}} \\ + \frac {1}{N} \sum_ {\iota} \Lambda_ {k, t, b} \eta_ {k, t, b} \sum_ {n _ {x}} e ^ {- j \frac {2 \pi}{\lambda} \left[ n _ {x} d _ {S} \left(f _ {1} \left(\theta_ {k, b} ^ {\mathrm {L}, e}, \phi_ {k, b} ^ {\mathrm {L}, e}\right) + f _ {1} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right)\right) + \mathcal {Y} _ {2} \right] e ^ {j 2 \pi \chi_ {k, b, n _ {x}, N _ {y} - 1}} \\ + \frac {1}{N} \sum_ {\iota} \Lambda_ {k, \iota , b} \eta_ {k, \iota , b} \sum_ {n _ {\mathrm {y}}} e ^ {- j \frac {2 \pi}{\lambda} \left[ \mathcal {Y} _ {3} + n _ {\mathrm {y}} d _ {\mathrm {S}} \left(f _ {2} \left(\theta_ {k, b} ^ {\mathrm {L}, \mathrm {e}}, \phi_ {k, b} ^ {\mathrm {L}, \mathrm {e}}\right) + f _ {2} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right)\right) \right]} e ^ {j 2 \pi \chi_ {k, b, 0, n _ {\mathrm {y}}}} \\ + \frac {1}{N} \sum_ {\iota} \Lambda_ {k, \iota , b} \eta_ {k, \iota , b} \sum_ {n _ {y}} e ^ {- j \frac {2 \pi}{\lambda} \left[ \mathcal {Y} _ {4} + n _ {y} d _ {S} \left(f _ {2} \left(\theta_ {k, b} ^ {\mathrm {L}, e}, \phi_ {k, b} ^ {\mathrm {L}, e}\right) + f _ {2} \left(\theta_ {k, t, b} ^ {\mathrm {L}}, \phi_ {k, t, b} ^ {\mathrm {L}}\right)\right) \right]} e ^ {j 2 \pi \chi_ {k, b, 0, N _ {y} - 1}} \\ - \frac {1}{N} \sum_ {\iota} \Lambda_ {k, \iota , b} \eta_ {k, \iota , b} e ^ {- j \frac {2 \pi}{\lambda} \mathcal {Y} _ {1}} e ^ {j 2 \pi \chi_ {k, b, 0, 0}} - \frac {1}{N} \sum_ {\iota} \Lambda_ {k, \iota , b} \eta_ {k, \iota , b} e ^ {- j \frac {2 \pi}{\lambda} \mathcal {Y} _ {2}} e ^ {j 2 \pi \chi_ {k, b, 0, N _ {y} - 1}} \\ - \sum_ {\iota} \Lambda_ {k, \iota , b} \eta_ {k, \iota , b} e ^ {- j \frac {2 \pi}{\lambda} \mathcal {Y} _ {3}} e ^ {j 2 \pi \chi_ {k, b, N _ {x} - 1, 0}} - \frac {1}{N} \sum_ {\iota} \Lambda_ {k, \iota , b} \eta_ {k, \iota , b} e ^ {- j \frac {2 \pi}{\lambda} \mathcal {Y} _ {4}} e ^ {j 2 \pi \chi_ {k, b, N _ {x} - 1, N _ {y} - 1}} \\ \stackrel {N \rightarrow \infty} {=} \frac {1}{N} \sum_ {\iota} \Lambda_ {k, \iota , b} \eta_ {k, \iota , b} \mathbb {E} \left[ e ^ {- j \frac {2 \pi}{\lambda} \mathcal {Y} _ {4}} e ^ {j 2 \pi \chi_ {k, b, n _ {x} - 1, n _ {y} - 1}} \right] \\ = 0. \tag {79} \\ \end{array}
$$

where we have

$$
\begin{array}{l} \mathcal {Z} _ {1} = \frac {- N _ {x} + 1}{2} d _ {\mathrm {S}} \left(- f _ {1} \left(\theta^ {\mathrm {L}, \mathrm {e}}, \phi^ {\mathrm {L}, \mathrm {e}}\right) + f _ {1} \left(\theta_ {0} ^ {\mathrm {L}}, \phi_ {0} ^ {\mathrm {L}}\right)\right) \\ + \frac {- N _ {y} + 1}{2} d _ {\mathrm {S}} \left(- f _ {2} \left(\theta^ {\mathrm {L}, \mathrm {e}}, \phi^ {\mathrm {L}, \mathrm {e}}\right) + f _ {2} \left(\theta_ {0} ^ {\mathrm {L}}, \phi_ {0} ^ {\mathrm {L}}\right)\right), \tag {84} \\ \end{array}
$$

$$
\begin{array}{l} \mathcal {Z} _ {2} = \frac {- N _ {x} + 1}{2} d _ {\mathrm {S}} \left(- f _ {1} \left(\theta^ {\mathrm {L}, \mathrm {e}}, \phi^ {\mathrm {L}, \mathrm {e}}\right) + f _ {1} \left(\theta_ {0} ^ {\mathrm {L}}, \phi_ {0} ^ {\mathrm {L}}\right)\right) \\ + \frac {N _ {\mathrm {y}} - 1}{2} d _ {\mathrm {S}} \left(- f _ {2} \left(\theta^ {\mathrm {L}, \mathrm {e}}, \phi^ {\mathrm {L}, \mathrm {e}}\right) + f _ {2} \left(\theta_ {0} ^ {\mathrm {L}}, \phi_ {0} ^ {\mathrm {L}}\right)\right), \tag {85} \\ \end{array}
$$

$$
\begin{array}{l} \mathcal {Z} _ {3} = \frac {N _ {x} - 1}{2} d _ {\mathrm {S}} \left(- f _ {1} \left(\theta^ {\mathrm {L}, \mathrm {e}}, \phi^ {\mathrm {L}, \mathrm {e}}\right) + f _ {1} \left(\theta_ {0} ^ {\mathrm {L}}, \phi_ {0} ^ {\mathrm {L}}\right)\right) \\ + \frac {- N _ {y} + 1}{2} d _ {\mathrm {S}} \left(- f _ {2} \left(\theta^ {\mathrm {L}, \mathrm {e}}, \phi^ {\mathrm {L}, \mathrm {e}}\right) + f _ {2} \left(\theta_ {0} ^ {\mathrm {L}}, \phi_ {0} ^ {\mathrm {L}}\right)\right), \tag {86} \\ \end{array}
$$

$$
\begin{array}{l} \mathcal {Z} _ {4} = \frac {N _ {x} - 1}{2} d _ {\mathrm {S}} \left(- f _ {1} \left(\theta^ {\mathrm {L}, \mathrm {e}}, \phi^ {\mathrm {L}, \mathrm {e}}\right) + f _ {1} \left(\theta_ {0} ^ {\mathrm {L}}, \phi_ {0} ^ {\mathrm {L}}\right)\right) \\ + \frac {N _ {\mathrm {y}} - 1}{2} d _ {\mathrm {S}} \left(- f _ {2} \left(\theta^ {\mathrm {L}, \mathrm {e}}, \phi^ {\mathrm {L}, \mathrm {e}}\right) + f _ {2} \left(\theta_ {0} ^ {\mathrm {L}}, \phi_ {0} ^ {\mathrm {L}}\right)\right). \tag {87} \\ \end{array}
$$

The ambiguity function (83) achieves the maximal value if and only if  $-f_{1}(\theta^{\mathrm{L,e}},\phi^{\mathrm{L,e}}) + f_{1}(\theta_{0}^{\mathrm{L}},\phi_{0}^{\mathrm{L}}) = p^{\prime}\frac{\lambda}{ds}$  and  $-f_{2}(\theta^{\mathrm{L,e}},\phi^{\mathrm{L,e}}) +$ $f_{2}(\theta_{0}^{\mathrm{L}},\phi_{0}^{\mathrm{L}}) = q^{\prime}\frac{\lambda}{ds}$ , where  $p^{\prime},q^{\prime}\in \mathbb{Z}$ . Observing Eq. (83), if we want to correctly estimate the direction vector  $\mathbf{f}(\theta_0^{\mathrm{L}},\phi_0^{\mathrm{L}})$  by adopting the FFT-based detection method and avoid the angle aliasing, we have  $-1 < p^{\prime} < 1$  and  $-1 < q^{\prime} < 1$ , while the following conditions should be satisfied

$$
\left| \frac {d _ {\mathrm {S}} \left(- f _ {1} \left(\theta^ {\mathrm {L} , \mathrm {e}} , \phi^ {\mathrm {L} , \mathrm {e}}\right) + f _ {1} \left(\theta_ {0} ^ {\mathrm {L}} , \phi_ {0} ^ {\mathrm {L}}\right)\right)}{\lambda} \right| <   1, \tag {88}
$$

$$
\left| \frac {d _ {\mathbb {S}} \left(- f _ {2} \left(\theta^ {\mathrm {L} , \mathrm {e}} , \phi^ {\mathrm {L} , \mathrm {e}}\right) + f _ {2} \left(\theta_ {0} ^ {\mathrm {L}} , \phi_ {0} ^ {\mathrm {L}}\right)\right)}{\lambda} \right| <   1. \tag {89}
$$

Since we have  $\left| - f_{1}(\theta^{\mathrm{L},\mathrm{e}},\phi^{\mathrm{L},\mathrm{e}}) + f_{1}(\theta_{0}^{\mathrm{L}},\phi_{0}^{\mathrm{L}})\right| < 2$  and  $\left| - f_{2}(\theta^{\mathrm{L},\mathrm{e}},\phi^{\mathrm{L},\mathrm{e}}) + f_{2}(\theta_{0}^{\mathrm{L}},\phi_{0}^{\mathrm{L}})\right| < 2$ , the constraints (88) and (89) can be reformulated as

$$
d _ {\mathrm {S}} <   \frac {\lambda}{\left| - f _ {1} \left(\theta^ {\mathrm {L}, \mathrm {e}}, \phi^ {\mathrm {L}, \mathrm {e}}\right) + f _ {1} \left(\theta_ {0} ^ {\mathrm {L}}, \phi_ {0} ^ {\mathrm {L}}\right) \right|}, \tag {90}
$$

$$
d _ {\mathrm {S}} <   \frac {\lambda}{\left| - f _ {2} \left(\theta^ {\mathrm {L} , \mathrm {e}}, \phi^ {\mathrm {L} , \mathrm {e}}\right) + f _ {2} \left(\theta_ {0} ^ {\mathrm {L}}, \phi_ {0} ^ {\mathrm {L}}\right) \right|}. \tag {91}
$$

Moreover, we have

$$
\frac {\lambda}{\left| - f _ {1} \left(\theta^ {\mathrm {L , e}} , \phi^ {\mathrm {L , e}}\right) + f _ {1} \left(\theta_ {0} ^ {\mathrm {L}}, \phi_ {0} ^ {\mathrm {L}}\right) \right|} > \frac {\lambda}{2}, \tag {92}
$$

$$
\frac {\lambda}{\left| - f _ {2} \left(\theta^ {\mathrm {L , e}} , \phi^ {\mathrm {L , e}}\right) + f _ {2} \left(\theta_ {0} ^ {\mathrm {L}}, \phi_ {0} ^ {\mathrm {L}}\right) \right|} > \frac {\lambda}{2}. \tag {93}
$$

Finally, the constraints (88)-(93), we have

$$
d _ {\mathrm {S}} \leq \frac {\lambda}{2}. \tag {94}
$$

Therefore, arrange the sensing element with the half-wavelength spacing is sufficient and laying the sensing elements with sub-wavelength is not required. The proof of Theorem 2 is completed.

# REFERENCES

[1] K. Ahmadi and W. A. Serdijn, “Advancements in laser and led-based optical wireless power transfer for IoT applications: A comprehensive review,” IEEE Internet of Things Journal, pp. 1–1, 2025.  
[2] Y. Zheng, J. Hu, Y. Zhao, and K. Yang, "Average age of sensing in wireless powered sensor networks," IEEE Transactions on Wireless Communications, vol. 23, no. 8, pp. 9265-9281, 2024.

[3] X. Li, Z. Han, G. Zhu, Y. Shi, J. Xu, Y. Gong, Q. Zhang, K. Huang, and K. B. Letaief, “Integrating sensing, communication, and power transfer: From theory to practice,” IEEE Communications Magazine, vol. 62, no. 9, pp. 122–127, 2024.  
[4] R. Ma, J. Tang, X. Y. Zhang, W. Feng, D. K. C. So, C.-B. Chae, K.-K. Wong, and J. A. Chambers, "Simultaneous wireless information and power transfer in IoT-based scenarios: Architectures, challenges, and prototype validation," IEEE Wireless Communications, vol. 31, no. 5, pp. 272-278, 2024.  
[5] G. Kwon, H. Park, and M. Z. Win, "Joint beamforming and power splitting for wideband millimeter wave swipt systems," IEEE Journal of Selected Topics in Signal Processing, vol. 15, no. 5, pp. 1211-1227, 2021.  
[6] J. S. Herd and M. D. Conway, "The evolution to modern phased array architectures," Proceedings of the IEEE, vol. 104, no. 3, pp. 519-529, 2016.  
[7] R. Deng, B. Di, H. Zhang, Y. Tan, and L. Song, "Reconfigurable holographic surface-enabled multi-user wireless communications: Amplitude-controlled holographic beamforming," IEEE Transactions on Wireless Communications, vol. 21, no. 8, pp. 6003-6017, 2022.  
[8] X. Shao, R. Zhang, Q. Jiang, and R. Schober, "6d movable antenna enhanced wireless network via discrete position and rotation optimization," ArXiv, vol. abs/2403.17122, 2024. [Online]. Available: https://api_semanticscholar.org/CorpusID:268691662  
[9] L. R. Varshney, “Transporting information and energy simultaneously,” in 2008 IEEE International Symposium on Information Theory, 2008, pp. 1612–1616.  
[10] B. Clerckx, "Wireless information and power transfer: Nonlinearity, waveform design, and rate-energy tradeoff," IEEE Transactions on Signal Processing, vol. 66, no. 4, pp. 847-862, 2018.  
[11] S. Na, F. Hu, Z. Ling, H. Dong, and X. Yao, "Joint frequency-time allocation and phase-shift optimization in intelligent reflecting surface assisted multigroup wpcn," IEEE Internet of Things Journal, vol. 12, no. 17, pp. 35249-35260, 2025.  
[12] R. Deng, B. Di, H. Zhang, and L. Song, "Hdma: Holographic-pattern division multiple access," IEEE Journal on Selected Areas in Communications, vol. 40, no. 4, pp. 1317-1332, 2022.  
[13] B. Di, “Reconfigurable holographic metasurface aided wideband ofdm communications against beam squint,” IEEE Transactions on Vehicular Technology, vol. 70, no. 5, pp. 5099-5103, 2021.  
[14] A. Azar Bahram, O. L. A. López, and M. Latva-Aho, "Waveform optimization and beam focusing for near-field wireless power transfer with dynamic metasurface antennas and non-linear energy harvesters," 2024. [Online]. Available: https://arxiv.org/abs/2307.01081  
[15] Q. Huang, J. Hu, Y. Zhao, and K. Yang, "Holographic integrated data and energy transfer," IEEE Transactions on Wireless Communications, vol. 23, no. 12, pp. 18987-19002, 2024.  
[16] K.-K. Wong, A. Shojaeifard, K.-F. Tong, and Y. Zhang, "Fluid antenna systems," IEEE Transactions on Wireless Communications, vol. 20, no. 3, pp. 1950-1962, 2021.  
[17] X. Lin, Y. Zhao, H. Yang, J. Hu, and K.-K. Wong, "Fluid antenna multiple access assisted integrated data and energy transfer: Outage and multiplexing gain analysis," IEEE Transactions on Wireless Communications, pp. 1-1, 2025.  
[18] X. Shao, Q. Jiang, and R. Zhang, "6d movable antenna based on user distribution: Modeling and optimization," IEEE Transactions on Wireless Communications, pp. 1-1, 2024.  
[19] W. Zhonglun, Z. Yizhe, H. Jie, and Y. Kun, "6dma-assisted integrated data and energy transfer: Joint spatial orientation and beamforming design," in 2025 IEEE 102th Vehicular Technology Conference (VTC2025-Fall), 2025, pp. 1-1.  
[20] X. Shao, R. Zhang, Q. Jiang, J. Park, T. Q. S. Quek, and R. Schober, "Distributed channel estimation and optimization for 6d movable antenna: Unveiling directional sparsity," 2024. [Online]. Available: https://api-semanticscholar.org/CorpusID:272880808  
[21] H. Li, S. Shen, and B. Clerckx, "Beyond diagonal reconfigurable intelligent surfaces: A multi-sector mode enabling highly directional full-space wireless coverage," IEEE Journal on Selected Areas in Communications, vol. 41, no. 8, pp. 2446-2460, 2023.  
[22] H. Zhang, H. Zhang, B. Di, M. D. Renzo, Z. Han, H. V. Poor, and L. Song, "Holographic integrated sensing and communication," IEEE Journal on Selected Areas in Communications, vol. 40, no. 7, pp. 2114-2130, 2022.  
[23] J. Hu, Y. Zheng, and K. Yang, "Multi-domain resource scheduling for full-duplex aided wireless powered communication network," IEEE Transactions on Vehicular Technology, vol. 71, no. 10, pp. 10849-10862, 2022.

[24] R. Deng, Y. Zhang, H. Zhang, B. Di, H. Zhang, H. V. Poor, and L. Song, "Reconfigurable holographic surfaces for ultra-massive mimo in 6g: Practical design, optimization and implementation," IEEE Journal on Selected Areas in Communications, vol. 41, no. 8, pp. 2367-2379, 2023.  
[25] R. Deng, B. Di, H. Zhang, H. V. Poor, and L. Song, "Holographic mimo forleo satellite communications aided by reconfigurable holographic surfaces," IEEE Journal on Selected Areas in Communications, vol. 40, no.10,pp.3071-3085,2022.  
[26] Y. Zhao, Y. Wu, J. Hu, and K. Yang, "A general analysis and optimization framework of time index modulation for integrated data and energy transfer," IEEE Transactions on Wireless Communications, vol. 22, no. 6, pp. 3657-3670, 2023.  
[27] S. Liu, X. Yu, Z. Gao, J. Xu, D. W. K. Ng, and S. Cui, "Sensing-enhanced channel estimation for near-field xl-mimo systems," IEEE Journal on Selected Areas in Communications, pp. 1-1, 2025.  
[28] Interior-Point Methods for Nonlinear Programming. New York, NY: Springer New York, 2006, pp. 563-597. [Online]. Available: https://doi.org/10.1007/978-0-387-40065-5 19  
[29] Q. Yang, A. Guerra, F. Guidi, N. Shlezinger, H. Zhang, D. Dardari, B. Wang, and Y. C. Eldar, "Near-field localization with dynamic metasurface antennas," in ICASSP 2023 - 2023 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), 2023, pp. 1-5.  
[30] M. Chen, Y. Miao, Y. Hao, and K. Hwang, "Narrow band internet of things," IEEE Access, vol. 5, pp. 20557-20577, 2017.  
[31] C. Hu, L. Dai, S. Han, and X. Wang, "Two-timescale channel estimation for reconfigurable intelligent surface aided wireless communications," IEEE Transactions on Communications, vol. 69, no. 11, pp. 7736-7747, 2021.