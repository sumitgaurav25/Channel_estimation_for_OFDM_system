<b><h2>TOPIC - Channel Estimation for OFDM Systems Using MMSE and LS Algorithms</h2></b>

<b><h3>1. INTRODUCTION</h3></b>

# Channel Estimation using LS and MMSE for MIMO-OFDM Systems

## 📖 Overview
This repository provides a study and simulation of **Least Squares (LS)** and **Minimum Mean Square Error (MMSE)** channel estimation techniques applied to **MIMO-OFDM systems**.  
The research investigates pilot-assisted block-type training for estimating channel state information (CSI), focusing on the trade-off between **bit error rate (BER)**, **signal-to-noise ratio (SNR)**, and computational complexity.  

The work is based on the IEEE paper:  
**"Channel Estimation using LS and MMSE Channel Estimation Techniques for MIMO-OFDM Systems" (2022).**

---

## ✨ Key Contributions
- Implementation of **pilot-aided channel estimation** for OFDM.
- Comparative study of **LS vs. MMSE estimators**.
- Analysis of **BER vs. SNR** under different modulation schemes (BPSK, QPSK, 16-QAM, 32-QAM, 64-QAM).
- Evaluation of **mean square error (MSE)** and **symbol error rate (SER)** performance.
- Application scenarios for **MIMO wireless systems, LTE/LTE-A, and 6G candidates**.

---

## 🏗 System Model
The baseband OFDM system includes:
- **IFFT/FFT blocks** for modulation/demodulation  
- **Cyclic prefix insertion/removal**  
- **Pilot-based channel estimation**  
- **Constellation mapping/demapping**  

### Basic OFDM Structure
![OFDM Structure](images/ofdm_structure.png)

---

## 📊 Simulation Results

### Pilot Patterns
Block-type vs. comb-type pilot insertions are tested:
![Pilot Patterns](images/pilot_patterns.png)

### QPSK Transmission Example
- Transmitted image with QPSK modulation  
- Received image with noise and channel distortions  

| Transmitted | Received |
|-------------|----------|
| ![QPSK Tx](images/qpsk_tx.png) | ![QPSK Rx](images/qpsk_rx.png) |

### Without Channel Estimation
![No Estimation](images/no_estimation.png)

### BER vs SNR
Comparison of LS and MMSE estimators:
![BER vs SNR](images/ber_snr.png)

### MSE and SER Comparison
![MSE vs SNR](images/mse_snr.png)  
![SER vs SNR](images/ser_snr.png)

---

## 🔍 Observations
- **LS Estimator**: Simple, lower computational load, but higher MSE.  
- **MMSE Estimator**: Better accuracy, reduced BER and SER, but computationally more complex.  
- **Higher-order QAM schemes**: Increase data rates but introduce greater sensitivity to noise and fading.  
- **Pilot design** is critical for achieving reliable estimation.  

---

## 🚀 Future Work
- Integration of **adaptive algorithms (LMS, RLS, Kalman filters)**.  
- Use of **compressive sensing** for sparse channel estimation.  
- Application of **deep learning-based channel estimators** for massive MIMO and mmWave systems.  
- Multi-band evaluation for **5G/6G heterogeneous networks**.  

---

## 📎 Reference
Ahmed, A.S., Khaleel, A.M., Hamdi, M.M., Fathy, M., Abood, M.S., & Khaleefah, S.H.  
*Channel Estimation using LS and MMSE Channel Estimation Techniques for MIMO-OFDM Systems.*  
IEEE, 2022. DOI: [10.1109/HORA55278.2022.9799887](https://doi.org/10.1109/HORA55278.2022.9799887)

---

