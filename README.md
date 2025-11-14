# OFDM-TxRx-Baseband-Project

## 📘 Overview
This project implements a clean **OFDM transmitter–receiver (Tx/Rx)** simulation in **MATLAB**.  
The system includes QAM modulation, IFFT/FFT, cyclic prefix (CP), Rayleigh fading channel, equalization, and performance visualization.

A compact RF/DSP project by **Brian Rono**.

---

## ⚙️ Features
- Random bit generation  
- 16-QAM mapping  
- Serial–parallel → **IFFT** → CP insertion  
- Rayleigh fading frequency-domain channel  
- Equalization in frequency domain  
- Simple BER computation  
- Clean time-domain waveform visualization  

---

## 🖼 Included Output
- `q1.png` → OFDM time-domain waveform (first 500 samples)

---

## ▶️ How to Run
1. Open MATLAB  
2. Load `ofdm_baseband_project.m`  
3. Press **Run**  
4. View the generated OFDM waveform  

---

## 🔮 Future Work
- Add AWGN + detailed multipath fading (EPA/EVA/ETU models)  
- Pilot-based channel estimation (LS/MMSE)  
- BER vs SNR performance curves  
- Adaptive modulation (4QAM/16QAM/64QAM/256QAM)  
- Timing + CFO synchronization (Schmidl-Cox)  
- Full bit→bit Tx/Rx chain with constellation plots  
- Extend to **MIMO-OFDM** (2×2 or 4×4)  

---

## 👤 Author
**Brian Rono**  
RF & Wireless • DSP • MATLAB  
🔗 https://github.com/ronobrian-eng
