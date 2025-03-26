# 🛰️ Non-contact Unrestrained Activity Analysis using Radar Sensor

> 비접촉 · 무구속 방식의 FMCW 레이더 기반 스트레스 정량화 시스템  
> 최태훈, 숭실대학교 전자정보공학부 졸업논문 (2024)

## 📌 Introduction

This project proposes a **non-contact, unrestrained method** of measuring stress levels using an FMCW radar sensor. Traditional approaches (e.g., ECG, PPG sensors) require body attachment, which causes discomfort and limits their use in daily environments.  
Instead, this study utilizes **FMCW radar to extract heart rate and heart rate variability (HRV)** metrics remotely and estimates stress scores based on those values.

## 🧠 Key Features

- 📶 **FMCW Radar-Based Measurement**  
  Utilizes beat frequency analysis to extract micro-vibrations caused by heartbeats.

- 🧮 **Custom Stress Scoring Algorithm**  
  Calculates normalized stress scores (0–100) using:
  - Heart Rate (HR)
  - SDNN (Standard Deviation of NN intervals)
  - RMSSD (Root Mean Square of Successive Differences)

- 🧹 **Signal Preprocessing**  
  - High-pass + band-pass filtering  
  - Moving Average filter for outlier removal  
  - FFT and phase difference analysis for fine-grained measurement

- 🧪 **Experiment Validation**  
  Compared radar-derived stress scores with **Samsung Galaxy Watch** measurements across 3 stress-inducing scenarios.

## 📊 Experimental Scenarios

1. **Meditative State**  
2. **Cognitive Stress** (solving aptitude test)  
3. **Physical Stress** (holding push-up position)

### 📉 Result Highlights

| Condition         | Radar Stress Score | Galaxy Watch Score |
|------------------|--------------------|--------------------|
| Meditation        | 50–56              | (similar trend)    |
| Aptitude Test     | 50–56              | (similar trend)    |
| Physical Overload | 79–84              | (similar trend)    |

→ Radar algorithm showed **consistent stress trends** with a commercial device.

## 🛠️ Technologies

- **Hardware**: MOD620 FMCW Radar
- **Signal Processing**: Python / MATLAB (assumed for algorithm implementation)
- **Metrics**:
  - HR: Beats per Minute
  - SDNN / RMSSD: HRV indicators
- **Data Analysis**: FFT, filtering, normalization

## 💡 Future Work

- Broader dataset collection in various environments
- Real-time implementation with embedded systems
- Extension to multi-user or dynamic movement scenarios

## 📚 References

- M.Jankiraman, *FMCW Radar Design*
- Rémi Grisot et al., *Monitoring of Heart Movements Using an FMCW Radar*
- Xiangyu Han et al., *A Real-Time Evaluation Algorithm for Noncontact HRV Monitoring*
- 기타 논문 참고는 원문 `참고문헌` 참조

