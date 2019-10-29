# Speech Recognition System

Aim is to accurately recognize the user who is speaking on the basis of information available in the audio signal. <br/>
One obvious challenge faced was that the audio signals are not stationary and the algorithm has to be robust to noise, variation over time and different speaking rates. The problem is tackled in two stages: 
1. Extracting the audio features using **Mel-Frequency Cepstrum Coefficients (MFCC)**
2. Pattern Recognition using **LBG (Linde, Buzo and Gray)** Algorithm

