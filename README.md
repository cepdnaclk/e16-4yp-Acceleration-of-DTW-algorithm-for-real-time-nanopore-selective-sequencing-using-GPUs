___
# Acceleration of Dynamic Time Warping algorithm for real-time nanopore selective sequencing using GPUs
___

## Aim of the project: 

The proposed project aims to optimize the Dynamic Time Warping (DTW) algorithm and accelerate it using Graphics Processing Units (GPUs), So that algorithm can be executed 
in a GPU-equipped laptop or a GPU-equipped embedded device like NVIDIA Jetson, rather 
than connecting to a massive server. 

## Background of the project: 

Nanopore sequencing is a technique that is used to real-time analyze long DNA and 
RNA fragments. Its work principle involves monitoring changes in an electric current passing 
through a protein nanopore (a nanopore is a nano-size pore or cavity) as a nucleic acid. The 
resulting signals can be used to identify specific DNA or RNA fragments. One of the most 
advantageous factors of nanopore sequencing is that it can be used to get real-time data 
sequencing which provides immediate access to the results.

The modern nanopore sequencers offer selective sequencing capability which allows 
one to select a DNA strand sequence programmatically in real-time by using the sequencer. 
This method is becoming more popular due to its real-time data access, cost when compared 
to other existing sequencing methods, applicability in rapid diagnostics, etc… An algorithm is 
used to determine whether a particular piece of DNA belongs to the expected region or not. 
This algorithm is called Dynamic Time Warping (DTW) algorithm. 
The dynamic Time Warping (DTW) algorithm is a programming algorithm that is 
used to identify the optimal alignment of two sequences of values. The main idea of the DTW 
is to measure the distance between the matching of similar elements of a time series. This 
algorithm is a highly computationally demanding algorithm that requires high processing 
power. 

There is a wide range of devices available to analyze DNA and RNA using nanopore 
technologies. Because of the high computational demand of the DTW algorithm, portable 
MinION sequencers are required to connect to a large server to do the analyses. Consequently,
it will reduce the widespread adaptation of selective sequencing in a portable setting.
The main objective of the project is to optimize the DTW algorithm and accelerate it 
using GPUs so that it will be able to run in GPU-equipped devices like laptops or NVIDIA 
Jetson. A Computer Unified Device Architecture (CUDA) which is developed by NVIDIA 
will be used to optimize the DTW algorithm. CUDA is a parallel computing platform and an 
application programming interface that allows the use of GPUs as general-purpose processing 
units.
