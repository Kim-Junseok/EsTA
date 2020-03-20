# Estimation of Timing Advance (EsTA)

## 1. Data Files
### - Location and file name
      - data-#(subcarrier spacing)/node-#.txt
### - Column information (unit)
      - 1: time (ns)
      - 2: Position.x (m), 3: Position.y (m), 4: Position.z (m)
      - 5: Velocity.x (m/s), 6: Velocity.y (m/s), 7: Velocity.z (m/s)
      - 8: 3D distance between UE and gNB (m)
      - 9: Channel Env. (LOS:1, NLOS:0)
      - 10: Reference Signal Received Power (dBm)
## 2. Python Files
### - EsTA_ML_#(subcarrier spacing)kHz.ipynb
      - Develop a Deep Neural Network model (DNN) using training and validation set
      - Test TA estimation accuracy using test set
      - Accuracy (in 'accuracy' function)
            - DNN: comment out the codes
            - DNN with timing error limit values: default
## 3. MATLAB Files
### - SumoTraceIntp.m
      - Read SUMO traces from the MATLAB data file named 'berlin-trace-100.mat'
      - Organize the traces by each node
      - Interpolate the traces in 5 ms increments
### - SumoTraceRsrp.m
      - Calculate RSRP and TA command
      - Store a set of information to the Data Files (for each node)
      - Use CalcDist.m, DetermineLos.m, and UeLocation.m functions
### - SumoTraceTraining.m
      - Divide data set into train, validation, and test sets
### - SumoTraceTest.m
      - Test based on the probability characteristics (mean, standard deviation) of RSRP values according to each TA command
