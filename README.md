# object-tracking-sample

1. Analysis:
* There are 2 sources of data: model prediction and measurements. Model prediction is based on dynamic equation, while measurements are updated in realtime. An assumption is the both model and sensor noise is Gaussian noise. However, the accurate covariance is not provided. 
* Tracking boxes are independent, therefore they can be handled seperately.

2. Solution:
The idea to identify target location is to fuse the prediction and measurements. For linear system, Kalman filter is one of the solutions.
    - Prediction: xhat & yhat is estimated based on prior estimated position with system model.
    - Measurement: x & y is assumed to be the center of the box as given.

The main steps include:
* Step 1: Rewrite dynamic equations into state space model.
    [dx1; dx2] = [-gamma 1; 0 0] * [x1; x2] + [0; 1] * u
    y = [1 0] * [x1; x2]
* Step 2: Discretise the above continuous state space model.
    [x1(k); x2(k)] = Phi * [x1(k-1); x2(k-1)] + Gamma * u(k-1)
    y = C * [x1(k-1); x2(k-1)]
* Step 3: Check if system is observable.
    Check if det(transpose[C CA]) != 0;
* Step 4: Apply Kalman filter.
    Follows formular introduced in De Silva, C. W. (2016). Sensor systems: Fundamentals and applications.
* Step 5: Iterate over a number of time steps and for each object.

Note: For convenience, the final outcome of state space and Kalman filter matrices is not fully described here. Please refer to source code "object_tracking.cc", for detailed result.

3. Discussion:
* Limitations:
    - Box ordering is not handled.
    - In estimation phase, one of the assumption is that measured target position is center of the box, which reduces the accuracy of algorithm.

* Future improvements:
    - Covariance of noises should be measured for better estimation of Kalman filter. 
    - Uncertain box ordering: a distance function such as L2 norm (d) can be applied to measure the distance between predicted position and measured position. The measurement is of an object iff the distance is less than a certain threshold. However, it might not work if there are many nearby objects.

