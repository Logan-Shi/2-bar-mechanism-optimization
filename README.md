# Problem
![img](media/mechanism.jpg)

For the above mechanism, find the best joint position, so that the system swing for the longest time under gravity with zero initial speed at given position. 

# Result

followings are some of the results under different damping

|![1](results/0.05-72.183.gif) C=0.05Nm/(rad/s) and optimized mechanism|![2](results/result.gif) C=0.05Nm/(rad/s) and random mechanism|
|-|-|
|![3](results/0.01-32.539.gif) **C=0.01Nm/(rad/s) and optimized mechanism**|![4](results/0.04.gif) **C=0.04Nm/(rad/s) and optimized mechanism**|

# Detailed discription

[report](media/report.pdf)

# Usage

1. run kinematic analysis for velocity equation and acceleration equation using Symbolic Toolbox
```
kinematics_analysis.m
```

2. run dynamic analysis with 4th Adams form and plot the animation under given dynamics parameters
```
dynamics_analysis.m
```

3. optimize for longest swing time with ga (genetic algorithm)
```
time_optimize.m
```
