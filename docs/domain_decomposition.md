# Domain Decomposition

Domain decomposition splits the simulation up in 3D space and allocates each processor a region to oversee.

## Scaling

### Strong

```mermaid
xychart-beta
  title "Strong Scaling"
  x-axis "nprocs" [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]
  y-axis "average elapsed time"
  line [10.7746, 4.2124, 2.9203, 2.3514, 2.1331, 1.9305, 1.9184, 1.8239, 1.8147, 1.8457, 1.9195, 1.8777, 2.0219, 1.9451, 1.9438, 1.9596, 2.0091, 1.931, 2.0265, 2.1824, 2.192, 2.2465, 2.2823, 2.2453, 2.2316, 2.2524, 2.2411]
```

> 100k particles, 0.05 cutoff, O3 optimisation, average of 10 results

```mermaid
xychart-beta
  title "Strong Scaling"
  x-axis "nprocs" [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]
  y-axis "average speed up"
  line [1, 2.557829266, 3.689552443, 4.582206345, 5.051146219, 5.581248381, 5.616451209, 5.907451066, 5.937400121, 5.837676762, 5.613232613, 5.738190339, 5.328948019, 5.539355303, 5.543059986, 5.498367014, 5.36289881, 5.579803211, 5.316851715, 4.937041789, 4.915419708, 4.796171823, 4.720939403, 4.798735136, 4.828195017, 4.783608595, 4.807728348]
  line [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]
```

> 100k particles, 0.05 cutoff, O3 optimisation, average of 10 results

### Weak

```mermaid
xychart-beta
  title "Weak Scaling"
  x-axis "nprocs" [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]
  y-axis "average elapsed time"
  line [1.8504, 4.0399, 4.7627, 5.1429, 5.6638, 5.9113, 7.4946, 7.4222, 8.4419, 9.4464, 11.3624, 9.2301, 12.8243, 11.4252, 12.6102, 12.3567, 18.9139, 12.9938, 22.4954, 16.1379, 20.5822, 24.4191, 30.5043, 19.4974, 22.6823, 29.9772, 22.7041]
  line [1.8504, 4.0399, 4.7627, 5.1429, 5.6638, 5.9113, 7.4946, 7.4222, 8.4419, 9.4464, 9.33825, 9.2301, 10.32765, 11.4252, 12.6102, 12.3567, 12.67525, 12.9938, 14.56585, 16.1379, 20.5822, 20.2206, 19.858999999999998, 19.4974, 22.6823, 22.6932, 22.7041]
```

> 50k particles per processor, 0.05 cutoff, O3 optimisation, average of 10 results. The straighter line ignores any `nprocs` value in which the cutoff was wider than any of the domains.