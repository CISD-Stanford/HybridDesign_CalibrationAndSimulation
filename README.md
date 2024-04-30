### Calibration for 2-step hybrid design

The 2-step hybrid trial design is a Frequentist approach to borrowing real-world data for clinical trial. 
It utilized propensity score matching to ensure the equivalence assumption is met and an equivalence test to validate it.
However, the original 2-step hybrid design did not consider the type I error inflation due to the change of distribution of the overall test statistics.
This R-package provides three calibration methods to ensure the type I error is under controlled when the equivalence assumption is met.
The simulation results show that the three calibration methods improved the power while maintaining the type I error under control, given that the distributions
of the external control data are similar to that of the current control data. Approach 1 is the most conservative method among the three. It has a larger control area, 
which makes it capable of controlling the type I error even if a slight mean difference exists between the current control arm and the external control arm. 
On the contrary, Approaches 2 and 3 can only control the type I error if the mean difference is 0. Therefore, if there is a strong belief that the real-world data shares
the same distribution as the randomized control arm, Approach 2 or Approach 3 will be preferred. However, if one wants to control the type I error more strictly, 
it is recommended to use approach 1 in borrowing real-world data. In addition, when there is a significant mean difference, Approach 1 can also ensure that the type I error 
is commensurate with the pre-specified level and that there will be no loss in power.  In contrast, Approaches 2 and 3 set the adjusted critical value based on the null 
hypothesis of the equivalence test, which results in a lower type I error than the pre-specified level when significant mean differences exist. Moreover, there will also 
be a power loss when the alternative holds.
A brif description, usage, arguments, and returned values of each method is provided in the R documentation. For detailed explaination, please refer to 
the following URL: toBeAdded
