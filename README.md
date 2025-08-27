# Forecasting with Log-Linear (S)VAR Models: Incorporating Annual Growth Conditions      
      
This repository provides MATLAB code for conditional forecasting using (S)VAR models. It allows for several special cases such as annual growth rate conditions. Some of these special cases are described in Mokinski and Roth (2025).      
      
Authors: Frieder Mokinski and Markus Roth    
      
## Files      
    
| File                               | Description                                                                                      |      
|------------------------------------|--------------------------------------------------------------------------------------------------|      
| `conditionalForecast.m`            | Main function: Computes conditional forecasts. See `conditionsTable.xlsx` for details on specifying conditions. |      
| `conditionsTable.xlsx`             | Documentation: Explains specification of conditions for `conditionalForecast.m` function. |      
| `checkVarianceMatrix.m`            | Function: Checks if a given matrix is a valid variance-covariance matrix (square, symmetric, positive semi-definite). Used internally for validation. |      
| `fvar.m`                           | Function: Estimates a frequentist VAR(p) model with or without a constant. |      
| `lagn.m`                           | Function: Creates lagged (and/or leaded) versions of a data matrix or table for given lags. |      
| `simvar.m`                         | Function: Simulates time series data from a VAR model with given parameters and initial values. |      
| `vec.m`                            | Function: Vectorizes a matrix (row-wise or column-wise). |      
| `isPositiveInteger.m`              | Function: Checks if a value is a positive integer (helper for function argument validation). |      
| `conditionalForecast_examples.m`   | MATLAB script: Simulates macroeconomic time series, estimates VAR models, and demonstrates conditional forecasting for a model specified in levels and one in first differences. Code uses the examples of `conditionsTable.xlsx`.|    
| `conditionalForecast_applicationNote.m` | MATLAB script: Reproduces Table 1 of the Mokinski and Roth (2025) that demonstrates how the linearization of the annual growth rate condition leads to an approximation error that can be reduced through the algorithm of Section 3 of the paper. |  
    
## References      
      
- Camba-Mendez, G. (2012): *Conditional Forecasts on SVAR models using the Kalman filter*, Economics Letters, Elsevier, volume 115, issue 3, pages 376-378.      
- Mokinski, F. & Roth, M. (2025): *Forecasting with Log-Linear (S)VAR Models: Incorporating Annual Growth Rate Conditions*, Unpublished manuscript.      

## License  
  
This project is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html).  
  
You are free to use, modify, and distribute this software under the terms of the GPL.    
See the [LICENSE](LICENSE) file for more details.  
  
## Disclaimer  
  
This software is provided "as is", without warranty of any kind, express or implied.    
See the GNU General Public License for more details.  
