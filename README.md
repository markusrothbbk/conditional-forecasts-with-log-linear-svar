# conditional-forecasts-with-log-linear-svar
Forecasting with Log-Linear (S)VAR Models: Incorporating Annual Growth Conditions

## Forecast condition types
Conditions are specified through the conditions argument and can have various types defined by the type field. Below is a summary of the available condition types:

---

### Type 1: Simple condition on variables' single-period values
Imposes that variable $j$ takes a fixed value in periods $T + h_0$ to $T + h_1$:

$y_{T+h}(j) = \text{value}, \quad \text{for} \quad h = h_0, \dots, h_1$

---

### Type 2: Simple condition on shocks' single-period values
Imposes that the structural shock $\varepsilon_{T+h}(j)$ takes a fixed value:

$\varepsilon_{T + h}(j) = \text{value}, \quad \text{for} \quad h = h_0, \dots, h_1$

---

### Type 3: Condition on average (e.g. annual) level if variables are in levels
Also: condition on multiple-period growth rate if variables are in log first differences

Imposes that the average of variable $j$ over periods $h_0$ to $h_1$ equals a value:

$\frac{1}{n} \sum_{h = h_0}^{h_1} y_{T+h}(j) = \text{value}, \quad n = h_1 - h_0 + 1$

---

### Type 4: Condition on change in the average (e.g. annual) level if variables are in levels
Imposes that the average change of variable $j$ between two periods equals a value:

$\frac{1}{n} \left( \sum_{h = h_0}^{h_1} y_{T+h}(j) - \sum_{h = h_0 - n}^{h_1 - n} y_{T+h}(j) \right) = \text{value}, \quad n = h_1 - h_0 + 1$

---

### Type 5: Condition on change in the average (e.g. annual) level if variables are in first differences
Similar to Type 4, but for variables specified in first differences:

$\frac{1}{n} \left( \sum_{h = h_0}^{h_1} y_{T+h}(j) - \sum_{h = h_0 - n}^{h_1 - n} y_{T+h}(j) \right) = \text{value}$

---

### Type 6: Condition on growth rate of the average (e.g. annual) level if variables are in log levels
Imposes that the average growth rate over periods equals a value:

$\left( \frac{\sum_{h = h_0}^{h_1} z_{T+h}(j)}{\sum_{h = h_0 - n}^{h_1 - n} z_{T+h}(j)} \right) - 1 = \text{value}$

where $y_t(j) = \ln(z_t(j))$, and $n = h_1 - h_0 + 1$.

---

### Type 7: Condition on growth rate of the average (e.g. annual) level if variables are in log first differences
Similar to Type 6, but for variables in first differences of logs:

$\left( \frac{\sum_{h = h_0}^{h_1} z_{T+h}(j)}{\sum_{h = h_0 - n}^{h_1 - n} z_{T+h}(j)} \right) - 1 = \text{value}$

where $y_t(j) = \ln(z_t(j)) - \ln(z_{t-1}(j))$.

---

### Type 8: Condition on average (e.g. annual)  level if variables are in log levels
Imposes that the average level of $z_t(j)$ equals a value:

$\frac{1}{n} \sum_{h = h_0}^{h_1} z_{T+h}(j) = \text{value}$

where $y_t(j) = \ln(z_t(j)) ), and ( n = h_1 - h_0 + 1$.

### Type 9: Condition on multiple-period change if variables are in levels
Imposes that the change in variable $j$ between periods equals a value:

$y_{T+h_1}(j) - y_{T+h_0}(j) = \text{value}$

---
