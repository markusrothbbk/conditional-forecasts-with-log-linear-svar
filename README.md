# conditional-forecasts-with-log-linear-svar
Forecasting with Log-Linear (S)VAR Models: Incorporating Annual Growth Conditions

## Forecast condition types
Conditions are specified through the conditions argument and can have various types defined by the type field. Below is a summary of the available condition types:

---

### Type 1: Simple conditions on variable values
Imposes that variable $j$ takes a fixed value in periods $T + h_0$ to $T + h_1$:

$y_{T+h}(j) = \text{value}, \quad \text{for} \quad h = h_0, \dots, h_1$

---

### Type 2: Simple conditions on shock values
Imposes that the structural shock $\varepsilon_{T+h}(j)$ takes a fixed value:

$\varepsilon_{T + h}(j) = \text{value}, \quad \text{for} \quad h = h_0, \dots, h_1$

---

### Type 3: Conditions on multiple-period averages (variables in levels)
Imposes that the average of variable $j$ over periods $h_0$ to $h_1$ equals a value:

$\frac{1}{n} \sum_{h = h_0}^{h_1} y_{T+h}(j) = \text{value}, \quad n = h_1 - h_0 + 1$

---

### Type 4: Conditions on average changes (variables in levels)
Imposes that the average change of variable $j$ between two periods equals a value:

$\frac{1}{n} \left( \sum_{h = h_0}^{h_1} y_{T+h}(j) - \sum_{h = h_0 - n}^{h_1 - n} y_{T+h}(j) \right) = \text{value}, \quad n = h_1 - h_0 + 1$

---

### Type 5: Conditions on average changes (first-difference variables)
Similar to Type 4, but for variables specified in first differences:

$\frac{1}{n} \left( \sum_{h = h_0}^{h_1} y_{T+h}(j) - \sum_{h = h_0 - n}^{h_1 - n} y_{T+h}(j) \right) = \text{value}$

---

### Type 6: Conditions on average growth rates (variables in log levels)
Imposes that the average growth rate over periods equals a value:

$\left( \frac{\sum_{h = h_0}^{h_1} z_{T+h}(j)}{\sum_{h = h_0 - n}^{h_1 - n} z_{T+h}(j)} \right) - 1 = \text{value}$

where $y_t(j) = \ln(z_t(j))$, and $n = h_1 - h_0 + 1$.

---

### Type 7: Conditions on average growth rates (first-difference log variables)
Similar to Type 6, but for variables in first differences of logs:

$\left( \frac{\sum_{h = h_0}^{h_1} z_{T+h}(j)}{\sum_{h = h_0 - n}^{h_1 - n} z_{T+h}(j)} \right) - 1 = \text{value}$

where $y_t(j) = \ln(z_t(j)) - \ln(z_{t-1}(j))$.

---

### Type 8: Conditions on multiple-period averages (variables in log levels)
Imposes that the average level of $z_t(j)$ equals a value:

$\frac{1}{n} \sum_{h = h_0}^{h_1} z_{T+h}(j) = \text{value}$

where $y_t(j) = \ln(z_t(j)) ), and ( n = h_1 - h_0 + 1$.

### Type 9: Conditions on multiple-period changes (variables in levels)
Imposes that the change in variable $j$ between periods equals a value:

$y_{T+h_1}(j) - y_{T+h_0}(j) = \text{value}$

---
