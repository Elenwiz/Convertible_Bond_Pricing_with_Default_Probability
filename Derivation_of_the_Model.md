## Derivation of the Model

Our model is based on a one-state variable, $S_t$, which describes the convertible bond's underlying stock price. Traditional methods assume that $S_t$ follows the following stochastic process:
$$
\mathrm{d}S_t = r_t S_t \mathrm{d}t + \sigma S_t \mathrm{d}W_t
$$


where $\mathrm{d}S_t$ and $\mathrm{d}W_t$ are the increments for the infinitesimal time period of the stock price and Wiener process. Different from that, we want to take the default risk of the company and the dividends of the stock into consideration. According to Clark and Weinstein (1983), the stock price drops on average by about 30% upon default. For this reason, we assume that the return of the underlying stock follows a process that is a combination of a diffusion process and a jump process, which means  $S_t$ follows:
$$
\mathrm{d}S_t = (r_t + \lambda \theta- v_t) S_t \mathrm{d}t + \sigma S_t \mathrm{d}W_t - \theta S_t \mathrm{d}q_t \label{1} \tag{1}
$$
where $\mathrm{d}q_t$ is the increment for the infinitesimal time period of the homogeneous Poisson process with intensity $\lambda$. In addition, we assume that there is no correlation between the Wiener process and the Poisson process. In Equation $\eqref{1}$,  $r_t$ represents the risk-free interest rate,  $v_t$ the dividend ratio of the stock, and $\theta$ the percentage that the stock's price falls immediately after default ($r>0$, $\sigma_t >0$, $\lambda \ge 0$, $0\le\theta\le1$). For the last one, we can prove it by using the following statement.

Let the arrival time for the Poisson process be $\tau$, $\tau^+$ be the forward time the instant after that arrival, and $\tau ^-$ be the backward time the instant before that arrival.  Since Ito's lemma for a process which is the sum of a drift-diffusion process and a jump process is just the sum of Ito's lemma for the individual parts, so we can get
$$
\ln (S_{\tau^+}) - \ln (S_{\tau^-}) = (r_t + \lambda \theta- v_t-\frac{1}{2}\sigma_t^2)\mathrm{d}t + \sigma_t \mathrm{d}W_t + \ln({1-\theta})\mathrm{d}q_t
$$


This implies the following equality,
$$
\ln (S_{\tau^+}) - \ln (S_{\tau^-}) = \ln (S_{\tau^+}^s) - \ln (S_{\tau^-}) + \ln({1-\theta})
$$


where $S_{\tau^+}$ refers to the value of the process in just one arrival ($\mathrm{d}q_t=1$) and $S_{\tau^-}^s$ describes the value of the process without an arrival ($\mathrm{d}q_t=0$).  So we can get
$$
\begin{align*}
& \ln (S_{\tau^+}) -  \ln (S_{\tau^+}^s)= \ln({1-\theta}) \\
& S_{\tau^+} = S_{\tau^+}^s (1-\theta)

\end{align*}
$$
Therefore, we know that through the moment of exactly one arrival of the Poisson process, the stock value drops by exactly $\theta$ percent compared to the case without an arrival.

## Binomial Tree of the Stock

Now, we need to calculate several parameters to build a binomial random walk of the underlying stock. To do this,  let us consider Equation (2) in discrete time with step $\delta t$. In this discretization, we obtain
$$
\begin{align*}
\delta S_t & = (r_t + \lambda \theta- v_t) S_t \delta t + \sigma S_t \delta W_t - \theta S_t \delta q_t \\
\frac{\delta S_t}{S_t} & = (r_t + \lambda \theta- v_t) \delta t + \sigma \delta W_t - \theta \delta q_t

\end{align*} \tag{2}
$$
Then we can have the approximations for the stock returns, with the expectation and the variance below:
$$
\begin{align*}
& \mathrm{E}(\frac{\delta S_t}{S_t}) = (r_t + \lambda \theta- v_t) \delta t  - \lambda\theta \delta q_t = (r_t-v_t)\delta t \\
& \mathrm{Var}(\frac{\delta S_t}{S_t}) = \sigma^2 \delta t + \lambda \theta^2 \delta t = (\sigma^2 + \lambda \theta^2)\delta t
\end{align*} \tag{3}
$$
Since $S_{t+\delta t} = S_{t} + \delta S_t$, we can further get:
$$
\begin{align*}
& \mathrm{E}(\frac{S_{t+\delta t}}{S_t}) = \mathrm{E}(\frac{\delta S_t}{S_t}) + 1 = 1 + (r_t - v_t)\delta t \\
& \mathrm{Var}(\frac{S_{t+\delta t}}{S_t}) = \mathrm{Var}(\frac{\delta S_t}{S_t}) = (\sigma^2 + \lambda \theta^2)\delta t
\end{align*} \tag{4}
$$
Given these two conditions, now we can represent random walks of a stock by using a binomial tree. Before we start, it should be noted that for the given time $\delta t$, the Poisson increment $\delta q_t$ may count more than one arrival. In this article, we model the event of default in the period $\delta t$ as the condition $\{  \delta q_t > 0 \} $. So, throughout every time step of the binomial tree, the probability of default $p_0$ is equal to $1 - e^{-\lambda\delta t}$. There is another point that should be mentioned in our model, which is once the stock suffers a $\theta$ percent decrease in its price after default, the binomial tree never moves further. Hence, there is an imaginary free node that states the value of the stock price after default in the structure of the binomial tree, as shown below.

![Figure 1](Tree_Structure.png)

All the possible values of $\frac{S_{t+\delta t}}{S_t}$ along the binomial tree with time $\delta t$ are $u$, $d$, and $(1-\theta)$ with probabilities $p_u$, $p_d$, and $p_0$ ($0\le p_u,p_d,p_0 \le 1$). We can get the following equations:
$$
\begin{cases}
p_0 + p_u + p_d = 1\\
\mathrm{E}(\frac{S_{t+\delta t}}{S_t}) = p_u u + p_d d + p_0 (1-\theta) \\
\mathrm{Var}(\frac{S_{t+\delta t}}{S_t})  = \mathrm{E}[(\frac{S_{t+\delta t}}{S_t})^2] - \mathrm{E}^2(\frac{S_{t+\delta t}}{S_t}) = p_u u^2 + p_d d^2 + p_0 (1-\theta)^2 - \mathrm{E}^2(\frac{S_{t+\delta t}}{S_t})
\end{cases}
$$


where $p_0 = 1-e^{-\lambda \delta t}$, $p_d = e^{-\lambda \delta t} - p_u$. Regarding $\mathrm{E}(\frac{S_{t+\delta t}}{S_t})$, we obtain:
$$
\begin{align*}
& p_u u + (e^{-\lambda \delta t} - p_u) d + (1-e^{-\lambda \theta})(1-\theta) = 1 + (r_t- v_t)\delta t \approx e^{(r_t - v_t)\delta t} \\
& p_u (u - d) = e^{(r_t - v_t)\delta t} - e^{-\lambda \delta t} d - (1-\theta) (1- e^ {-\lambda \delta t}) \\
& p_u = \frac{e^{(r_t - v_t)\delta t} - e^{-\lambda \delta t} d - (1-\theta) (1- e^ {-\lambda \delta t})}{u - d} \tag{5}\label{5}
\end{align*}
$$
Regarding  $\mathrm{Var}(\frac{S_{t+\delta t}}{S_t})$, we obtain:
$$
p_u (u^2 -d^2) + e^{-\lambda \delta t} d^2 + (1-e^{-\lambda \delta t})(1-\theta)^2 - e^{2(r_t - v_t)\delta t} = (\sigma_t^2 + \lambda \theta^2)\delta t \tag{6} \label{6}
$$
To satisfy the assumptions of the binomial tree, we look for $u$ and $d$ satisfying $ud = 1$. Moreover, we write $u$ in the form $u = e^{\sqrt{A\delta t}}$, which can help us to proceed with Equation $\eqref{6}$, ignoring all terms containing powers of $\delta t$ greater than 2. After substituting Equation $\eqref{5}$, we get:
$$
\begin{align*}
(\sigma_t^2 + \lambda \theta^2)\delta t = & [1 + (r_t - v_t)\delta t -(1 - \lambda\delta t)d - (1 - \theta)\lambda \delta t]( u + d) \\
 & + (1 - \lambda \delta t)d^2 + (1-\theta)^2 \lambda \delta t - [1+2(r_t - v_t)\delta  t] \\
 = & [1 + (r_t - v_t)\delta t - (1 - \theta)\lambda \delta t]( u + d) - (1 - \lambda\delta t)ud \\
 & + (1-\theta)^2 \lambda \delta t - [1+2(r_t - v_t)\delta  t]
\end{align*}
$$
Then, by following the approximation
$$
u \approx 1 + \sqrt{A\delta t} + \frac{A\delta t}{2} + \frac{\sqrt{A\delta t}^{1.5}}{6} + \frac{(A\delta t)^2}{24} \\
d \approx 1 - \sqrt{A\delta t} + \frac{A\delta t}{2} - \frac{\sqrt{A\delta t}^{1.5}}{6} + \frac{(A\delta t)^2}{24} \\
u + d \approx 2 + A\delta t + \frac{(A\delta t)^2}{12}
$$
the equation for the variance can be written as
$$
\begin{align*}
(\sigma_t^2 + \lambda \theta^2)\delta t = & [1 + (r_t - v_t)\delta t - (1 - \theta)\lambda \delta t](2 + A \delta t) - (1 - \lambda\delta t) + (1 + \theta^2 - 2 \theta ) \lambda \delta t - 1 - 2(r_t - v_t)\delta  t \\
= & [1 + (r_t - v_t)\delta t - (1 - \theta)\lambda \delta t](2 + A \delta t) - 1 + \lambda \delta t + \lambda \delta t + \theta^2 \lambda \delta t - 2\theta \lambda \delta t - 1 - 2(r_t - v_t) \delta t \\
= & 2 + 2(r_t - v_t)\delta t - 2\lambda \delta t + 2\theta \lambda \delta t + A\delta t - 2 + 2\lambda\delta t - 2\theta \lambda \delta t + \theta ^2 \lambda \delta t - 2(r_t - v_t)\delta t \\
= & (A+\lambda \theta^2)\delta t
\end{align*}
$$
which means,
$$
A = \sigma_t^2
$$
Now, we can get all the expressions of parameters of the binomial tree for a defaultable stock that follows process $\eqref{1}$,  as shown in ==Table 1==.

| Parameter                   | Expression                                                   |
| --------------------------- | ------------------------------------------------------------ |
| Length of tree step         | $\delta t$                                                   |
| Multiplier for moving up    | $u = e^{\sigma_t\sqrt{\delta t}}$                            |
| Multiplier for moving down  | $d = e^{-\sigma_t \sqrt{\delta t}}$                          |
| Probability for moving up   | $p_u = \frac{e^{(r_t - v_t)\delta t} - e^{-\lambda \delta t} d - (1-\theta) (1- e^ {-\lambda \delta t})}{u - d}$ |
| Probability for moving down | $p_d = - \frac{e^{(r_t - v_t)\delta t} - e^{-\lambda \delta t} u - (1-\theta) (1- e^ {-\lambda \delta t})}{u - d} $ |
| Probability for default     | $p_0 = 1 - e ^ {-\lambda \delta t}$                          |

## Convertible Bond Pricing

### Considering Hold to Maturity

Now we can turn to the valuation of convertible bonds within the binomial tree model for the price of the underlying stock derived above. We begin with a portfolio, which consists of ==one convertible bond (first we only consider the situation where holders choose to hold all the time)== and a short proportion of the underlying stock with quantity $\Delta$. At time $t$, this portfolio has a value of
$$
\Pi_t = H_t - \Delta S_t
$$
where $H_t$ denotes the convertible bond's value, and $S_t$ is the stock's value. In time $\delta t$, the convertible bond takes one of three values, $H_t^u$, $H^d_t$, or $H^0_t = \max(RF, n_t(1-\theta)S_t)$, depending on whether the underlying stock rises to a value $S_t u$, falls to a value $S_t d$, or drops to $(1-\theta)S_t$. Here, $R$ refers to the recovery rate, $F$ the face value of the convertible bond, and $n_t$ the conversion ratio valid at time $t+\delta t$ (which equals to $F/S_c$, $S_c$ stands for the conversion price). Note that according to Ayache et al. (2003), we incorporate the convertible bond holder's option to convert the bond after the bankruptcy. So, the possible values of the portfolio at time $t+\delta t$ are:
$$
\Pi_{t+\delta t} =
\begin{cases}
H^u_t-\Delta S_t u \\
H^d_t-\Delta S_t d \\
H^0_t-\Delta S_t (1- \theta)
\end{cases}
$$
For a risk-neutral portfolio, there is no risk of diffusion in it. So, we get the hedge quantity $\Delta$ from the following equation:
$$
H^u_t-\Delta S_t u = H^d_t-\Delta S_t d
$$
So we get 
$$
\Delta = \frac{H_t^u - H_t^d}{S_t(u-d)}
$$
Combined with the probabilities we get in the prior section, the expectation of $\Pi_{t+\delta t}$ can be written as
$$
\begin{align*}
\mathrm{E}(\Pi_{t+\delta t}) & = [H_t^u - \Delta S_t u] (p_u + p_d) + [H_t^0 - \Delta S_t (1-\theta)]p_0 \\
& = \frac{H_t^d u - H_t^u d}{u-d}e^{-\lambda \delta t} + [H_t^0 - \frac{H_t^u - H_t^d}{(u-d)}(1-\theta)](1-e^{-\lambda \delta t}) \\
& = -\frac{e^{-\lambda \delta t}d + (1-\theta)(1-e^{-\lambda \delta t})}{u-d}H_t^u + \frac{e^{-\lambda \delta t}u + (1-\theta)(1-e^{-\lambda \delta t})}{u-d} H_t^d + (1-e^{-\lambda \delta t})H_t^0
\end{align*}
$$
Further, assuming that the default risk is diversifiable (the portfolio return is based on a risk-free interest rate no matter whether the default occurs or not), which implies
$$
\mathrm{E}(\Pi_{t+\delta t}) = \Pi_{t}e^{(r_t-v_t)\delta t}
$$
where $r_t$ is the continuously compounded interest rate. After substituting with the equation above and proceeding with grouping, we obtain
$$
\begin{align*}
e^{(r_t - v_t)\delta t} H_t = & \frac{e^{(r_t - v_t)\delta t} - e^{-\lambda \delta t} d - (1-\theta) (1- e^ {-\lambda \delta t})}{u - d} H_t^u \\
& - \frac{e^{(r_t - v_t)\delta t} - e^{-\lambda \delta t} u - (1-\theta) (1- e^ {-\lambda \delta t})}{u - d}H_t^d \\
& + (1-e^{-\lambda \delta t})H_t^0 \\
\end{align*}
$$
Rewrite the above equation in a more precise way, which is
$$
H_t = e^{-(r_t-v_t)\delta t} (p_u H_t^u + p_d H_t^d + p_o H_t^0)
$$
indicating that convertible bond valuation within our binomial tree model is just as same as that in the classical framework. Moreover, it can be easily proved that when we take coupon cash flow into account. If we denote $t_i^c$ the moment at which the convertible bond paid out a coupon amount $c_i$ and assume that this happens within the time step $\delta t$, then we can get the following equation of $H_t$:
$$
H_t = e^{-(r_t-v_t)\delta t} (p_u H_t^u + p_d H_t^d + p_o H_t^0) + e^{-(r_t-v_t)(t_i^c-t)-\lambda\delta t}c_i
$$
showing that the risk-neutral present value of the coupon cash flow at the moment $t$ has to be adjusted by the probability for non-default ($e^{-\lambda \delta t}$) during the period $\delta t$.

### Including Early Conversion and Conversion Probability

In this section, we further consider those complex features of convertible bonds, including call provision,  put provision, and early conversion (to simplify, we do not take those trigger conditions into account in the modeling part). Our valuation of convertible bonds now follows:
$$
V_t = \max(n_tS_t, P_t, \min(H_t + c_t, C_t)) \label{7} \tag{7}
$$
In Equation $\eqref{7}$,  $n_t$ denotes the conversion ratio valid at time $t+\delta t$, $P_t$ the value of the convertible bond when the bondholders use the put provision, $c_t$ the coupon amount paid out within the time step, $H_t$ the value of the convertible bond when bondholders continue holding, and $C_t$  the value of the bond when the issuer uses the call provision. At time point $t$, $H_t$ follows:
$$
H_t = p_u V_t^u e^{-(y_t^u-v_t)\delta t} + p_d V_t^d e^{-(y_t^d-v_t)\delta t} + p_0 V_t^0 e^{-(y_t^0-v_t)\delta t}\label{8}\tag{8}
$$
in which $V_t^u$ refers to the value of the convertible bond when the price of the underlying stock rises, $V_t^d$ the case when the stock price drops, and $V_t^d$ the case when default happens, respectively. $y_t^u$ refers to the discount rate when the stock price rises, $y_t^d$ the case when the price drops, and $y_t^0$ the case when default happens. We do this modification to consider the probability for holders to choose to convert the bond into the stock at every time point $t$. In every node, the discount rate equals:
$$
y_t = p_t \cdot r_t + (1-p_t)\cdot (r_t+Spread_t) \tag{9}\label{x}
$$
where $p_t$ stands for the conversion probability and $Spread_t$ is the credit spread of the bond. Specifically, the calculation of $p_t$ takes one of three conditions: (1) if the conversion occurs at $t$, $p_t = 1$; (2) if the redemption or maturity occurs at $t$, $p_t = 0$; (3) if neither condition (1) nor (2), $p_t = p_u p_{t+\delta t}^u + p_d p_{t+\delta t}^d + p_0 p_{t+\delta t}^0$, where $p_{t+\delta t}^u$ stands for the conversion probability at $t+\delta t$ when the stock rises and $p_{t+\delta t}^d$ the case when the stock falls. By doing so, we are able to consider the convertible bond as a pure bond when the option part is deeply out of money by using the interest rate including credit spread, as an equity asset when the option part is in the money respectively by using risk-free interest rate as a discount factor, and give a smooth transition between these two rates with the help of conversion probability.

