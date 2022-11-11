# Assignment 6

## Question 1

We can use a Dirac function with a convolution to shift a signal around. More precisely, we know that

$$ f(t-k)=f(t) \circledast \delta(t-k) $$

Hence, using this property, we can make a function that simply convolves an array with a shifted dirac. In this case, we have the `shift` argument which corresponds to a shift in the array's index:

```python
def convo_shift(x, shift):
    dirac = np.zeros(len(x))

    dirac[int(len(x)/2)+shift] = 1

    return np.convolve(x, dirac, 'same')
```

The function works by creating an array of the same length as the input array, then adds the dirac's peak at a position shifted from the center by the inputted amount. Then these two arrays are convolved with each other and returned.

We can see this in action, by shifting a gaussian, by 200 in the following:

![param2power](figs/a6q1_convo.jpg)

## Question 2

### Part a)

We can apply the given theorem to create a convolution using fourier space. This is done using the following code:

```python
def cor_func(a, b, shift = False):
    if not shift:
        return np.fft.ifft(np.fft.fft(a) * np.conjugate(np.fft.fft(b)))
    if shift:
        return np.fft.fftshift(np.fft.ifft(np.fft.fft(a) * np.conjugate(np.fft.fft(b))))
```

This code simply takes the discrete FFTs and IFFTs of the arrays `a` and `b` and outputs the resulting convolution. We added an extra argument to allow to use the `fftshift` method so as to recenter the peak. We can see in the following figure that we obtain a Gaussian in our Fourier space:

![a6q2_a_gauss_cor](figs/a6q2_a_gauss_cor.jpg)

Notice that we also plotted the reshifted Gaussian.

### Part b)

When shifting the Gaussian, the fourier space result ends up looking the same (up to a phase).

![a6q2_b_gauss_cor](figs/a6q2_b_gauss_cor.jpg)

We can see the shifted Gaussian and the unchanged Fourier space function (up to a phase). This is because, in fourier space, the math basically "wraps" the array around.

## Question 3

## Question 4

### Part a)

We know that for geometric series we have:

$$ 
\sum_{k=0}^{n-1} a x^k = \begin{cases}
a \left(\frac{1 - x^n}{1-x}\right)n
,\\
x(n-1)\\
x(n-1)
\end{cases}
$$

In our case we have the sum:

$$\sum_{x=0}^{N-1}\exp\left(-\frac{2\pi ik x}{N}\right)$$

Consequently, reading off the terms in the equation we can substitute in the result from the geometric series:

$$\sum_{x=0}^{N-1}\exp\left(-\frac{2\pi ik x}{N}\right) = \sum_{k=0}^{N-1}  \exp\left(-\frac{2\pi i k}{N}\right)^x = \left(\frac{1 - \exp\left(-\frac{2\pi i k}{N}\right)^{N}}{1-\exp\left(-\frac{2\pi i k}{N}\right)}\right)$$

$$= \left(\frac{1 - \exp\left(-2\pi i k\right)}{1-\exp\left(-\frac{2\pi i k}{N}\right)}\right)$$

### Part b)

## Question 5

### Part a)

To model the noise we smooth out the data using a convolution with a fat gaussian. This basically averages out neighboring points. We can see this effect in the following picture:

![a6q5_comp_smooth_ps](figs/a5q6_comp_smooth_ps.jpg)

The orange curve is the smoothed data in Fourier space and the blue curve is the raw data in Fourier space. The raw time-domain data was put in Fourier space using a cosine windowing to prevent artefact in the frequency domain. After playing with this and being able to find a peak for the GW event, I didn't find useful to use an alternate windowing scheme. Maybe there was some noise reduction in some of the options I tried, but I did not notice anything significant, hence I decided to only use the cosine window.

Using this smoothed out curve, we can now whiten the data. Intuitively, we want to flatten out to frequency response to make it closer to white noise. To do this, we can simply divide the raw Fourier spectrum by the smoothed out version. This does indeed flatten out our spectrum as shown in the following plot:

![a6q5_comp_smooth_ps](figs/a6q5_whitened_ps.jpg)

We indeed see the power spectrum being flattened, except possibly for the very begining. This is not really a problem since that range (on a log scale) does not represent many points in comparison to the rest of the points and also it is not in the range where our sought-after signal exists.

### Part b) 

