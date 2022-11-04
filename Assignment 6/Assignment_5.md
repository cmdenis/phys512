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


