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
