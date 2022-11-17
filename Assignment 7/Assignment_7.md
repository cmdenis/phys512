# Assignment 7

By Christian Denis
For PHYS-512

## Question 1

## Question 2

For this question I initially took a recursive approach, but for the following reasons:

1. Speed
2. Difficulty/inelegance of counting the rejection proportion

I decided to go with a vectorizing approach. In a different language than python, the recursive version might have worked better. I also think it's slightly more elegant, but anyways... the code to do this looks like this:

```python
def rand(pdf, bounds):
    # Generating 2 random numbers
    r = np.random.rand(2)

    # Rescaling the uniform distribution
    dist = bounds[1] - bounds[0]
    r[0] = r[0]*dist + bounds[0]

    # Check if they are accepted, else try again recursively
    if r[1] < pdf(r[0], mean, sigma):
        return r[0]
    else:
        return rand(pdf, bounds)
```

In any case, I went with the vectorizing method, which works as follows:

```python
def rand_np(n, bounds, pdf):
    # Generating two random arrays of numbers
    r1 = np.random.rand(n)
    r2 = np.random.rand(n)

    # Rescaling the uniform distribution
    dist = bounds[1] - bounds[0] 
    r1 = r1*dist + bounds[0]
    
    # Create an "accept" array
    accept = r2 < pdf(r1)

    # Return array containing only accepted values
    return r1[accept]
```

The way this works is that the function `rand_np` takes as input

1. `n`: the number of "trials" for the random numbers. To find thee rejection rate we can easily compare `n` with the length of the outputted array.
2. `bounds`: the delimitations of the PDF.
3. `pdf` the desired PDF function.

Inside the function, we create two arrays of random numbers with length `n`. One of which we rescale its value so they are uniformly distributed within the input bounds. We then compare the values of the uniformly distributed between 0 and 1 second array to the values of the rescaled arrays when they are put into the PDF. This produces a Boolean array which we can use to select the values of the rescaled arrays which will "pass through" and will be distributed according to our PDF.

In the previous case, we used a 


## Question 3

