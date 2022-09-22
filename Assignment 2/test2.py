import numpy as np

# Compute stuff initially
nb_sub = 11
start_val = 0
end_val = 2

first_array = np.linspace(start_val, end_val , nb_sub)

mid = (max(first_array)+min(first_array))/2



# Move on to second recursion

new_start_val = start_val
new_end_val = mid
new_int = new_end_val - new_start_val # New interval size


# Make "half filled array with holes to be filled in between elements"
second_array = np.linspace(
    new_start_val + new_int/(nb_sub-1), 
    new_end_val - new_int/(nb_sub-1), 
    int((nb_sub - 1)/2)
    )

# Fill in the blanks in the final array
second_array = np.insert(
    second_array, 
    np.arange(int((nb_sub+1)/2)).astype(int), 
    first_array[:int((nb_sub - 1)/2)+1]
    )


#farr = np.array([1, 2, 3, 4])
#sarr = np.array([5, 6, 7, 8, 9])
#test_array = np.insert(farr, [0, 1, 2, 3,4], sarr )
#print(test_array)

print


print(first_array)
print(second_array)

objective = np.linspace(start_val, end_val/2, nb_sub)
print(objective)

# if (num % 2) == 0