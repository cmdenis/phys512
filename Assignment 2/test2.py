import numpy as np

nb_sub = 11
start_val = 0
end_val = 2

first_array = np.linspace(start_val, end_val , nb_sub)

mid = (max(first_array)+min(first_array))/2


new_int = mid - start_val

second_array = np.linspace(start_val + new_int/(nb_sub-1), mid - new_int/(nb_sub-1), int((nb_sub - 1)/2))

print(np.arange(int((nb_sub+1)/2)))

second_array = np.insert(
    second_array, 
    np.arange(int((nb_sub+1)/2)).astype(int), 
    first_array[:int((nb_sub - 1)/2)+1]
    )


#farr = np.array([1, 2, 3, 4])
#sarr = np.array([5, 6, 7, 8, 9])
#test_array = np.insert(farr, [0, 1, 2, 3,4], sarr )
#print(test_array)


print(first_array)
print(second_array)

objective = np.linspace(start_val, end_val/2, nb_sub)
print(objective)

# if (num % 2) == 0