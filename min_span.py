# given a genetic map and a length cutoff in cM, calculate the minimum span (in bp)
import numpy as np
from tqdm import tqdm

f = open("example_data/genetic_map_GRCH38_chr20.txt")
data = f.readlines()
for i in range(len(data)):
    data[i] = data[i].strip().split()
    data[i][1] = int(data[i][1])
    data[i][3] = float(data[i][3])
print(data[0])
max_bp = data[-1][1]
sites = [i for i in range(max_bp)]
xp = [i[1] for i in data]
yp = [i[3] for i in data]
interp = np.interp(sites, xp, yp)
min_span = max_bp
def find_min_int_diff_for_float_diff(int_list, float_list, target_float_diff):
    n = len(int_list)
    min_int_diff = float('inf')
    
    # Initialize two pointers
    i, j = 0, 0
    
    while i < n and j < n:
        float_diff = abs(float_list[i] - float_list[j])
        
        if float_diff == target_float_diff:
            int_diff = abs(int_list[i] - int_list[j])
            if int_diff < min_int_diff:
                min_int_diff = int_diff
            # Move the pointer that is on the smaller index to find other potential pairs
            if i < j:
                i += 1
            else:
                j += 1
        elif float_diff < target_float_diff:
            j += 1
        else:
            i += 1
    
    return min_int_diff if min_int_diff != float('inf') else None
    

find_min_int_diff_for_float_diff(sites, interp, 2)