import numpy as np

def Half_fill(w_base, w_add, N):
    A = np.ones((N, N))*w_base
    for i in range(N):
        
        if (i % 2 == 0):
            j_start = 1
        else:
            j_start = 0
            
        for j in range(j_start, N, 2):
            A[i, j] = w_base + w_add
    
    return A

def Random_swap(W, N):
    A = W.copy()

    N_i = A.shape[0]
    N_j = A.shape[1]

    i_a = np.random.randint(low=0, high=N_i, size=N)
    j_a = np.random.randint(low=0, high=N_j, size=N)

    i_b = np.random.randint(low=0, high=N_i, size=N)
    j_b = np.random.randint(low=0, high=N_j, size=N)

    for k in range(N):
        tmp = A[i_a[k], j_a[k]]
        A[i_a[k], j_a[k]] = A[i_b[k], j_b[k]]
        A[i_b[k], j_b[k]] = tmp
    
    return A


W = Half_fill(2.0, 0.5, 5)
print(W)
A = Random_swap(W, 1000)
print('')
print(A)