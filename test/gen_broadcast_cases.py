import numpy as np
import random

class array_generator:
    def __init__(self, file, n_dim_lim, dim_lim, one_prob, resize_prob):
        self.file = file
        self.seed = 0
        self.n_dim_lim      = n_dim_lim
        self.dim_lim        = dim_lim
        self.one_prob       = one_prob
        self.resize_prob    = resize_prob

    def rand(self, limit):
        return random.randrange(1,limit)

    def rand_ones(self, array, prob):
        # Set randomly values within array to one
        choixe = np.random.rand(*array.shape)
        array[choixe<prob] = 1

    def with_shape_limit(self):
        # n_shape_lim:  limit number of dimensions
        # dim_lim:    limit of each dimension
        shape = np.array([self.rand(self.dim_lim) for _ in range(random.randrange(1, self.n_dim_lim))])
        return np.random.randint(-1000, 1000, np.prod(shape)).reshape(shape)

    def compatible_with(self, array):
        #array: array the generated array must be compatible with.
        # Will get resized (with prob. resize_prob), and have random dimensions set to 1 (with prob. one_prob)
        start_pos = 0
        if self.rand(100) < int((100*self.resize_prob)) and len(array.shape) > 1:
            r = random.randrange(1, len(array.shape))
            start_pos = r
        shape = np.array(array.shape)[start_pos:]
        self.rand_ones(shape, self.one_prob)
        return np.random.randint(-1000, 1000, np.prod(shape)).reshape(shape)

    
    def write_to(self, array):
        if isinstance(array, list):
            for a in array:
                s = ""
                for i in a.shape:
                    s = s + str(i) + " "
                s = s + " content "
                for i in np.nditer(a):
                    s = s + str(i) + " "
                self.file.write(s + "\n")
            
            f.write("\n")
        else:
            print("Invalid input")

def generate_test_addition(n_samples, a_gen):
    with open("test_cases.txt", "w") as f:
        for i in range(n_samples):
            lhs = a_gen.with_shape_limit()
            rhs = a_gen.compatible_with(lhs)
            a_gen.write_to([lhs, rhs, lhs + rhs])

if __name__ == "__main__":
    print("Generating test cases...")
    with open("test_cases.txt", "w+") as f:
        # Generate small arrays
        a_gen = array_generator(f, 5, 10, 0.2, 0.8)
        generate_test_addition(500, a_gen)
