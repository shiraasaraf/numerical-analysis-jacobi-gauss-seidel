# Numerical Analysis - Solving Linear Systems using Jacobi and Gauss-Seidel Methods

# Set global parameters for maximum iterations and tolerance
MAX_ITERATIONS = 1000
TOLERANCE = 0.00001

# Define the system of equations
matrixA = [
    [1, 3, -2],
    [4, 1, -1],
    [2, -1, 1]
]
# Define the right-hand side vector
vectorB = [5, 6, 7]

# Define the initial guess for the solution
initial_guess = [0 for _ in range(len(vectorB))]

# Check if a matrix is diagonally dominant
def is_diagonally_dominant(matrix):
    for i in range(len(matrix)):
        row_sum = sum(abs(matrix[i][j]) for j in range(len(matrix)) if i != j)
        if abs(matrix[i][i]) < row_sum:
            return False
    return True

# Try to rearrange the matrix to make it diagonally dominant
def rearrange_to_diagonal_dominance(matrix, vector):
    n = len(matrix)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(matrix[j][i]) > abs(matrix[i][i]):
                matrix[i], matrix[j] = matrix[j], matrix[i]
                vector[i], vector[j] = vector[j], vector[i]
    return is_diagonally_dominant(matrix)

# Calculate the infinity norm between two vectors
def vector_norm_inf(v1, v2):
    return max(abs(a - b) for a, b in zip(v1, v2))

# Solve the system using the Jacobi method
def jacobi_method(matrix, vector, initial_guess, tol=TOLERANCE, max_iterations=MAX_ITERATIONS):
    x = initial_guess[:]
    n = len(matrix)

    if not is_diagonally_dominant(matrix):
        if rearrange_to_diagonal_dominance(matrix, vector):
            print("Matrix rearranged to achieve diagonal dominance.")
        else:
            print("Matrix is not diagonally dominant and cannot be rearranged.")
            return

    for iteration in range(max_iterations):
        x_new = [0] * n
        for i in range(n):
            s = sum(matrix[i][j] * x[j] for j in range(n) if j != i)
            x_new[i] = (vector[i] - s) / matrix[i][i]

        print(f"Iteration {iteration + 1}: {x_new}")

        if vector_norm_inf(x_new, x) < tol:
            print(f"Results are: {x_new}")
            print(f"Converged in {iteration + 1} iterations.")
            return
        x = x_new

    print("The Jacobi method did not converge.")

# Solve the system using the Gauss-Seidel method
def gauss_seidel_method(matrix, vector, initial_guess, tol=TOLERANCE, max_iterations=MAX_ITERATIONS):
    x = initial_guess[:]
    n = len(matrix)

    if not is_diagonally_dominant(matrix):
        if rearrange_to_diagonal_dominance(matrix, vector):
            print("Matrix rearranged to achieve diagonal dominance.")
        else:
            print("Matrix is not diagonally dominant and cannot be rearranged.")
            return

    for iteration in range(max_iterations):
        x_new = x[:]
        for i in range(n):
            s1 = sum(matrix[i][j] * x_new[j] for j in range(i))
            s2 = sum(matrix[i][j] * x[j] for j in range(i + 1, n))
            x_new[i] = (vector[i] - s1 - s2) / matrix[i][i]

        print(f"Iteration {iteration + 1}: {x_new}")

        if vector_norm_inf(x_new, x) < tol:
            print(f"Results are: {x_new}")
            print(f"Converged in {iteration + 1} iterations.")
            return
        x = x_new

    print("The Gauss-Seidel method did not converge.")

# Main program: Allow user to choose a method or exit
def solve_system(matrixA, vectorB, initial_guess):
    while True:
        method = input("\nChoose a method (jacobi/gauss_seidel) or type 'exit' to quit: ").strip().lower()
        if method == "jacobi":
            jacobi_method(matrixA, vectorB, initial_guess)
        elif method == "gauss_seidel":
            gauss_seidel_method(matrixA, vectorB, initial_guess)
        elif method == "exit":
            print("Exiting the solver. Goodbye!")
            break
        else:
            print("Invalid choice. Please type 'jacobi', 'gauss_seidel', or 'exit'.")

# Run the solver
solve_system(matrixA, vectorB,initial_guess)