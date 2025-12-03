import numpy as np
import matplotlib.pyplot as plt

def generate_chord_function():
    # 1. Define the Data Points from NACA TR 594 (Nose 6-Propeller C)
    # Radius (m)
    radius_data = np.array([0.525, 0.675, 0.825, 0.975, 1.125, 1.275, 1.425])
    # Chord (m)
    chord_data = np.array([0.18, 0.225, 0.225, 0.21, 0.1875, 0.1425, 0.12])

    # 2. Fit a Polynomial (Degree 4 usually fits blade distributions well)
    # We look for coefficients p such that chord ~ p[0]*r^4 + p[1]*r^3 + ...
    degree = 4
    coeffs = np.polyfit(radius_data, chord_data, degree)

    # 3. Create the Lambda Function String for the user
    # We format the coefficients into a Python-executable string
    terms = []
    for i, c in enumerate(coeffs):
        power = degree - i
        if power == 0:
            terms.append(f"{c:.6f}")
        elif power == 1:
            terms.append(f"{c:.6f}*r")
        else:
            terms.append(f"{c:.6f}*r**{power}")
    
    # Join terms with ' + ' and fix any '+ -' issues
    poly_str = " + ".join(terms).replace("+ -", "- ")
    lambda_str = f"chord_lambda = lambda r: {poly_str}"

    print("-" * 30)
    print("Generated Lambda Function:")
    print("-" * 30)
    print(lambda_str)
    print("-" * 30)

    # 4. Define the actual function for plotting
    chord_func = np.poly1d(coeffs)

    # 5. Plotting to Verify the Approximation
    r_smooth = np.linspace(min(radius_data), max(radius_data), 100)
    c_smooth = chord_func(r_smooth)

    plt.figure(figsize=(8, 5))
    plt.plot(radius_data, chord_data, 'ro', label='Original Data (TR 594)')
    plt.plot(r_smooth, c_smooth, 'b-', label=f'Polynomial Approx (Deg {degree})')
    plt.title('Chord Length Distribution Approximation')
    plt.xlabel('Radius (m)')
    plt.ylabel('Chord Length (m)')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.show()

if __name__ == "__main__":
    generate_chord_function()